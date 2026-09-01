#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

THREADS="${THREADS:-32}"
EXPECTED_MAGS="${EXPECTED_MAGS:-6505}"
ANI="${ANI:-0.95}"
P_ANI="${P_ANI:-0.90}"
COV="${COV:-0.50}"
MIN_COMP="${MIN_COMP:-50}"
MAX_CON="${MAX_CON:-10}"

# Primary MAG FASTA files from the HiFi-MAG-Pipeline summary output:
#   8-summary/<sample>/MAGs/*.fa
#
# For example, if pb-metagenomics-tools is installed under
#   $HOME/projects/pb-metagenomics-tools
# then MAG_ROOT would typically be:
#   $HOME/projects/pb-metagenomics-tools/HiFi-MAG-Pipeline/8-summary
#
# Replace the placeholder below or export MAG_ROOT before running.
MAG_ROOT="${MAG_ROOT:-path/to/8-summary}"

# Choose any output directory with sufficient space/inodes.
# Replace the placeholder below or export OUT_ROOT before running.
OUT_ROOT="${OUT_ROOT:-path/to/ANI95_output}"

# MAG metadata released with the associated Data Descriptor.
# Download 2-02_mag_metadata.tsv from the published dataset and provide its path.
MAG_METADATA="${MAG_METADATA:-path/to/2-02_mag_metadata.tsv}"

GENOME_DIR="$OUT_ROOT/input_genomes"
DREP_OUT="$OUT_ROOT/drep"
SUMMARY_DIR="$OUT_ROOT/_summary"
REP_DIR="$OUT_ROOT/representative_MAGs_ANI95"
LOG_DIR="$OUT_ROOT/logs"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

for exe in dRep fastANI mash python3; do
    command -v "$exe" >/dev/null 2>&1 || die "$exe not found."
done

[[ -d "$MAG_ROOT" ]] || die "MAG_ROOT not found: $MAG_ROOT"

if ! dRep dereplicate -h 2>&1 | grep -q -- "--genomeInfo"; then
    die "This dRep installation does not expose --genomeInfo."
fi

mkdir -p "$GENOME_DIR" "$SUMMARY_DIR" "$REP_DIR" "$LOG_DIR"

if [[ ! -s "$MAG_METADATA" ]]; then
    cat >&2 <<'EOF'
ERROR: 2-02_mag_metadata.tsv was not found.

Please download 2-02_mag_metadata.tsv from the dataset published with the
associated Data Descriptor, then provide its path with MAG_METADATA.

Example:
  MAG_ROOT=path/to/8-summary \
  MAG_METADATA=path/to/2-02_mag_metadata.tsv \
  OUT_ROOT=path/to/ANI95_output \
  bash run_drep_MAGs_ANI95_catalogue.sh
EOF
    exit 1
fi

if [[ -e "$DREP_OUT/data_tables/Cdb.csv" ]]; then
    die "Existing dRep output detected: $DREP_OUT ; move/remove it or change OUT_ROOT."
fi

log "[INFO] dRep version  = $(dRep --version 2>&1 | head -1 || true)"
log "[INFO] MAG_ROOT      = $MAG_ROOT"
log "[INFO] MAG_METADATA  = $MAG_METADATA"
log "[INFO] OUT_ROOT      = $OUT_ROOT"
log "[INFO] THREADS       = $THREADS"
log "[INFO] ANI/P_ANI/COV = $ANI / $P_ANI / $COV"
log "[INFO] dRep QC filter= completeness >= $MIN_COMP ; contamination <= $MAX_CON"

# ----------------------------------------------------------------------
# Collect the 6,505 primary MAG nucleotide FASTA files directly from
# 8-summary/<sample>/MAGs/*.fa and create stable dRep input names.
# ----------------------------------------------------------------------

MAP_TSV="$SUMMARY_DIR/input_MAG_map.tsv"
DREP_GENOME_LIST="$SUMMARY_DIR/drep_genomes.list"

printf "MAG_ID\tsample_id\tbase\tsource_fa\tdrep_fna\n" > "$MAP_TSV"
: > "$DREP_GENOME_LIST"

count=0

for sample_dir in "$MAG_ROOT"/*; do
    [[ -d "$sample_dir" ]] || continue

    sample_id="$(basename "$sample_dir")"
    mag_dir="$sample_dir/MAGs"
    [[ -d "$mag_dir" ]] || continue

    for fa in "$mag_dir"/*.fa; do
        [[ -f "$fa" ]] || continue

        base="$(basename "$fa" .fa)"
        mag_id="${sample_id}__${base}"
        target="$GENOME_DIR/${mag_id}.fna"

        src="$(readlink -f "$fa")"
        [[ -n "$src" && -f "$src" ]] || die "Source FASTA not found: $fa"

        ln -sfn "$src" "$target"

        printf "%s\t%s\t%s\t%s\t%s\n" \
            "$mag_id" "$sample_id" "$base" "$src" "$target" >> "$MAP_TSV"
        printf "%s\n" "$target" >> "$DREP_GENOME_LIST"

        count=$((count + 1))
    done
done

(( count == EXPECTED_MAGS )) || die "Expected $EXPECTED_MAGS MAG FASTAs, found $count"

broken_links="$(find "$GENOME_DIR" -xtype l | wc -l)"
(( broken_links == 0 )) || die "Broken genome symlinks detected: $broken_links"

bad_paths=0
while IFS= read -r genome_path; do
    if [[ ! -f "$genome_path" ]]; then
        echo "[BAD_PATH] $genome_path" >&2
        bad_paths=$((bad_paths + 1))
    fi
done < "$DREP_GENOME_LIST"
(( bad_paths == 0 )) || die "Genome list contains $bad_paths invalid paths"

list_count="$(wc -l < "$DREP_GENOME_LIST")"
(( list_count == EXPECTED_MAGS )) || die "Genome list has $list_count paths; expected $EXPECTED_MAGS"

dup_count="$(tail -n +2 "$MAP_TSV" | cut -f1 | sort | uniq -d | wc -l)"
(( dup_count == 0 )) || die "Duplicate MAG_IDs detected."

log "[INFO] Genome path validation passed: $list_count/$EXPECTED_MAGS valid FASTAs"

# ----------------------------------------------------------------------
# Build dRep --genomeInfo from CheckM2 completeness and contamination
# values in the released 2-02_mag_metadata.tsv metadata table.
# ----------------------------------------------------------------------

GENOME_INFO="$SUMMARY_DIR/drep_genomeInfo.csv"
UNMATCHED="$SUMMARY_DIR/unmatched_MAG_metadata.tsv"

MAP_TSV="$MAP_TSV" \
MAG_METADATA="$MAG_METADATA" \
GENOME_INFO="$GENOME_INFO" \
UNMATCHED="$UNMATCHED" \
python3 - <<'PY'
from pathlib import Path
import csv
import os

map_tsv = Path(os.environ["MAP_TSV"])
meta_tsv = Path(os.environ["MAG_METADATA"])
out_csv = Path(os.environ["GENOME_INFO"])
unmatched_tsv = Path(os.environ["UNMATCHED"])

with meta_tsv.open(newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    required = {"mag_id", "completeness", "contamination"}
    missing = required - set(r.fieldnames or [])
    if missing:
        raise SystemExit(f"Missing metadata columns: {sorted(missing)}")
    meta = list(r)

by_id = {}
for x in meta:
    if x["mag_id"] in by_id:
        raise SystemExit(f'Duplicate mag_id: {x["mag_id"]}')
    by_id[x["mag_id"]] = x

by_base, dup = {}, set()
for mid, x in by_id.items():
    base = mid.split("__", 1)[-1]
    if base in by_base:
        dup.add(base)
    else:
        by_base[base] = x

for base in dup:
    by_base.pop(base, None)

rows, unmatched = [], []

with map_tsv.open(newline="") as f:
    for x in csv.DictReader(f, delimiter="\t"):
        hit = by_id.get(x["MAG_ID"]) or by_id.get(x["base"]) or by_base.get(x["base"])

        if hit is None:
            unmatched.append(x)
            continue

        comp = float(hit["completeness"])
        con = float(hit["contamination"])

        if not (0 <= comp <= 100 and 0 <= con <= 100):
            raise SystemExit(f'Out-of-range quality for {x["MAG_ID"]}')

        rows.append({
            "genome": Path(x["drep_fna"]).name,
            "completeness": f"{comp:.6f}",
            "contamination": f"{con:.6f}",
            "MAG_ID": x["MAG_ID"],
            "metadata_mag_id": hit["mag_id"],
        })

if unmatched:
    with unmatched_tsv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=unmatched[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(unmatched)
    raise SystemExit(f"{len(unmatched)} FASTAs unmatched; see {unmatched_tsv}")

if len(rows) != len(meta):
    raise SystemExit(
        f"FASTA/metadata mismatch: {len(rows)} matched vs {len(meta)} metadata rows"
    )

with out_csv.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=["genome", "completeness", "contamination"],
    )
    w.writeheader()
    for x in rows:
        w.writerow({
            "genome": x["genome"],
            "completeness": x["completeness"],
            "contamination": x["contamination"],
        })

audit = out_csv.with_name("drep_genomeInfo_match_audit.tsv")
with audit.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=[
            "genome",
            "MAG_ID",
            "metadata_mag_id",
            "completeness",
            "contamination",
        ],
        delimiter="\t",
    )
    w.writeheader()
    w.writerows(rows)

print(f"[INFO] genomeInfo rows: {len(rows)}")
print(f"[INFO] genomeInfo: {out_csv}")
print(f"[INFO] audit: {audit}")
PY

n_info=$(( $(wc -l < "$GENOME_INFO") - 1 ))
(( n_info == EXPECTED_MAGS )) || die "genomeInfo has $n_info rows; expected $EXPECTED_MAGS"

# ----------------------------------------------------------------------
# dRep ANI95 dereplication.
# ----------------------------------------------------------------------

log "[INFO] Running dRep ANI95 with existing CheckM2 quality estimates."

dRep dereplicate "$DREP_OUT" \
    -g "$DREP_GENOME_LIST" \
    -p "$THREADS" \
    -comp "$MIN_COMP" \
    -con "$MAX_CON" \
    --genomeInfo "$GENOME_INFO" \
    -pa "$P_ANI" \
    -sa "$ANI" \
    -nc "$COV" \
    --S_algorithm fastANI

# ----------------------------------------------------------------------
# Convert dRep cluster IDs to stable public ANI95 IDs and generate
# membership, cluster-summary, representative-metadata, presence/absence,
# and sample-level cluster-count tables.
# ----------------------------------------------------------------------

POSTPROCESS="$OUT_ROOT/postprocess_drep_ani95.py"

cat > "$POSTPROCESS" <<'PY'
#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import os
import re

DREP_OUT = Path(os.environ["DREP_OUT"])
SUMMARY_DIR = Path(os.environ["SUMMARY_DIR"])
MAP_TSV = Path(os.environ["MAP_TSV"])
MAG_METADATA = Path(os.environ["MAG_METADATA"])
REP_DIR = Path(os.environ["REP_DIR"])

cdb = DREP_OUT / "data_tables" / "Cdb.csv"
wdb = DREP_OUT / "data_tables" / "Wdb.csv"

for p in (cdb, wdb):
    if not p.exists():
        raise SystemExit(f"Missing dRep table: {p}")

mag_map = {}
genome_to_mag = {}

with MAP_TSV.open(newline="") as f:
    for x in csv.DictReader(f, delimiter="\t"):
        mag_map[x["MAG_ID"]] = x
        genome_to_mag[Path(x["drep_fna"]).name] = x["MAG_ID"]

def natural_key(x):
    return [int(v) if v.isdigit() else v for v in re.split(r"(\d+)", str(x))]

def find_col(fields, preferred):
    for name in preferred:
        if name in fields:
            return name
    return None

cluster_members = defaultdict(list)

with cdb.open(newline="") as f:
    r = csv.DictReader(f)
    fields = r.fieldnames or []

    gcol = find_col(fields, ["genome", "genome_name", "genome_id"])
    ccol = find_col(
        fields,
        ["secondary_cluster", "secondary_cluster_id", "cluster", "cluster_id"],
    )

    if ccol is None:
        ccol = next(
            (
                c for c in fields
                if "secondary" in c.lower() and "cluster" in c.lower()
            ),
            None,
        )

    if gcol is None or ccol is None:
        raise SystemExit(f"Cannot parse Cdb columns: {fields}")

    for x in r:
        gname = Path(x[gcol]).name
        mag_id = genome_to_mag.get(gname) or Path(gname).stem

        if mag_id not in mag_map:
            raise SystemExit(f"Cdb genome not in input map: {gname}")

        cluster_members[str(x[ccol])].append(mag_id)

winner_by_cluster = {}
score_by_cluster = {}

with wdb.open(newline="") as f:
    r = csv.DictReader(f)
    fields = r.fieldnames or []

    gcol = find_col(fields, ["genome"])
    ccol = find_col(fields, ["cluster", "secondary_cluster"])
    scol = find_col(fields, ["score"])

    if gcol is None or ccol is None:
        raise SystemExit(f"Cannot parse Wdb columns: {fields}")

    for x in r:
        gname = Path(x[gcol]).name
        mag_id = genome_to_mag.get(gname) or Path(gname).stem

        if mag_id not in mag_map:
            raise SystemExit(f"Wdb genome not in input map: {gname}")

        cluster = str(x[ccol])
        winner_by_cluster[cluster] = mag_id
        score_by_cluster[cluster] = x.get(scol, "") if scol else ""

missing = set(cluster_members) - set(winner_by_cluster)
if missing:
    raise SystemExit(f"Clusters without winner: {sorted(missing)[:10]}")

clusters = sorted(cluster_members, key=natural_key)
public_id = {
    cluster: f"ANI95_{i:05d}"
    for i, cluster in enumerate(clusters, 1)
}

cluster_samples = {
    cluster: sorted({
        mag_map[mag_id]["sample_id"]
        for mag_id in members
    })
    for cluster, members in cluster_members.items()
}

membership = SUMMARY_DIR / "ani95_cluster_membership.tsv"
membership_fields = [
    "mag_id",
    "sample_id",
    "ani95_cluster_id",
    "drep_secondary_cluster",
    "representative_mag_id",
    "is_representative",
    "drep_winner_score",
    "cluster_size",
    "n_samples",
]

with membership.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=membership_fields, delimiter="\t")
    w.writeheader()

    for cluster in clusters:
        rep = winner_by_cluster[cluster]

        for mag_id in sorted(cluster_members[cluster]):
            w.writerow({
                "mag_id": mag_id,
                "sample_id": mag_map[mag_id]["sample_id"],
                "ani95_cluster_id": public_id[cluster],
                "drep_secondary_cluster": cluster,
                "representative_mag_id": rep,
                "is_representative": 1 if mag_id == rep else 0,
                "drep_winner_score": (
                    score_by_cluster[cluster] if mag_id == rep else ""
                ),
                "cluster_size": len(cluster_members[cluster]),
                "n_samples": len(cluster_samples[cluster]),
            })

summary = SUMMARY_DIR / "ani95_cluster_summary.tsv"
summary_fields = [
    "ani95_cluster_id",
    "drep_secondary_cluster",
    "representative_mag_id",
    "drep_winner_score",
    "cluster_size",
    "n_samples",
    "member_mag_ids",
    "member_samples",
]

with summary.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=summary_fields, delimiter="\t")
    w.writeheader()

    for cluster in clusters:
        w.writerow({
            "ani95_cluster_id": public_id[cluster],
            "drep_secondary_cluster": cluster,
            "representative_mag_id": winner_by_cluster[cluster],
            "drep_winner_score": score_by_cluster[cluster],
            "cluster_size": len(cluster_members[cluster]),
            "n_samples": len(cluster_samples[cluster]),
            "member_mag_ids": ";".join(sorted(cluster_members[cluster])),
            "member_samples": ";".join(cluster_samples[cluster]),
        })

with MAG_METADATA.open(newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    metadata_fields = r.fieldnames or []
    meta_rows = list(r)

by_id = {x["mag_id"]: x for x in meta_rows}
by_base = {}
duplicate_bases = set()

for x in meta_rows:
    base = x["mag_id"].split("__", 1)[-1]

    if base in by_base:
        duplicate_bases.add(base)
    else:
        by_base[base] = x

for base in duplicate_bases:
    by_base.pop(base, None)

def metadata_for(mag_id):
    base = mag_map[mag_id]["base"]
    return (
        by_id.get(mag_id)
        or by_id.get(base)
        or by_base.get(base)
        or {}
    )

repmeta = SUMMARY_DIR / "ani95_representative_metadata.tsv"
prefix_fields = [
    "ani95_cluster_id",
    "drep_secondary_cluster",
    "representative_mag_id",
    "drep_winner_score",
    "cluster_size",
    "n_samples",
]

with repmeta.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        fieldnames=prefix_fields + metadata_fields,
        delimiter="\t",
    )
    w.writeheader()

    for cluster in clusters:
        rep = winner_by_cluster[cluster]

        row = {
            "ani95_cluster_id": public_id[cluster],
            "drep_secondary_cluster": cluster,
            "representative_mag_id": rep,
            "drep_winner_score": score_by_cluster[cluster],
            "cluster_size": len(cluster_members[cluster]),
            "n_samples": len(cluster_samples[cluster]),
        }

        row.update(metadata_for(rep))
        w.writerow(row)

REP_DIR.mkdir(parents=True, exist_ok=True)

for p in REP_DIR.glob("*.fna"):
    p.unlink()

for cluster in clusters:
    rep = winner_by_cluster[cluster]
    src = Path(mag_map[rep]["drep_fna"]).resolve()
    dst = REP_DIR / f"{public_id[cluster]}.fna"
    dst.symlink_to(src)

samples = sorted({
    x["sample_id"]
    for x in mag_map.values()
})

presence = SUMMARY_DIR / "ani95_cluster_by_sample_presence.tsv"

with presence.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ani95_cluster_id", "representative_mag_id"] + samples)

    for cluster in clusters:
        present_samples = set(cluster_samples[cluster])

        w.writerow(
            [
                public_id[cluster],
                winner_by_cluster[cluster],
            ]
            + [
                1 if sample in present_samples else 0
                for sample in samples
            ]
        )

sample_clusters = defaultdict(set)

for cluster, members in cluster_members.items():
    for mag_id in members:
        sample = mag_map[mag_id]["sample_id"]
        sample_clusters[sample].add(public_id[cluster])

counts = SUMMARY_DIR / "ani95_sample_cluster_counts.tsv"

with counts.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["sample_id", "n_ani95_clusters"])

    for sample in samples:
        w.writerow([sample, len(sample_clusters[sample])])

print(f"[DONE] membership: {membership}")
print(f"[DONE] cluster summary: {summary}")
print(f"[DONE] representative metadata: {repmeta}")
print(f"[DONE] representative FASTAs: {REP_DIR}")
print(f"[DONE] presence matrix: {presence}")
print(f"[DONE] sample cluster counts: {counts}")
print(f"[INFO] n_MAGs={sum(map(len, cluster_members.values()))}")
print(f"[INFO] n_clusters={len(clusters)}")
print(f"[INFO] n_representatives={len(winner_by_cluster)}")
PY

DREP_OUT="$DREP_OUT" \
SUMMARY_DIR="$SUMMARY_DIR" \
MAP_TSV="$MAP_TSV" \
MAG_METADATA="$MAG_METADATA" \
REP_DIR="$REP_DIR" \
python3 "$POSTPROCESS" | tee "$LOG_DIR/postprocess_drep_ani95.log"

membership_n=$(( $(wc -l < "$SUMMARY_DIR/ani95_cluster_membership.tsv") - 1 ))
(( membership_n == EXPECTED_MAGS )) || \
    die "Membership has $membership_n rows; expected $EXPECTED_MAGS"

n_clusters=$(( $(wc -l < "$SUMMARY_DIR/ani95_cluster_summary.tsv") - 1 ))
n_reps=$(( $(wc -l < "$SUMMARY_DIR/ani95_representative_metadata.tsv") - 1 ))

(( n_clusters == n_reps )) || \
    die "Cluster/representative mismatch: $n_clusters vs $n_reps"

log "[DONE] dRep ANI95 catalogue pipeline finished"
log "[INFO] clusters/representatives = $n_clusters"
log "[SUMMARY] $SUMMARY_DIR"
log "[REP FASTA] $REP_DIR"
