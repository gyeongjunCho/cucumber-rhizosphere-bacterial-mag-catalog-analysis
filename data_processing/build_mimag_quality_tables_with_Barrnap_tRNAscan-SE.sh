#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################################################################
# Build MIMAG quality tables using Barrnap, tRNAscan-SE, and CheckM2 values
#
# Requirements:
#   - Barrnap
#   - tRNAscan-SE
#   - Python 3
#
# Primary MAG input:
#   8-summary/<sample>/MAGs/*.fa
#
# For example, if pb-metagenomics-tools is installed under:
#   $HOME/projects/pb-metagenomics-tools
#
# then MAG_ROOT would typically be:
#   $HOME/projects/pb-metagenomics-tools/HiFi-MAG-Pipeline/8-summary
#
# MAG metadata:
#   Download 2-02_mag_metadata.tsv from the dataset published with the
#   associated Data Descriptor and provide its path using MAG_METADATA.
#
# Example:
#   MAG_ROOT=path/to/8-summary \
#   MAG_METADATA=path/to/2-02_mag_metadata.tsv \
#   OUT_ROOT=path/to/mimag_quality_output \
#   bash build_mimag_quality_tables_with_barrnap_trnascan.sh
###############################################################################

THREADS="${THREADS:-32}"
EXPECTED_MAGS="${EXPECTED_MAGS:-6505}"
FORCE="${FORCE:-0}"

MAG_ROOT="${MAG_ROOT:-path/to/8-summary}"
MAG_METADATA="${MAG_METADATA:-path/to/2-02_mag_metadata.tsv}"
OUT_ROOT="${OUT_ROOT:-path/to/mimag_quality_output}"

BARRNAP_DIR="$OUT_ROOT/barrnap"
TRNASCAN_DIR="$OUT_ROOT/trnascan"
LOG_DIR="$OUT_ROOT/logs"
DONE_DIR="$OUT_ROOT/done"
TMP_DIR="$OUT_ROOT/tmp"
SUMMARY_DIR="$OUT_ROOT/_summary"

MANIFEST="$SUMMARY_DIR/input_mag_manifest.tsv"
FASTA_LIST="$SUMMARY_DIR/input_fasta.list"
VERSION_TSV="$SUMMARY_DIR/software_versions.tsv"

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

for exe in barrnap tRNAscan-SE python3; do
    command -v "$exe" >/dev/null 2>&1 || die "$exe not found."
done

[[ -d "$MAG_ROOT" ]] || die "MAG_ROOT not found: $MAG_ROOT"

if [[ ! -s "$MAG_METADATA" ]]; then
    cat >&2 <<'EOF'
ERROR: 2-02_mag_metadata.tsv was not found.

Please download 2-02_mag_metadata.tsv from the dataset published with the
associated Data Descriptor, then provide its path with MAG_METADATA.

Example:
  MAG_ROOT=path/to/8-summary \
  MAG_METADATA=path/to/2-02_mag_metadata.tsv \
  OUT_ROOT=path/to/mimag_quality_output \
  bash build_mimag_quality_tables_with_barrnap_trnascan.sh
EOF
    exit 1
fi

MAG_ROOT_ABS="$(readlink -f "$MAG_ROOT")"
MAG_METADATA_ABS="$(readlink -f "$MAG_METADATA")"

mkdir -p \
    "$BARRNAP_DIR" \
    "$TRNASCAN_DIR" \
    "$LOG_DIR" \
    "$DONE_DIR" \
    "$TMP_DIR" \
    "$SUMMARY_DIR"

###############################################################################
# Record software versions
###############################################################################

BARRNAP_VERSION="$(barrnap --version 2>&1 | head -n 1)"
TRNASCAN_VERSION="$(
    tRNAscan-SE -h 2>&1 \
    | grep -m1 '^tRNAscan-SE ' \
    || true
)"

printf "software\tversion\n" > "$VERSION_TSV"
printf "Barrnap\t%s\n" "$BARRNAP_VERSION" >> "$VERSION_TSV"
printf "tRNAscan-SE\t%s\n" "$TRNASCAN_VERSION" >> "$VERSION_TSV"

log "[INFO] MAG_ROOT      = $MAG_ROOT_ABS"
log "[INFO] MAG_METADATA  = $MAG_METADATA_ABS"
log "[INFO] OUT_ROOT      = $(readlink -m "$OUT_ROOT")"
log "[INFO] THREADS       = $THREADS"
log "[INFO] EXPECTED_MAGS = $EXPECTED_MAGS"
log "[INFO] FORCE         = $FORCE"
log "[INFO] $BARRNAP_VERSION"
log "[INFO] $TRNASCAN_VERSION"

###############################################################################
# Discover primary MAG FASTAs
#
# Exact expected layout:
#   8-summary/<sample>/MAGs/*.fa
###############################################################################

log "[INFO] Discovering primary MAG FASTAs from 8-summary/<sample>/MAGs/*.fa ..."

: > "$FASTA_LIST"

for sample_dir in "$MAG_ROOT_ABS"/*; do
    [[ -d "$sample_dir" ]] || continue

    mag_dir="$sample_dir/MAGs"
    [[ -d "$mag_dir" ]] || continue

    for fa in "$mag_dir"/*.fa; do
        [[ -f "$fa" ]] || continue
        printf "%s\n" "$fa" >> "$FASTA_LIST"
    done
done

sort -o "$FASTA_LIST" "$FASTA_LIST"

N_FASTA="$(wc -l < "$FASTA_LIST")"

(( N_FASTA > 0 )) || die "No MAG FASTAs found under $MAG_ROOT_ABS"
(( N_FASTA == EXPECTED_MAGS )) \
    || die "Expected $EXPECTED_MAGS MAG FASTAs, found $N_FASTA"

log "[INFO] Found $N_FASTA MAG FASTAs"

###############################################################################
# Build stable input manifest
#
# sample_id = first directory below MAG_ROOT
# base      = FASTA basename without .fa
# mag_id    = sample_id + "__" + base
###############################################################################

printf "mag_id\tsample_id\tbase\tinput_relpath\n" > "$MANIFEST"

declare -A SEEN_MAG_ID=()

while IFS= read -r fa; do
    [[ -f "$fa" ]] || die "Input FASTA vanished: $fa"

    rel="${fa#"$MAG_ROOT_ABS"/}"
    sample_id="${rel%%/*}"

    filename="$(basename "$fa")"
    base="${filename%.fa}"
    mag_id="${sample_id}__${base}"

    if [[ -n "${SEEN_MAG_ID[$mag_id]:-}" ]]; then
        die "Duplicate derived mag_id: $mag_id"
    fi
    SEEN_MAG_ID["$mag_id"]=1

    printf "%s\t%s\t%s\t%s\n" \
        "$mag_id" "$sample_id" "$base" "$rel" >> "$MANIFEST"
done < "$FASTA_LIST"

N_MANIFEST=$(( $(wc -l < "$MANIFEST") - 1 ))
(( N_MANIFEST == EXPECTED_MAGS )) \
    || die "Manifest contains $N_MANIFEST MAGs, expected $EXPECTED_MAGS"

log "[INFO] Manifest validated: $N_MANIFEST unique MAG IDs"

###############################################################################
# Per-MAG Barrnap + tRNAscan-SE worker
###############################################################################

run_one() {
    local fa="$1"

    local rel sample_id filename base mag_id
    local barrnap_out trna_out log_file done_file
    local tmp_barrnap tmp_trna

    rel="${fa#"$MAG_ROOT_ABS"/}"
    sample_id="${rel%%/*}"
    filename="$(basename "$fa")"
    base="${filename%.fa}"
    mag_id="${sample_id}__${base}"

    barrnap_out="$BARRNAP_DIR/${mag_id}.gff"
    trna_out="$TRNASCAN_DIR/${mag_id}.tsv"
    log_file="$LOG_DIR/${mag_id}.log"
    done_file="$DONE_DIR/${mag_id}.done"

    if [[ "$FORCE" != "1" && -f "$done_file" ]]; then
        return 0
    fi

    tmp_barrnap="$TMP_DIR/${mag_id}.barrnap.$$.tmp"
    tmp_trna="$TMP_DIR/${mag_id}.trnascan.$$.tmp"

    rm -f "$tmp_barrnap" "$tmp_trna"

    {
        echo "[START] $(date '+%F %T')"
        echo "[MAG_ID] $mag_id"
        echo "[INPUT]  $fa"

        echo "[RUN] barrnap --kingdom bac --threads 1"
        barrnap \
            --quiet \
            --kingdom bac \
            --threads 1 \
            "$fa" \
            > "$tmp_barrnap"

        echo "[RUN] tRNAscan-SE -B"
        tRNAscan-SE \
            -B \
            -o "$tmp_trna" \
            "$fa"

        mv -f "$tmp_barrnap" "$barrnap_out"
        mv -f "$tmp_trna" "$trna_out"

        touch "$done_file"

        echo "[DONE] $(date '+%F %T')"
    } > "$log_file" 2>&1
}

export -f run_one
export MAG_ROOT_ABS
export BARRNAP_DIR
export TRNASCAN_DIR
export LOG_DIR
export DONE_DIR
export TMP_DIR
export FORCE

###############################################################################
# Parallel annotation
###############################################################################

log "[INFO] Running Barrnap + tRNAscan-SE on $EXPECTED_MAGS MAGs..."

set +e
xargs -r -P "$THREADS" -n 1 \
    bash -c 'run_one "$1"' _ \
    < "$FASTA_LIST"
XARGS_STATUS=$?
set -e

if (( XARGS_STATUS != 0 )); then
    die "One or more MAG annotation jobs failed. Check $LOG_DIR"
fi

###############################################################################
# Annotation integrity checks
###############################################################################

N_DONE="$(find "$DONE_DIR" -maxdepth 1 -type f -name "*.done" | wc -l)"
N_BARRNAP="$(find "$BARRNAP_DIR" -maxdepth 1 -type f -name "*.gff" | wc -l)"
N_TRNA="$(find "$TRNASCAN_DIR" -maxdepth 1 -type f -name "*.tsv" | wc -l)"

log "[INFO] completed markers = $N_DONE"
log "[INFO] Barrnap GFF files = $N_BARRNAP"
log "[INFO] tRNAscan files    = $N_TRNA"

(( N_DONE == EXPECTED_MAGS )) \
    || die "Only $N_DONE/$EXPECTED_MAGS MAGs completed"

(( N_BARRNAP == EXPECTED_MAGS )) \
    || die "Barrnap outputs: $N_BARRNAP/$EXPECTED_MAGS"

(( N_TRNA == EXPECTED_MAGS )) \
    || die "tRNAscan outputs: $N_TRNA/$EXPECTED_MAGS"

###############################################################################
# Build MIMAG quality tables
#
# Criteria implemented here:
#   - CheckM2 completeness > 90%
#   - CheckM2 contamination < 5%
#   - at least one non-partial 5S, 16S, and 23S rRNA hit
#   - tRNAs for >=18 of 20 canonical amino-acid types
#
# tRNAscan-SE handling:
#   - rows annotated "pseudo" are excluded
#   - Ile2 -> Ile
#   - fMet -> Met
#   - SeC, Sup, and Undet are excluded
#
# Barrnap handling:
#   - hits containing "partial" in GFF attributes are not counted toward HQ
###############################################################################

POSTPROCESS="$OUT_ROOT/summarize_mimag_quality.py"

cat > "$POSTPROCESS" <<'PY'
#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


CANONICAL_AA20 = {
    "Ala", "Arg", "Asn", "Asp", "Cys",
    "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tyr", "Val",
}

TRNA_NORMALIZATION = {
    "Ile2": "Ile",
    "fMet": "Met",
}

EXCLUDED_TRNA_TYPES = {"SeC", "Sup", "Undet"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Classify MAGs against MIMAG high-quality criteria."
    )
    p.add_argument("--mag-metadata", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--barrnap-dir", required=True, type=Path)
    p.add_argument("--trnascan-dir", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--expected-mags", type=int, default=6505)
    return p.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def require_columns(
    rows: list[dict[str, str]],
    required: set[str],
    label: str,
) -> None:
    if not rows:
        raise SystemExit(f"ERROR: {label} is empty.")

    available = set(rows[0].keys())
    missing = required - available

    if missing:
        raise SystemExit(
            f"ERROR: {label} missing required columns: {sorted(missing)}; "
            f"available={sorted(available)}"
        )


def parse_barrnap_gff(path: Path) -> dict:
    total = Counter({"5S": 0, "16S": 0, "23S": 0})
    complete = Counter({"5S": 0, "16S": 0, "23S": 0})
    partial = Counter({"5S": 0, "16S": 0, "23S": 0})
    unexpected_names = Counter()

    name_re = re.compile(r"(?:^|;)Name=([^;]+)")

    with path.open(encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 9:
                raise ValueError(f"{path}:{lineno}: malformed GFF line")

            attrs = fields[8]
            match = name_re.search(attrs)

            if not match:
                continue

            name = match.group(1)

            rrna_type = {
                "5S_rRNA": "5S",
                "16S_rRNA": "16S",
                "23S_rRNA": "23S",
            }.get(name)

            if rrna_type is None:
                unexpected_names[name] += 1
                continue

            total[rrna_type] += 1

            if "partial" in attrs.lower():
                partial[rrna_type] += 1
            else:
                complete[rrna_type] += 1

    return {
        "total": total,
        "complete": complete,
        "partial": partial,
        "unexpected_names": unexpected_names,
    }


def parse_trnascan_table(path: Path) -> dict:
    canonical_counts = Counter()
    raw_type_counts = Counter()
    excluded_type_counts = Counter()
    unknown_type_counts = Counter()

    total_rows = 0
    nonpseudo_rows = 0
    pseudo_rows = 0

    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            raw = line.strip()

            if not raw:
                continue

            fields = raw.split()

            if len(fields) < 9:
                continue

            try:
                int(fields[1])
                int(fields[2])
                int(fields[3])
                float(fields[8])
            except (ValueError, IndexError):
                continue

            total_rows += 1

            raw_type = fields[4]
            raw_type_counts[raw_type] += 1

            if "pseudo" in raw.lower():
                pseudo_rows += 1
                continue

            nonpseudo_rows += 1

            if raw_type in EXCLUDED_TRNA_TYPES:
                excluded_type_counts[raw_type] += 1
                continue

            normalized = TRNA_NORMALIZATION.get(raw_type, raw_type)

            if normalized in CANONICAL_AA20:
                canonical_counts[normalized] += 1
            else:
                unknown_type_counts[raw_type] += 1

    canonical_types = sorted(canonical_counts)

    return {
        "total_rows": total_rows,
        "nonpseudo_rows": nonpseudo_rows,
        "pseudo_rows": pseudo_rows,
        "canonical_counts": canonical_counts,
        "canonical_types": canonical_types,
        "raw_type_counts": raw_type_counts,
        "excluded_type_counts": excluded_type_counts,
        "unknown_type_counts": unknown_type_counts,
    }


def counter_string(counter: Counter) -> str:
    if not counter:
        return ""
    return ";".join(f"{k}:{counter[k]}" for k in sorted(counter))


def main() -> None:
    args = parse_args()

    for path, label in [
        (args.mag_metadata, "MAG metadata"),
        (args.manifest, "annotation manifest"),
    ]:
        if not path.is_file():
            raise SystemExit(f"ERROR: {label} not found: {path}")

    if not args.barrnap_dir.is_dir():
        raise SystemExit(
            f"ERROR: Barrnap directory not found: {args.barrnap_dir}"
        )

    if not args.trnascan_dir.is_dir():
        raise SystemExit(
            f"ERROR: tRNAscan-SE directory not found: {args.trnascan_dir}"
        )

    args.outdir.mkdir(parents=True, exist_ok=True)

    metadata_rows = read_tsv(args.mag_metadata)
    require_columns(
        metadata_rows,
        {"mag_id", "sample_id", "completeness", "contamination"},
        "MAG metadata",
    )

    manifest_rows = read_tsv(args.manifest)
    require_columns(
        manifest_rows,
        {"mag_id", "sample_id", "base", "input_relpath"},
        "annotation manifest",
    )

    if args.expected_mags > 0:
        if len(metadata_rows) != args.expected_mags:
            raise SystemExit(
                f"ERROR: MAG metadata has {len(metadata_rows)} rows; "
                f"expected {args.expected_mags}."
            )

        if len(manifest_rows) != args.expected_mags:
            raise SystemExit(
                f"ERROR: manifest has {len(manifest_rows)} rows; "
                f"expected {args.expected_mags}."
            )

    metadata = {}
    metadata_order = []

    for row in metadata_rows:
        mag_id = row["mag_id"].strip()

        if not mag_id:
            raise SystemExit("ERROR: empty mag_id in MAG metadata")

        if mag_id in metadata:
            raise SystemExit(
                f"ERROR: duplicate mag_id in MAG metadata: {mag_id}"
            )

        try:
            completeness = float(row["completeness"])
            contamination = float(row["contamination"])
        except ValueError:
            raise SystemExit(
                f"ERROR: invalid completeness/contamination for {mag_id}"
            )

        metadata[mag_id] = {
            "sample_id": row["sample_id"].strip(),
            "completeness": completeness,
            "contamination": contamination,
        }
        metadata_order.append(mag_id)

    manifest_by_base = {}

    for row in manifest_rows:
        base = row["base"].strip()
        internal_id = row["mag_id"].strip()
        sample_id = row["sample_id"].strip()

        if base in manifest_by_base:
            raise SystemExit(f"ERROR: duplicate manifest base: {base}")

        if base not in metadata:
            raise SystemExit(
                f"ERROR: manifest base {base!r} not present in MAG metadata"
            )

        if sample_id != metadata[base]["sample_id"]:
            raise SystemExit(
                f"ERROR: sample mismatch for {base}: "
                f"manifest={sample_id}, metadata={metadata[base]['sample_id']}"
            )

        manifest_by_base[base] = {
            "internal_id": internal_id,
            "sample_id": sample_id,
        }

    missing_manifest = [
        mag_id
        for mag_id in metadata_order
        if mag_id not in manifest_by_base
    ]

    if missing_manifest:
        raise SystemExit(
            f"ERROR: {len(missing_manifest)} metadata MAGs missing from manifest; "
            f"examples={missing_manifest[:10]}"
        )

    results = []

    global_raw_trna = Counter()
    global_excluded_trna = Counter()
    global_unknown_trna = Counter()
    global_unexpected_barrnap = Counter()
    total_partial_rrna = Counter()

    for mag_id in metadata_order:
        md = metadata[mag_id]
        ann = manifest_by_base[mag_id]
        internal_id = ann["internal_id"]

        barrnap_path = args.barrnap_dir / f"{internal_id}.gff"
        trna_path = args.trnascan_dir / f"{internal_id}.tsv"

        if not barrnap_path.is_file():
            raise SystemExit(
                f"ERROR: missing Barrnap file: {barrnap_path}"
            )

        if not trna_path.is_file():
            raise SystemExit(
                f"ERROR: missing tRNAscan-SE file: {trna_path}"
            )

        rrna = parse_barrnap_gff(barrnap_path)
        trna = parse_trnascan_table(trna_path)

        global_raw_trna.update(trna["raw_type_counts"])
        global_excluded_trna.update(trna["excluded_type_counts"])
        global_unknown_trna.update(trna["unknown_type_counts"])
        global_unexpected_barrnap.update(rrna["unexpected_names"])
        total_partial_rrna.update(rrna["partial"])

        completeness = md["completeness"]
        contamination = md["contamination"]

        pass_completeness = completeness > 90.0
        pass_contamination = contamination < 5.0

        pass_5s = rrna["complete"]["5S"] >= 1
        pass_16s = rrna["complete"]["16S"] >= 1
        pass_23s = rrna["complete"]["23S"] >= 1
        pass_all_rrna = pass_5s and pass_16s and pass_23s

        trna_aa_types = len(trna["canonical_types"])
        pass_trna = trna_aa_types >= 18

        mimag_hq = (
            pass_completeness
            and pass_contamination
            and pass_all_rrna
            and pass_trna
        )

        failures = []

        if not pass_completeness:
            failures.append("completeness<=90")
        if not pass_contamination:
            failures.append("contamination>=5")
        if not pass_5s:
            failures.append("5S_rRNA_missing_or_partial")
        if not pass_16s:
            failures.append("16S_rRNA_missing_or_partial")
        if not pass_23s:
            failures.append("23S_rRNA_missing_or_partial")
        if not pass_trna:
            failures.append("tRNA_aa_types<18")

        results.append({
            "mag_id": mag_id,
            "sample_id": md["sample_id"],
            "completeness": f"{completeness:.6g}",
            "contamination": f"{contamination:.6g}",

            "rrna_5S_hits": rrna["total"]["5S"],
            "rrna_5S_complete_hits": rrna["complete"]["5S"],
            "rrna_5S_partial_hits": rrna["partial"]["5S"],

            "rrna_16S_hits": rrna["total"]["16S"],
            "rrna_16S_complete_hits": rrna["complete"]["16S"],
            "rrna_16S_partial_hits": rrna["partial"]["16S"],

            "rrna_23S_hits": rrna["total"]["23S"],
            "rrna_23S_complete_hits": rrna["complete"]["23S"],
            "rrna_23S_partial_hits": rrna["partial"]["23S"],

            "trna_total_predictions": trna["total_rows"],
            "trna_nonpseudo_predictions": trna["nonpseudo_rows"],
            "trna_pseudogene_predictions": trna["pseudo_rows"],
            "trna_aa_types": trna_aa_types,
            "trna_aa_type_names": ",".join(trna["canonical_types"]),
            "trna_excluded_special_types": counter_string(
                trna["excluded_type_counts"]
            ),
            "trna_unknown_types": counter_string(
                trna["unknown_type_counts"]
            ),

            "pass_completeness_gt90": int(pass_completeness),
            "pass_contamination_lt5": int(pass_contamination),
            "pass_5S_rRNA": int(pass_5s),
            "pass_16S_rRNA": int(pass_16s),
            "pass_23S_rRNA": int(pass_23s),
            "pass_all_5S_16S_23S_rRNA": int(pass_all_rrna),
            "pass_trna_at_least_18_of_20_aa_types": int(pass_trna),

            "mimag_high_quality": int(mimag_hq),
            "hq_failure_reasons": (
                "PASS" if mimag_hq else ";".join(failures)
            ),
        })

    if global_unknown_trna:
        raise SystemExit(
            "ERROR: unexpected noncanonical tRNAscan-SE Type values were found: "
            + counter_string(global_unknown_trna)
        )

    detailed_fields = list(results[0].keys())

    quality_path = args.outdir / "mimag_quality.tsv"
    write_tsv(
        quality_path,
        results,
        detailed_fields,
    )

    total = len(results)

    hq = sum(
        int(row["mimag_high_quality"])
        for row in results
    )

    checkm_hq = sum(
        int(row["pass_completeness_gt90"])
        and int(row["pass_contamination_lt5"])
        for row in results
    )

    all_rrna = sum(
        int(row["pass_all_5S_16S_23S_rRNA"])
        for row in results
    )

    trna18 = sum(
        int(row["pass_trna_at_least_18_of_20_aa_types"])
        for row in results
    )

    summary_rows = [
        {
            "metric": "total_MAGs",
            "value": total,
        },
        {
            "metric": "CheckM2_completeness_gt90_and_contamination_lt5",
            "value": checkm_hq,
        },
        {
            "metric": "all_5S_16S_23S_rRNA_present",
            "value": all_rrna,
        },
        {
            "metric": "tRNA_at_least_18_of_20_amino_acid_types",
            "value": trna18,
        },
        {
            "metric": "MIMAG_high_quality_MAGs",
            "value": hq,
        },
        {
            "metric": "MIMAG_high_quality_percent",
            "value": f"{100.0 * hq / total:.3f}",
        },
    ]

    summary_path = args.outdir / "mimag_summary.tsv"
    write_tsv(
        summary_path,
        summary_rows,
        ["metric", "value"],
    )

    by_sample = defaultdict(list)

    for row in results:
        by_sample[row["sample_id"]].append(row)

    sample_rows = []

    for sample_id in sorted(by_sample):
        sample_results = by_sample[sample_id]
        n = len(sample_results)

        sample_hq = sum(
            int(row["mimag_high_quality"])
            for row in sample_results
        )

        sample_rows.append({
            "sample_id": sample_id,
            "MAGs": n,
            "MIMAG_high_quality_MAGs": sample_hq,
            "MIMAG_high_quality_percent": (
                f"{100.0 * sample_hq / n:.3f}"
            ),
            "CheckM2_gt90_lt5": sum(
                int(row["pass_completeness_gt90"])
                and int(row["pass_contamination_lt5"])
                for row in sample_results
            ),
            "all_5S_16S_23S_rRNA_present": sum(
                int(row["pass_all_5S_16S_23S_rRNA"])
                for row in sample_results
            ),
            "tRNA_at_least_18_of_20_aa_types": sum(
                int(row["pass_trna_at_least_18_of_20_aa_types"])
                for row in sample_results
            ),
        })

    sample_path = args.outdir / "mimag_by_sample.tsv"
    write_tsv(
        sample_path,
        sample_rows,
        list(sample_rows[0].keys()),
    )

    qc_rows = [
        {
            "metric": "mag_metadata_rows",
            "value": len(metadata_rows),
        },
        {
            "metric": "manifest_rows",
            "value": len(manifest_rows),
        },
        {
            "metric": "classified_MAGs",
            "value": len(results),
        },
        {
            "metric": "Barrnap_partial_5S_hits",
            "value": total_partial_rrna["5S"],
        },
        {
            "metric": "Barrnap_partial_16S_hits",
            "value": total_partial_rrna["16S"],
        },
        {
            "metric": "Barrnap_partial_23S_hits",
            "value": total_partial_rrna["23S"],
        },
        {
            "metric": "Barrnap_unexpected_Name_values",
            "value": counter_string(global_unexpected_barrnap),
        },
        {
            "metric": "tRNAscan_raw_Type_counts",
            "value": counter_string(global_raw_trna),
        },
        {
            "metric": "tRNAscan_excluded_noncanonical_Type_counts_nonpseudo",
            "value": counter_string(global_excluded_trna),
        },
        {
            "metric": "tRNAscan_unknown_Type_counts_nonpseudo",
            "value": counter_string(global_unknown_trna),
        },
        {
            "metric": "tRNA_normalization_Ile2",
            "value": "Ile",
        },
        {
            "metric": "tRNA_normalization_fMet",
            "value": "Met",
        },
        {
            "metric": "tRNA_excluded_types",
            "value": "SeC,Sup,Undet",
        },
        {
            "metric": "tRNA_pseudogenes_counted_for_18_of_20",
            "value": "no",
        },
        {
            "metric": "Barrnap_partial_hits_counted_for_HQ",
            "value": "no",
        },
        {
            "metric": "MIMAG_completeness_rule",
            "value": ">90%",
        },
        {
            "metric": "MIMAG_contamination_rule",
            "value": "<5%",
        },
        {
            "metric": "MIMAG_rRNA_rule",
            "value": "5S>=1;16S>=1;23S>=1",
        },
        {
            "metric": "MIMAG_tRNA_rule",
            "value": ">=18_of_20_canonical_amino_acid_types",
        },
    ]

    qc_path = args.outdir / "mimag_qc.tsv"
    write_tsv(
        qc_path,
        qc_rows,
        ["metric", "value"],
    )

    print(f"[DONE] {quality_path}")
    print(f"[DONE] {summary_path}")
    print(f"[DONE] {sample_path}")
    print(f"[DONE] {qc_path}")
    print()
    print("metric\tvalue")

    for row in summary_rows:
        print(f"{row['metric']}\t{row['value']}")


if __name__ == "__main__":
    main()
PY

python3 "$POSTPROCESS" \
    --mag-metadata "$MAG_METADATA_ABS" \
    --manifest "$MANIFEST" \
    --barrnap-dir "$BARRNAP_DIR" \
    --trnascan-dir "$TRNASCAN_DIR" \
    --outdir "$SUMMARY_DIR" \
    --expected-mags "$EXPECTED_MAGS"

###############################################################################
# Compact run summary
###############################################################################

RUN_SUMMARY="$SUMMARY_DIR/run_summary.tsv"

printf "metric\tvalue\n" > "$RUN_SUMMARY"
printf "input_MAGs\t%s\n" "$N_FASTA" >> "$RUN_SUMMARY"
printf "completed_MAGs\t%s\n" "$N_DONE" >> "$RUN_SUMMARY"
printf "barrnap_outputs\t%s\n" "$N_BARRNAP" >> "$RUN_SUMMARY"
printf "trnascan_outputs\t%s\n" "$N_TRNA" >> "$RUN_SUMMARY"
printf "threads\t%s\n" "$THREADS" >> "$RUN_SUMMARY"

log "[DONE] MIMAG RNA annotation and quality tables completed"
log "[OUTPUT] $(readlink -m "$OUT_ROOT")"
log "[SUMMARY] $(readlink -m "$SUMMARY_DIR")"

echo
cat "$VERSION_TSV"
echo
cat "$RUN_SUMMARY"
echo
cat "$SUMMARY_DIR/mimag_summary.tsv"

