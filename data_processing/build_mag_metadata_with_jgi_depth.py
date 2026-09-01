#!/usr/bin/env python3
"""
Build the released 2-02_mag_metadata.tsv table from completed
pb-metagenomics-tools/HiFi-MAG-Pipeline outputs and NCBI accession tables.

This script reconstructs the final public MAG metadata table for the
6,505 primary bacterial MAGs.

Upstream inputs
---------------
1) pb-metagenomics-tools/HiFi-MAG-Pipeline outputs

    path/to/HiFi-MAG-Pipeline/
    ├── 2-bam/
    │   ├── <sample>.JGI.depth.txt
    │   └── ...
    └── 8-summary/
        ├── <sample>/
        │   └── <sample>.HiFi_MAG.summary.txt
        └── ...

The sample-level HiFi_MAG.summary.txt files provide:
    Name
    Completeness
    Contamination
    circular
    Contig_Number
    Contig_Names
    Genome_Size
    GC_Content
    classification

Thus, circularity is read directly from the corresponding 8-summary file for
each sample. No separately merged private metadata file is required.

The *.JGI.depth.txt files are upstream outputs from
jgi_summarize_bam_contig_depths (MetaBAT2) and provide contigLen and
totalAvgDepth for calculation of MAG-level average sequencing depth.

2) NCBI genome accession tables

    path/to/ncbi_genome_info/
    ├── genome-info-1.tsv
    ├── genome-info-2.tsv
    └── ...

Each genome-info-*.tsv file must contain:
    sample_name
    genome_acc

The sample_name values are expected to match the primary MAG identifiers.
Valid genome accessions are merged into the final table; MAGs without an
assigned accession are reported as NA.

Requirements
------------
- Python 3 (standard library only)
- Completed pb-metagenomics-tools/HiFi-MAG-Pipeline outputs
- NCBI genome-info-*.tsv accession tables

Relevant upstream software used by HiFi-MAG-Pipeline includes:
- MetaBAT2 v2.15 / jgi_summarize_bam_contig_depths
- CheckM2 v1.0.2
- GTDB-Tk v2.4.0 with GTDB release R220

MAG average depth
-----------------
For each MAG:

    sum(contigLen * totalAvgDepth) / sum(contigLen)

where the sum is taken across all contigs assigned to that MAG.

Released output
---------------
The resulting 2-02_mag_metadata.tsv contains:

    mag_id
    sample_id
    completeness
    contamination
    circular
    contig_count
    average_depth
    genome_size_bp
    gc_content
    gtdb_taxonomy
    ncbi_genome_accession

Example
-------
    python3 build_mag_metadata_with_jgi_depth.py \
        --pipeline-root path/to/HiFi-MAG-Pipeline \
        --genome-info-dir path/to/ncbi_genome_info \
        --out 2-02_mag_metadata.tsv
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


STUDY_SAMPLES = [
    "CW1",
    "CW2",
    "Con",
    "GY1",
    "MP1",
    "MP2",
    "NH2",
    "NH3",
    "PH1",
    "PH2",
    "US1",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build the released 2-02_mag_metadata.tsv table from "
            "HiFi-MAG-Pipeline outputs and NCBI accession tables."
        )
    )
    parser.add_argument(
        "--pipeline-root",
        type=Path,
        required=True,
        help=(
            "Path to HiFi-MAG-Pipeline containing 2-bam/ and 8-summary/."
        ),
    )
    parser.add_argument(
        "--genome-info-dir",
        type=Path,
        required=True,
        help="Directory containing NCBI genome-info-*.tsv files.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("2-02_mag_metadata.tsv"),
        help="Output TSV path. Default: 2-02_mag_metadata.tsv",
    )
    parser.add_argument(
        "--expected-mags",
        type=int,
        default=6505,
        help="Expected number of primary MAGs. Default: 6505",
    )
    parser.add_argument(
        "--expected-accessions",
        type=int,
        default=4031,
        help=(
            "Expected number of MAGs with NCBI accessions. "
            "Default: 4031. Set to 0 to disable this check."
        ),
    )
    return parser.parse_args()


def split_contigs(value: str) -> list[str]:
    """Split Contig_Names robustly across common delimiters."""
    value = (value or "").strip()
    if not value:
        return []
    return [x for x in re.split(r"[,;|\s]+", value) if x]


def load_jgi_depth(path: Path) -> dict[str, tuple[float, float]]:
    """Load contig length and totalAvgDepth from a JGI depth table."""
    if not path.is_file():
        raise FileNotFoundError(path)

    depth_by_contig: dict[str, tuple[float, float]] = {}

    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"contigName", "contigLen", "totalAvgDepth"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise RuntimeError(
                f"{path}: required columns missing: {sorted(missing)}; "
                f"found={reader.fieldnames}"
            )

        for row in reader:
            name = row["contigName"].strip()

            if not name:
                raise RuntimeError(f"{path}: empty contigName")
            if name in depth_by_contig:
                raise RuntimeError(f"{path}: duplicate contigName: {name}")

            try:
                contig_len = float(row["contigLen"])
                total_avg_depth = float(row["totalAvgDepth"])
            except ValueError as exc:
                raise RuntimeError(
                    f"{path}: invalid contigLen/totalAvgDepth for {name}"
                ) from exc

            if contig_len <= 0:
                raise RuntimeError(
                    f"{path}: non-positive contigLen for {name}: {contig_len}"
                )
            if total_avg_depth < 0:
                raise RuntimeError(
                    f"{path}: negative totalAvgDepth for {name}: "
                    f"{total_avg_depth}"
                )

            depth_by_contig[name] = (contig_len, total_avg_depth)

    return depth_by_contig


def genome_info_sort_key(path: Path) -> tuple[int, str]:
    """Sort genome-info files numerically when a number is present."""
    match = re.search(r"(\d+)", path.stem)
    if match:
        return int(match.group(1)), path.name
    return 10**12, path.name


def load_ncbi_accessions(
    genome_info_dir: Path,
) -> tuple[dict[str, str], list[Path]]:
    """Load MAG-to-GenBank accession mappings from NCBI tables."""
    files = sorted(
        genome_info_dir.glob("genome-info-*.tsv"),
        key=genome_info_sort_key,
    )

    if not files:
        raise RuntimeError(
            f"No genome-info-*.tsv files found in {genome_info_dir}"
        )

    accession_by_mag: dict[str, str] = {}

    for path in files:
        with path.open(
            newline="",
            encoding="utf-8-sig",
            errors="replace",
        ) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"sample_name", "genome_acc"}
            missing = required - set(reader.fieldnames or [])

            if missing:
                raise RuntimeError(
                    f"{path}: required columns missing: {sorted(missing)}; "
                    f"found={reader.fieldnames}"
                )

            for row in reader:
                mag_id = (row.get("sample_name") or "").strip()
                accession = (row.get("genome_acc") or "").strip()

                if not mag_id or not accession:
                    continue

                # Keep only accession-like values. Error/status strings in
                # submission tables are intentionally ignored.
                if re.search(
                    r"ERR|ERROR|FAIL|FAILED",
                    accession,
                    flags=re.IGNORECASE,
                ):
                    continue

                if not re.fullmatch(r"[A-Z]+[0-9]+", accession):
                    continue

                previous = accession_by_mag.get(mag_id)
                if previous is not None and previous != accession:
                    raise RuntimeError(
                        f"Conflicting NCBI accessions for {mag_id}: "
                        f"{previous} vs {accession}"
                    )

                accession_by_mag[mag_id] = accession

    return accession_by_mag, files


def normalize_circular(value: str, mag_id: str) -> str:
    """Normalize circularity to yes/no."""
    circular = (value or "").strip().lower()

    if circular not in {"yes", "no"}:
        raise RuntimeError(
            f"{mag_id}: unexpected circular value: {value!r}"
        )

    return circular


def main() -> None:
    args = parse_args()

    pipeline_root = args.pipeline_root.expanduser().resolve()
    analysis_dir = pipeline_root / "8-summary"
    jgi_dir = pipeline_root / "2-bam"
    genome_info_dir = args.genome_info_dir.expanduser().resolve()
    out_path = args.out.expanduser()

    for path, label in [
        (pipeline_root, "pipeline root"),
        (analysis_dir, "8-summary directory"),
        (jgi_dir, "2-bam directory"),
        (genome_info_dir, "NCBI genome-info directory"),
    ]:
        if not path.is_dir():
            raise SystemExit(f"ERROR: {label} not found: {path}")

    ncbi_map, genome_info_files = load_ncbi_accessions(genome_info_dir)

    rows_out: list[dict[str, str]] = []
    seen_mag_ids: set[str] = set()

    missing_contigs: list[tuple[str, list[str]]] = []
    count_mismatch: list[tuple[str, int, int]] = []
    genome_size_mismatch: list[tuple[str, int, int]] = []

    for sample in STUDY_SAMPLES:
        summary_path = (
            analysis_dir / sample / f"{sample}.HiFi_MAG.summary.txt"
        )
        jgi_path = jgi_dir / f"{sample}.JGI.depth.txt"

        if not summary_path.is_file():
            raise FileNotFoundError(summary_path)
        if not jgi_path.is_file():
            raise FileNotFoundError(jgi_path)

        print(
            f"[INFO] sample={sample} "
            f"summary={summary_path} "
            f"jgi_depth={jgi_path}"
        )

        jgi = load_jgi_depth(jgi_path)

        with summary_path.open(
            newline="",
            encoding="utf-8-sig",
            errors="replace",
        ) as fh:
            reader = csv.DictReader(fh, delimiter="\t")

            required = {
                "Name",
                "Completeness",
                "Contamination",
                "circular",
                "Contig_Number",
                "Contig_Names",
                "Genome_Size",
                "GC_Content",
                "classification",
            }
            missing_columns = required - set(reader.fieldnames or [])

            if missing_columns:
                raise RuntimeError(
                    f"{summary_path}: required columns missing: "
                    f"{sorted(missing_columns)}; "
                    f"found={reader.fieldnames}"
                )

            for row in reader:
                bin_name = (row["Name"] or "").strip()

                if not bin_name:
                    raise RuntimeError(
                        f"{summary_path}: empty Name field"
                    )

                # HiFi-MAG-Pipeline sample summary files contain bin names.
                # Prefix with the sample ID to produce the stable primary MAG ID.
                mag_id = f"{sample}_{bin_name}"

                if mag_id in seen_mag_ids:
                    raise RuntimeError(
                        f"Duplicate primary MAG ID: {mag_id}"
                    )
                seen_mag_ids.add(mag_id)

                circular = normalize_circular(row["circular"], mag_id)

                contigs = split_contigs(row["Contig_Names"])

                try:
                    expected_n = int(float(row["Contig_Number"]))
                    expected_genome_size = int(float(row["Genome_Size"]))
                except ValueError as exc:
                    raise RuntimeError(
                        f"{mag_id}: invalid Contig_Number or Genome_Size"
                    ) from exc

                if len(contigs) != expected_n:
                    count_mismatch.append(
                        (mag_id, expected_n, len(contigs))
                    )
                    continue

                weighted_sum = 0.0
                total_len = 0.0
                missing_here: list[str] = []

                for contig in contigs:
                    if contig not in jgi:
                        missing_here.append(contig)
                        continue

                    length, depth = jgi[contig]
                    weighted_sum += length * depth
                    total_len += length

                if missing_here:
                    missing_contigs.append((mag_id, missing_here))
                    continue

                if total_len <= 0:
                    raise RuntimeError(
                        f"{mag_id}: total MAG contig length <= 0"
                    )

                jgi_total_len = int(round(total_len))
                if jgi_total_len != expected_genome_size:
                    genome_size_mismatch.append(
                        (
                            mag_id,
                            expected_genome_size,
                            jgi_total_len,
                        )
                    )
                    continue

                average_depth = weighted_sum / total_len

                rows_out.append({
                    "mag_id": mag_id,
                    "sample_id": sample,
                    "completeness": row["Completeness"],
                    "contamination": row["Contamination"],
                    "circular": circular,
                    "contig_count": row["Contig_Number"],
                    "average_depth": f"{average_depth:.6f}",
                    "genome_size_bp": row["Genome_Size"],
                    "gc_content": row["GC_Content"],
                    "gtdb_taxonomy": row["classification"],
                    "ncbi_genome_accession": ncbi_map.get(
                        mag_id,
                        "NA",
                    ),
                })

    if count_mismatch:
        print(
            "ERROR: Contig_Number does not match parsed Contig_Names.",
            file=sys.stderr,
        )
        for item in count_mismatch[:10]:
            print(item, file=sys.stderr)
        raise SystemExit(1)

    if missing_contigs:
        print(
            "ERROR: Some MAG contigs were not found in JGI depth tables.",
            file=sys.stderr,
        )
        for mag_id, contigs in missing_contigs[:10]:
            print(mag_id, contigs[:10], file=sys.stderr)
        raise SystemExit(1)

    if genome_size_mismatch:
        print(
            "ERROR: Genome_Size differs from summed JGI contig lengths.",
            file=sys.stderr,
        )
        for item in genome_size_mismatch[:10]:
            print(item, file=sys.stderr)
        raise SystemExit(1)

    if len(rows_out) != args.expected_mags:
        raise RuntimeError(
            f"Expected {args.expected_mags} MAGs, "
            f"recovered {len(rows_out)}"
        )

    catalogue_ids = {row["mag_id"] for row in rows_out}
    unknown_ncbi = sorted(set(ncbi_map) - catalogue_ids)

    if unknown_ncbi:
        print(
            f"WARNING: {len(unknown_ncbi)} NCBI sample_name values do not "
            "match the primary MAG catalogue.",
            file=sys.stderr,
        )
        for mag_id in unknown_ncbi[:10]:
            print(f"  {mag_id}", file=sys.stderr)

    matched_accessions = sum(
        row["ncbi_genome_accession"] != "NA"
        for row in rows_out
    )

    if (
        args.expected_accessions > 0
        and matched_accessions != args.expected_accessions
    ):
        raise RuntimeError(
            f"Expected {args.expected_accessions} NCBI accessions, "
            f"matched {matched_accessions}"
        )

    fieldnames = [
        "mag_id",
        "sample_id",
        "completeness",
        "contamination",
        "circular",
        "contig_count",
        "average_depth",
        "genome_size_bp",
        "gc_content",
        "gtdb_taxonomy",
        "ncbi_genome_accession",
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with out_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows_out)

    circular_yes = sum(row["circular"] == "yes" for row in rows_out)

    print()
    print("Final primary MAG metadata")
    print("--------------------------")
    print(f"Output                  : {out_path}")
    print(f"Total MAGs              : {len(rows_out)}")
    print(f"Circular MAGs           : {circular_yes}")
    print(f"Non-circular MAGs       : {len(rows_out) - circular_yes}")
    print(f"NCBI accessions matched : {matched_accessions}")
    print(f"No NCBI accession       : {len(rows_out) - matched_accessions}")
    print(f"Genome-info files read  : {len(genome_info_files)}")
    print("Depth source            : JGI totalAvgDepth")
    print(
        "MAG average depth       : "
        "sum(contigLen * totalAvgDepth) / sum(contigLen)"
    )


if __name__ == "__main__":
    main()

