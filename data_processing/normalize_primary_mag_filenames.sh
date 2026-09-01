#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -------------------------------------------------------------------------
# Normalize primary MAG filenames
# -------------------------------------------------------------------------
#
# Rename MAG FASTA files generated under:
#
#   8-summary/<sample>/MAGs/*.fa
#
# so that every filename begins with its sample identifier.
#
# Example:
#
#   8-summary/CW1/MAGs/complete.105.fa
#       -> 8-summary/CW1/MAGs/CW1_complete.105.fa
#
# Files that already have the expected sample prefix are left unchanged.
#
# Usage
# -----
# Option 1: edit/export MAG_ROOT, then run:
#
#   MAG_ROOT=/path/to/HiFi-MAG-Pipeline/8-summary \
#     bash normalize_primary_mag_filenames.sh
#
# Option 2: if the default placeholder below has been replaced:
#
#   bash normalize_primary_mag_filenames.sh
#
# Optional dry run:
#
#   DRY_RUN=1 \
#   MAG_ROOT=/path/to/HiFi-MAG-Pipeline/8-summary \
#     bash normalize_primary_mag_filenames.sh
#
# Notes
# -----
# - The script is idempotent: already-normalized filenames are skipped.
# - Existing destination files are never overwritten.
# - Only files matching 8-summary/<sample>/MAGs/*.fa are considered.
# - This naming normalization is intended to be performed before downstream
#   MAG metadata, MIMAG, ANI95, and other catalogue-level analyses.

MAG_ROOT="${MAG_ROOT:-path/to/HiFi-MAG-Pipeline/8-summary}"
DRY_RUN="${DRY_RUN:-0}"

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { echo "ERROR: $*" >&2; exit 1; }

[[ -d "$MAG_ROOT" ]] || die "MAG_ROOT not found: $MAG_ROOT"

renamed=0
skipped=0
collisions=0
total=0

log "[INFO] MAG_ROOT = $MAG_ROOT"
log "[INFO] DRY_RUN  = $DRY_RUN"

for sample_dir in "$MAG_ROOT"/*; do
    [[ -d "$sample_dir" ]] || continue

    sample_id="$(basename "$sample_dir")"
    mag_dir="$sample_dir/MAGs"
    [[ -d "$mag_dir" ]] || continue

    for f in "$mag_dir"/*.fa; do
        [[ -f "$f" ]] || continue

        total=$((total + 1))

        base="$(basename "$f")"

        # Already normalized.
        if [[ "$base" == "${sample_id}_"* ]]; then
            echo "skip (already normalized): $f"
            skipped=$((skipped + 1))
            continue
        fi

        new="$mag_dir/${sample_id}_${base}"

        # Never overwrite an existing destination file.
        if [[ -e "$new" ]]; then
            echo "collision (target exists): $new" >&2
            collisions=$((collisions + 1))
            continue
        fi

        if [[ "$DRY_RUN" == "1" ]]; then
            echo "would rename: $f -> $new"
        else
            mv -- "$f" "$new"
            echo "renamed: $f -> $new"
        fi

        renamed=$((renamed + 1))
    done
done

echo
echo "============================================================"
echo "Primary MAG filename normalization summary"
echo "============================================================"
echo "MAG_ROOT              : $MAG_ROOT"
echo "MAG FASTAs examined   : $total"
echo "Renamed / would rename: $renamed"
echo "Already normalized    : $skipped"
echo "Collisions             : $collisions"
echo "============================================================"

if (( collisions > 0 )); then
    die "One or more destination filenames already existed; review collisions."
fi

if (( total == 0 )); then
    die "No MAG FASTA files found under $MAG_ROOT/<sample>/MAGs/*.fa"
fi

if [[ "$DRY_RUN" == "1" ]]; then
    log "[DONE] Dry run completed; no files were modified."
else
    log "[DONE] Primary MAG filename normalization completed."
fi

