library(dplyr)
library(tibble)
library(flextable)
library(officer)

# -------------------------------------------------------------------------
# Data
# -------------------------------------------------------------------------

table2 <- tribble(
  ~`Data group`, ~`File`, ~`Description`,

  "1. Sequencing and assembly",
  "1-01_sequencing_assembly_summary.tsv",
  "Sample-level PacBio HiFi sequencing yield, read N50, assembly statistics, primary MAG recovery, and MAG-associated mapped-base metrics",

  "1. Sequencing and assembly",
  "1-02_mag_associated_coverage_qc.tsv",
  "Quality-control summary supporting MAG-associated mapped-base calculations",

  "2. Primary bacterial MAG catalogue",
  "2-01_sample_metadata.tsv",
  "Source-soil, cultivation, rhizosphere sampling, and sequencing metadata for the 11 composite samples",

  "2. Primary bacterial MAG catalogue",
  "2-02_mag_metadata.tsv",
  "MAG-level genome-quality metrics, genome statistics, GTDB taxonomy, and GenBank accession mapping",

  "2. Primary bacterial MAG catalogue",
  "2-03_bacterial_MAGs_6505_by_sample.tar.gz",
  "Sample-resolved FASTA files for the 6,505 primary bacterial MAGs",

  "3. MIMAG quality assessment",
  "3-01_mimag_by_sample.tsv",
  "Sample-level summary of MIMAG quality classifications",

  "3. MIMAG quality assessment",
  "3-02_mimag_qc.tsv",
  "Quality-control summary supporting the MIMAG assessment",

  "3. MIMAG quality assessment",
  "3-03_mimag_quality.tsv",
  "MAG-level CheckM2, rRNA, tRNA, and MIMAG quality results",

  "3. MIMAG quality assessment",
  "3-04_mimag_summary.tsv",
  "Catalogue-level summary of MIMAG quality classifications",

  "4. ANI95 catalogue",
  "4-01_ani95_cluster_membership.tsv",
  "ANI95 cluster membership for all 6,505 primary bacterial MAGs",

  "4. ANI95 catalogue",
  "4-02_ani95_cluster_summary.tsv",
  "Cluster-level summary including representative MAG, cluster size, sample occupancy, member MAGs, and member samples",

  "4. ANI95 catalogue",
  "4-03_ani95_representative_metadata.tsv",
  "Metadata for the 4,032 ANI95 representative MAGs",

  "4. ANI95 catalogue",
  "4-04_ani95_cluster_by_sample_presence.tsv",
  "ANI95 cluster-by-sample presence/absence matrix",

  "4. ANI95 catalogue",
  "4-05_ani95_sample_cluster_counts.tsv",
  "Number of distinct ANI95 clusters detected in each composite rhizosphere sample",

  "4. ANI95 catalogue",
  "4-06_representative_MAGs_ANI95.tar.gz",
  "FASTA files for the 4,032 ANI95 representative MAGs"
)

# -------------------------------------------------------------------------
# Format
# -------------------------------------------------------------------------

# Rows ending each data group
group_end_rows <- cumsum(rle(table2$`Data group`)$lengths)
group_end_rows <- head(group_end_rows, -1)

ft <- flextable(table2) |>
  theme_booktabs() |>
  merge_v(j = 1) |>
  hline(
    i = group_end_rows,
    border = officer::fp_border(width = 0.5),
    part = "body"
  ) |>
  bold(part = "header") |>
  valign(j = 1, valign = "center", part = "body") |>
  valign(j = 2:3, valign = "top", part = "body") |>
  align(j = 1:3, align = "left", part = "all") |>
  width(j = 1, width = 1.45) |>
  width(j = 2, width = 2.65) |>
  width(j = 3, width = 3.60) |>
  fontsize(size = 9, part = "all") |>
  padding(padding = 3, part = "all") |>
  set_caption(
    caption = "Table 2. Data files deposited in Zenodo for the bacterial MAG catalogue."
  )

ft
# -------------------------------------------------------------------------
# Export
# -------------------------------------------------------------------------

dir.create("figure_and_table", showWarnings = FALSE)

doc <- read_docx() |>
  body_add_flextable(ft)

print(
  doc,
  target = "figure_and_table/Table2.docx"
)
