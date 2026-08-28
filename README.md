# Cucumber Rhizosphere Bacterial MAG Catalogue Analysis

This repository contains the R scripts used to reproduce Figures 1–4 and Tables 1–2 presented in the associated *Scientific Data* Data Descriptor describing a genome-resolved bacterial metagenomic catalogue from cucumber rhizospheres amended with coastal plant-habitat soils.

Genome reconstruction, MAG recovery, genome-quality assessment, taxonomic classification, and ANI95 clustering were performed using publicly available bioinformatics software as described in the manuscript. The resulting primary MAG catalogue and derived data products are deposited in Zenodo and are not duplicated in this repository.

## Data availability

The data files used by the scripts in this repository are available from Zenodo:

**DOI: [10.5281/zenodo.22077545](https://doi.org/10.5281/zenodo.22077545)**

Please refer to the Zenodo record for the primary and derived data products generated using the open-source bioinformatics workflows described in the associated Data Descriptor.

The Zenodo deposition includes:

- sequencing and assembly summary data;
- sample metadata;
- metadata for the 6,505 primary MAGs;
- the sample-resolved primary MAG FASTA archive;
- MIMAG quality-assessment results;
- ANI95 cluster membership and cluster summaries;
- metadata for the 4,032 ANI95 representative MAGs; and
- the ANI95 representative MAG FASTA archive.

## Repository structure

```text
cucumber-rhizosphere-bacterial-mag-catalog-analysis/
├── README.md
├── Fig 1.R
├── Fig 2.R
├── Fig 3.R
├── Fig 4.R
├── Table 1.R
├── Table 2.R
└── fig1_asset/
    └── Fig1C_workflow.png
```

The scripts correspond directly to Figures 1–4 and Tables 1–2 in the associated Data Descriptor:

- `Fig 1.R` — coastal source-soil sampling sites, experimental design, and genome-resolved metagenomic workflow;
- `Fig 2.R` — HiFi sequencing yield, assembly characteristics, primary MAG recovery, and MAG-associated HiFi bases;
- `Fig 3.R` — MAG completeness and contamination, contiguity, and MIMAG high-quality criteria;
- `Fig 4.R` — catalogue subset relationships and ANI95 cluster structure, occupancy, and sharing among samples;
- `Table 1.R` — sample and source-soil metadata table; and
- `Table 2.R` — data files deposited in Zenodo for the bacterial MAG catalogue.

`fig1_asset/` contains the workflow schematic used in Figure 1C. The `figure_and_table/` directory is used for generated outputs.

## Input data

Download the required TSV files from Zenodo and place them in the working directory used to run the R scripts.

The released data products follow a numbered naming scheme:

```text
1 — Sequencing and assembly
2 — Primary MAG catalogue
3 — MIMAG quality assessment
4 — ANI95 representative catalogue
```

Examples include:

```text
1-01_sequencing_assembly_summary.tsv
1-02_mag_associated_coverage_qc.tsv

2-01_sample_metadata.tsv
2-02_mag_metadata.tsv

3-01_mimag_by_sample.tsv
3-02_mimag_qc.tsv
3-03_mimag_quality.tsv
3-04_mimag_summary.tsv

4-01_ani95_cluster_membership.tsv
4-02_ani95_cluster_summary.tsv
4-03_ani95_representative_metadata.tsv
4-04_ani95_cluster_by_sample_presence.tsv
4-05_ani95_sample_cluster_counts.tsv
```

See the Zenodo `README.md` for the complete description of each deposited file.

## Requirements

The figure and table scripts were written in R and use standard CRAN packages. Depending on the script, required packages include:

```r
dplyr
tidyr
ggplot2
ggrepel
ggnewscale
ggpubr
ComplexUpset
VennDiagram
rnaturalearth
sf
png
showtext
sysfonts
officer
flextable
```

Install missing packages as needed, for example:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "ggnewscale",
  "ggpubr",
  "ComplexUpset",
  "VennDiagram",
  "rnaturalearth",
  "sf",
  "png",
  "showtext",
  "sysfonts",
  "officer",
  "flextable"
))
```

Some figures use `cairo_pdf()` for PDF export. A Cairo-enabled R installation is therefore recommended.

Figure 3 uses DejaVu Sans to support Unicode symbols used in figure labels. The script detects the installed font through the system `fc-match` command rather than using a user-specific font path.

## Reproducing the figures and tables

Run each script from the repository root after downloading the corresponding Zenodo input files.

For example:

```r
source("Fig 1.R")
source("Fig 2.R")
source("Fig 3.R")
source("Fig 4.R")
source("Table 1.R")
source("Table 2.R")
```

The scripts generate the manuscript figures and tables from the released data products. Output files are written to `figure_and_table/`.

The current manuscript organization is:

```text
Figure 1 — Sampling design and genome-resolved metagenomic workflow
Figure 2 — Sequencing, assembly, and primary MAG recovery
Figure 3 — MAG quality and genome characteristics
Figure 4 — ANI95 catalogue structure and relationships
Table 1  — Sample and source-soil metadata
Table 2  — Data files deposited in Zenodo
```

## Scope of this repository

This repository is intended to document and reproduce the analyses used to generate Figures 1–4 and Tables 1–2 in the Data Descriptor.

It does **not** provide a replacement implementation of the upstream bioinformatics software used to reconstruct and characterize the MAG catalogue. Genome assembly, binning, quality assessment, taxonomic classification, rRNA/tRNA assessment, and ANI95 clustering were performed with publicly available software, with software versions and key parameters reported in the associated manuscript.

For the resulting data products, use the Zenodo deposition:

**[10.5281/zenodo.22077545](https://doi.org/10.5281/zenodo.22077545)**

## Related resources

- **Zenodo dataset:** 10.5281/zenodo.22077545
- **NCBI BioProject:** PRJNA1272261
- **Associated article:** *Scientific Data* Data Descriptor (citation to be added after publication)

## Citation

If you use the MAG catalogue or derived data products, please cite the Zenodo dataset and the associated Data Descriptor once published.

For the deposited dataset:

> Cho, G. *Cucumber rhizosphere bacterial MAG catalogue from coastal plant-habitat soil amendments*. Zenodo (2026). https://doi.org/10.5281/zenodo.22077545

The final article citation will be added here after publication.

## License

The R analysis code in this repository is released under the MIT License.

The figure asset in `fig1_asset/` is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

The associated MAG catalogue and derived data products are distributed separately through Zenodo under CC BY 4.0:

https://doi.org/10.5281/zenodo.22077545

