# Cucumber Rhizosphere Bacterial MAG Catalogue Analysis

This repository contains the analysis and data-processing scripts used to reproduce Figures 1–4, Tables 1–2, and selected derived data products presented in the associated *Scientific Data* Data Descriptor describing a genome-resolved bacterial metagenomic catalogue from cucumber rhizospheres amended with coastal plant-habitat soils.

Genome assembly, MAG recovery, genome-quality assessment, and taxonomic classification were performed using publicly available bioinformatics software as described in the manuscript. This repository additionally provides scripts for selected downstream catalogue-processing steps, including construction of the ANI95 representative catalogue and reconstruction of the released primary-MAG metadata table.

## Data availability

The data files used by the scripts in this repository are available from Zenodo:

**DOI: [10.5281/zenodo.22077545](https://doi.org/10.5281/zenodo.22077545)**

The Zenodo deposition includes:

- sequencing and assembly summary data;
- sample metadata;
- metadata for the 6,505 primary MAGs;
- the sample-resolved primary MAG FASTA archive;
- MIMAG quality-assessment results;
- ANI95 cluster membership and cluster summaries;
- metadata for the 4,032 ANI95 representative MAGs; and
- the ANI95 representative MAG FASTA archive.

See the Zenodo `README.md` for the complete description of each deposited file.

## Repository structure

```text
cucumber-rhizosphere-bacterial-mag-catalog-analysis/
├── Fig 1.R
├── Fig 2.R
├── Fig 3.R
├── Fig 4.R
├── README.md
├── Table 1.R
├── Table 2.R
├── data_processing/
│   ├── build_ani95_catalogue_with_dRep.sh
│   ├── build_mag_metadata_with_jgi_depth.py
│   └── build_mimag_quality_tables_with_Barrnap_tRNAscan-SE.sh
└── fig1_asset/
    └── Fig1C_workflow.png
```

The R scripts correspond directly to Figures 1–4 and Tables 1–2 in the associated Data Descriptor:

- `Fig 1.R` — coastal source-soil sampling sites, experimental design, and genome-resolved metagenomic workflow;
- `Fig 2.R` — HiFi sequencing yield, assembly characteristics, primary MAG recovery, and MAG-associated HiFi bases;
- `Fig 3.R` — MAG completeness and contamination, contiguity, and MIMAG high-quality criteria;
- `Fig 4.R` — catalogue subset relationships and ANI95 cluster structure, occupancy, and sharing among samples;
- `Table 1.R` — sample and source-soil metadata table; and
- `Table 2.R` — data files deposited in Zenodo for the bacterial MAG catalogue.

The `data_processing/` directory contains scripts used to reconstruct selected released data products:

- `build_ani95_catalogue_with_dRep.sh` — constructs the ANI95 catalogue from the 6,505 primary MAG nucleotide FASTA files using dRep and generates stable ANI95 identifiers, cluster membership tables, representative metadata, sample-presence tables, and representative MAG FASTA links;
- `build_mag_metadata_with_jgi_depth.py` — reconstructs the released `2-02_mag_metadata.tsv` table from completed `pb-metagenomics-tools/HiFi-MAG-Pipeline` outputs and calculates MAG average depth as a contig-length-weighted mean of JGI `totalAvgDepth` values; and
- `build_mimag_quality_tables_with_Barrnap_tRNAscan-SE.sh` — runs Barrnap and tRNAscan-SE on the 6,505 primary MAGs, combines the RNA annotations with CheckM2 completeness and contamination values from the released `2-02_mag_metadata.tsv` table, applies the MIMAG high-quality criteria used in the study, and generates the released MIMAG quality, summary, QC, and per-sample tables.

`fig1_asset/` contains the workflow schematic used in Figure 1C. The `figure_and_table/` directory is used for generated manuscript outputs and is not tracked in this repository.

## Input data

### Figure and table scripts

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

### ANI95 catalogue processing

`data_processing/build_ani95_catalogue_with_dRep.sh` expects the primary MAG FASTA files in the final summary structure produced by `pb-metagenomics-tools/HiFi-MAG-Pipeline`:

```text
path/to/8-summary/
├── <sample>/
│   └── MAGs/
│       ├── *.fa
│       └── ...
└── ...
```

It also requires the released `2-02_mag_metadata.tsv` table, which provides the CheckM2 completeness and contamination values used by dRep through `--genomeInfo`.

Example:

```bash
MAG_ROOT=path/to/8-summary \
MAG_METADATA=path/to/2-02_mag_metadata.tsv \
OUT_ROOT=path/to/ANI95_output \
bash data_processing/build_ani95_catalogue_with_dRep.sh
```

### Primary MAG metadata reconstruction

`data_processing/build_mag_metadata_with_jgi_depth.py` is intended to be run on the completed output of `pb-metagenomics-tools/HiFi-MAG-Pipeline` after processing the study FASTA/FASTQ inputs.

The relevant upstream files are:

```text
path/to/HiFi-MAG-Pipeline/
├── 2-bam/
│   ├── <sample>.JGI.depth.txt
│   └── ...
└── 8-summary/
    ├── <sample>/
    │   └── <sample>.HiFi_MAG.summary.txt
    └── ...
```

The script reconstructs the released `2-02_mag_metadata.tsv` table. MAG average depth is calculated as:

```text
sum(contigLen × totalAvgDepth) / sum(contigLen)
```

where `contigLen` and `totalAvgDepth` are taken from the sample-specific JGI depth table for contigs assigned to each MAG.

Example:

```bash
python3 data_processing/build_mag_metadata_with_jgi_depth.py \
    --pipeline-root path/to/HiFi-MAG-Pipeline \
    --out 2-02_mag_metadata.tsv
```

### MIMAG quality processing

`data_processing/build_mimag_quality_tables_with_Barrnap_tRNAscan-SE.sh` uses the same primary MAG FASTA layout:

```text
path/to/8-summary/
├── <sample>/
│   └── MAGs/
│       ├── *.fa
│       └── ...
└── ...
```

It also requires the released `2-02_mag_metadata.tsv` table for CheckM2 completeness and contamination values. The script runs Barrnap and tRNAscan-SE, then applies the MIMAG high-quality criteria used in the study:

```text
CheckM2 completeness > 90%
CheckM2 contamination < 5%
complete 5S, 16S, and 23S rRNA genes present
tRNAs for at least 18 of the 20 canonical amino-acid types
```

For the tRNA criterion, pseudogene predictions are excluded, `Ile2` is normalized to `Ile`, `fMet` is normalized to `Met`, and `SeC`, `Sup`, and `Undet` are excluded from the canonical 20-amino-acid count. Barrnap hits annotated as partial are not counted toward the rRNA criterion.

Example:

```bash
MAG_ROOT=path/to/8-summary \
MAG_METADATA=path/to/2-02_mag_metadata.tsv \
OUT_ROOT=path/to/mimag_quality_output \
bash data_processing/build_mimag_quality_tables_with_Barrnap_tRNAscan-SE.sh
```

## Requirements

### Figure and table scripts

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

### Data-processing scripts

`build_ani95_catalogue_with_dRep.sh` requires:

```text
dRep
fastANI
Mash
Python 3
```

The script uses the existing CheckM2 completeness and contamination values provided in the released `2-02_mag_metadata.tsv` table.

`build_mimag_quality_tables_with_Barrnap_tRNAscan-SE.sh` requires:

```text
Barrnap
tRNAscan-SE
Python 3
```

The script uses the CheckM2 completeness and contamination values provided in the released `2-02_mag_metadata.tsv` table and does not rerun CheckM2.

`build_mag_metadata_with_jgi_depth.py` requires only Python 3 standard-library modules at execution time, together with completed `pb-metagenomics-tools/HiFi-MAG-Pipeline` outputs. The upstream pipeline outputs used by this script include results generated with `jgi_summarize_bam_contig_depths` from MetaBAT2, CheckM2, and GTDB-Tk.

Software versions and key parameters used for the study are reported in the associated Data Descriptor.

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

This repository supports reproducibility at two levels:

1. reconstruction of selected released catalogue products, including primary-MAG metadata, MIMAG quality tables, and the ANI95 representative catalogue, through scripts in `data_processing/`; and
2. reproduction of Figures 1–4 and Tables 1–2 from the released Zenodo data products.

The upstream genome assembly, binning, genome-quality assessment, and taxonomic-classification tools are publicly available software and are documented with versions and key parameters in the associated Data Descriptor.

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

The analysis and data-processing code in this repository, including the scripts in `data_processing/`, is released under the MIT License. External software invoked by these scripts (for example, dRep, fastANI, Mash, Barrnap, tRNAscan-SE, MetaBAT2, CheckM2, and GTDB-Tk) remains subject to its own respective license.

The figure asset in `fig1_asset/` is licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

The associated MAG catalogue and derived data products are distributed separately through Zenodo under CC BY 4.0:

https://doi.org/10.5281/zenodo.22077545

