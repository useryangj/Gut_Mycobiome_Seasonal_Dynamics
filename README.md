# ITS Mycobiome Seasonal Dynamics

This repository contains R workflows for analyzing fungal communities from ITS amplicon sequencing data, with additional analyses of seasonal variation, monthly dynamics, macro-fungi, diet composition, and climatic variables.

> **Public article URL:** To be added  
> **Public data URL:** To be added  
> **Project status:** Analysis scripts under preparation

## Overview

The main workflow includes:

1. Quality control, primer removal, and DADA2 denoising of paired-end FASTQ files.
2. ASV inference, chimera removal, and taxonomic assignment against the UNITE database.
3. Generation of raw-count, rarefied, and relative-abundance feature tables.
4. Aggregation of abundance at the phylum, family, genus, and species levels.
5. Alpha-diversity analysis by season and month using Kruskal–Wallis tests and Dunn tests with Benjamini–Hochberg correction.
6. PCoA, PERMANOVA, and PERMDISP based on Bray–Curtis and Robust Aitchison distances.
7. Visualization of community composition, rarefaction curves, Venn diagrams, heatmaps, and Sankey diagrams.
8. Screening and seasonal/monthly analysis of macro-fungi-related taxa.
9. RDA, variance partitioning, sensitivity analyses, and core-genus correlation analyses using monthly diet and climate data.

## Repository Structure

| File | Description |
| --- | --- |
| [`ITS.R`](ITS.R) | ITS preprocessing, taxonomic profiling, diversity statistics, community analyses, and visualizations |
| [`macro.R`](macro.R) | Screening and visualization of macro-fungi-related genera and species |
| [`Monthly_RDA.R`](Monthly_RDA.R) | Monthly RDA linking fungal communities with diet composition, temperature, and humidity |

## Analysis Workflow

### 1. ITS Processing

`ITS.R` uses `dada2`, `ShortRead`, and `Biostrings` to:

- Read paired-end FASTQ files.
- Remove reads containing invalid bases.
- Remove forward and reverse primers using `cutadapt`.
- Perform quality filtering, error-model learning, and DADA2 denoising.
- Merge paired reads.
- Remove chimeric sequences.
- Assign taxonomy using a UNITE reference database.
- Export ASV tables, taxonomic annotations, and tables aggregated at multiple taxonomic levels.

Representative intermediate files include:

- `quality_control.txt`
- `feature_table.txt`
- `feature_table_relative.txt`
- `feature_table_unflatted.txt`
- `phylum_table.txt`
- `family_table.txt`
- `genus_table.txt`
- `species_table.txt`

### 2. Alpha Diversity

The workflow calculates:

- Shannon index
- Simpson index
- Chao1
- ACE

Grouping variables include:

- Seasons: Spring, Summer, Autumn, and Winter
- Months: Jan–Dec

Overall differences are tested using Kruskal–Wallis tests. Pairwise differences are evaluated using Dunn tests with Benjamini–Hochberg multiple-testing correction.

### 3. Beta Diversity

Two distance-based analysis frameworks are included:

- **Bray–Curtis distance:** monthly and seasonal PERMANOVA, PERMDISP, and PCoA.
- **Robust Aitchison distance:** compositional-data analysis of community differences.

Typical outputs include:

- Overall PERMANOVA results
- Pairwise PERMANOVA results
- PERMDISP results
- PCoA coordinates and figures
- BH-adjusted `q` values

### 4. Macro-Fungi Analysis

`macro.R` uses the taxon list in `macro_fungi.txt` to extract relevant genera and species from the abundance tables. The script produces:

- Macro-fungi genus and species abundance tables
- Monthly and seasonal relative-abundance summaries
- Bar plots
- Sankey diagrams
- Monthly abundance-flow plots

### 5. Monthly RDA

`Monthly_RDA.R` treats the 12 months as analysis units and relates genus-level fungal communities to the following diet and environmental variables:

- Young leaves
- Mature leaves
- Fruits
- Flowers
- Invertebrate
- Seeds
- Roots
- Macro-fungi
- Barks
- Fallen leaves
- Others
- Temperature
- Humidity

The script generates single-factor RDA results, a reduced RDA model, marginal effects, VIF values, variance partitioning, leave-one-month-out sensitivity analyses, circular-shift sensitivity analyses, and a core-genus correlation heatmap.

The default output directory is:

```text
Monthly_RDA_final1/
```

## Software Requirements

R 4.x is recommended. Install the required R packages with:

```r
install.packages(c(
  "vegan", "ggplot2", "ggpubr", "FSA", "cowplot",
  "ggvenn", "reshape2", "dplyr", "readxl", "tibble",
  "pheatmap", "ggrepel", "permute", "jsonlite",
  "htmlwidgets", "networkD3", "ggalluvial", "microeco"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("dada2", "ShortRead", "Biostrings"))
```

The ITS preprocessing workflow also requires:

- `cutadapt`
- A UNITE ITS reference database

## Data and Path Configuration

The scripts currently contain absolute paths from the original analysis environment, including:

```text
C:/Users/xuzy0/Desktop/xzy/split
D:/split
```

Before running the workflow on another computer, update:

1. The raw FASTQ input directory.
2. The path to the `cutadapt` executable.
3. The UNITE reference database path.
4. Feature-table, macro-fungi, and environmental-data paths.
5. Output directories.

`Monthly_RDA.R` additionally requires:

- `genus_table1.txt`
- `每月食物结构.xlsx`
- `每月温度湿度.xlsx`

The two Excel filenames are retained in their original form so that they can be matched directly to the current analysis files.

## Recommended Execution Order

```text
1. ITS.R
2. macro.R
3. Monthly_RDA.R
```

After each step, inspect the quality-control table, sample counts, group assignments, and generated outputs before proceeding.

## Statistical Considerations

- The monthly RDA uses 12 months as observations, which limits statistical power.
- Month is temporally ordered; unrestricted permutations may not fully account for temporal autocorrelation.
- RDA, correlation, and variance-partitioning results should be interpreted as exploratory associations rather than causal effects.
- Significant PERMANOVA results should be interpreted together with PERMDISP results to assess the influence of differences in within-group dispersion.
- Raw counts, rarefied counts, and relative abundances must not be interchanged without checking the assumptions of each analysis.
- If samples represent repeated individuals, use a design-aware restricted permutation scheme rather than unrestricted permutations.

## Article and Data Availability

The public article and data links will be added later. The following information should be completed:

- Article title
- Authors
- Journal or preprint platform
- DOI or public URL
- Sequence-data accession number
- Mapping between repository outputs and figures or analytical sections in the article

## Citation

When using this repository, please cite the software and databases used in the analysis, including DADA2, UNITE, vegan, cutadapt, and all other packages used to generate the reported results.
