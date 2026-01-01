# breseqConverter
An R package for parsing, comparing, and visualizing genomic mutation data from [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) output files. Designed for tracking genomic changes during bacterial strain engineering and adaptive laboratory evolution (ALE) experiments.
![Fig  S10](https://github.com/user-attachments/assets/4107d550-f704-4bcb-9437-47bf8f9eb504)

## Overview

`breseqConverter` streamlines the analysis of breseq HTML output by:

- **Extracting** mutation tables from breseq HTML reports
- **Comparing** genomic changes between strains or evolution time points
- **Tracking** sequential mutations across multiple rounds of strain improvement
- **Exporting** results to Excel for easy sharing and visualization

## Installation

```r
# Install from GitHub (recommended)
devtools::install_github("yourusername/breseqConverter")

# Or install locally
install.packages("path/to/breseqConverter", repos = NULL, type = "source")
```

## Dependencies

```r
install.packages(c("rvest", "dplyr", "openxlsx"))
```

## Quick Start

```r
library(breseqConverter)

# 1. Extract tables from breseq HTML output
extract_named_tables(directory = "path/to/breseq/output")

# 2. Compare two strains
comparison <- compare_named_tables(strain_A, strain_B)

# 3. Save results to Excel
save_comparison_results_to_excel(comparison, "comparison_results.xlsx")
```

## Core Functions

### `extract_named_tables()`

Parses breseq HTML output files and extracts mutation data into R data frames.

```r
extract_named_tables(directory = ".", file_pattern = "index.*\\.html")
```

**Extracts three table types:**
| Table | Description |
|-------|-------------|
| `Mutation predictions` | Confirmed SNPs, insertions, deletions |
| `Unassigned missing coverage evidence` | Potential deletions or low-coverage regions |
| `Unassigned new junction evidence` | Novel junction sequences (potential rearrangements) |

**Example:**
```r
# Extract from breseq output directory
extract_named_tables(directory = "./breseq_results/")

# Tables are saved to global environment with file names
# e.g., index_WT, index_mutant1, index_mutant2
```

### `compare_named_tables()`

Compares mutation profiles between two strains, identifying new, lost, and unchanged mutations.

```r
compare_named_tables(list_A, list_B)
```

**Returns a data frame with status column:**
- `New in B` — Mutations present only in strain B
- `Missing in B` — Mutations lost in strain B  
- `Unchanged or Modified` — Mutations present in both

**Example:**
```r
# Compare wild-type to evolved strain
results <- compare_named_tables(index_WT, index_evolved)

# View new mutations
results$`Mutation predictions` %>% 
  filter(status == "New in index_evolved")
```

### `compare_multiple_tables()`

Sequentially compares multiple strains for tracking mutations across evolution lineages.

```r
compare_multiple_tables(list_of_tables)
```

**Example:**
```r
# Track mutations across sequential rounds
lineage <- list(
  round1 = index_R1,
  round2 = index_R2,
  round3 = index_R3
)

all_comparisons <- compare_multiple_tables(lineage)
# Returns: Comparison_1_to_2, Comparison_2_to_3
```

### `save_comparison_results_to_excel()`

Exports comparison results to a multi-sheet Excel workbook.

```r
save_comparison_results_to_excel(comparison_results, "output.xlsx")
```

**Sheet names:**
- `Mutations` — Mutation predictions comparison
- `Missing` — Missing coverage evidence comparison
- `New junction` — New junction evidence comparison

### `save_results_to_excel()`

Saves each comparison as a separate Excel file (useful for multiple pairwise comparisons).

```r
save_results_to_excel(results, output_dir = "output_folder")
```

## Typical Workflow

### Scenario: Tracking mutations during strain engineering

```r
library(breseqConverter)

# Step 1: Run breseq on each strain
# (done outside R - generates index.html files)

# Step 2: Extract mutation tables
extract_named_tables(directory = "./breseq_output/")

# Step 3: Organize strains
strains <- list(
  WT = index_WT,
  engineering_round1 = index_R1,
  engineering_round2 = index_R2,
  final_strain = index_final
)

# Step 4: Sequential comparison
all_comparisons <- compare_multiple_tables(strains)

# Step 5: Export results
for (name in names(all_comparisons)) {
  save_comparison_results_to_excel(
    all_comparisons[[name]], 
    paste0(name, ".xlsx")
  )
}
```

### Scenario: Comparing parallel evolution experiments

```r
# Compare independently evolved strains to ancestor
ancestor <- index_ancestor

evolved_strains <- list(
  lineage_A = index_A,
  lineage_B = index_B,
  lineage_C = index_C
)

# Pairwise comparisons to ancestor
for (name in names(evolved_strains)) {
  comparison <- compare_named_tables(ancestor, evolved_strains[[name]])
  save_comparison_results_to_excel(comparison, paste0("vs_ancestor_", name, ".xlsx"))
}
```

## Output Interpretation

### Status Column Guide

| Status | Meaning | Biological Interpretation |
|--------|---------|---------------------------|
| `New in B` | Mutation found only in strain B | Acquired mutation |
| `Missing in B` | Mutation in A but not B | Reversion or loss |
| `Unchanged or Modified` | Present in both | Stable mutation (check if modified) |

### Key Columns

- **position**: Genomic coordinate of mutation
- **mutation**: Type and details (e.g., `A→G`, `Δ1 bp`)
- **annotation**: Gene/feature affected
- **description**: Functional annotation

## Utility Functions

### `get_screen_size()`

Helper function for adaptive plot sizing (macOS).

```r
dims <- get_screen_size(default_width = 12, default_height = 7, dpi = 96)
pdf("plot.pdf", width = dims$width, height = dims$height)
```

## Notes

- Input files must be breseq HTML output (`index.html`)
- Comparison is based on genomic `position` column
- For `Unassigned missing coverage evidence`, comparison uses `start` column
- Excel sheet names are truncated to 31 characters (Excel limitation)

## Citation

If you use this package, please cite:

```
[DNMB] DNMB: Deletion of non-canonical defense systems enable domestication of non-model Geobacillus for thermophilic platform engineering.
             Jae-Yoon Sung, Mun Hoe Lee, Jungsoo Park, Hyungbin Kim, Dariimaa Ganbat, Donggyu Kim, Hyun-Woo Cho, Sang Jae Lee, Seong Bo Kim, and Dong-Woo Lee. 2024.
             XXX, XXX, https://doi.org/XXX

[breseq] breseq: Identification of mutations in laboratory-evolved microbes from next-generation sequencing data using breseq.
                 Deatherage DE, Barrick JE (2014) 
                 Methods Mol Biol 1151:165-188.
```


## License
MIT License

