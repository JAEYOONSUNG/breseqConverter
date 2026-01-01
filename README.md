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


## SnapGene Visualization
![45DF Map](https://github.com/user-attachments/assets/048184e4-4ed5-45e0-9fe4-3281383ce3d1)

### `breseq_for_snapgene.py`

Converts breseq mutation data to GenBank format for visualization in SnapGene. Mutations are annotated as `misc_feature` entries with color-coded labels indicating the stage of first detection.

### Usage

```bash
python breseq_for_snapgene.py \
  --dna_file reference.dna \
  --mutation_xlsx merged_results.xlsx \
  --output_dir output/
```

### Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `--dna_file` | Reference genome (.dna or .gbk) | Required |
| `--mutation_xlsx` | Excel file from R comparison | `merged_mutant_results.xlsx` |
| `--output_dir` | Output directory | `deletion_results/` |

### Features

**1. Automatic Format Conversion**
- Converts SnapGene `.dna` files to GenBank format using SnapGene CLI
- Preserves existing annotations from the reference genome

**2. Mutation Annotation**
- Parses mutation positions and calculates affected ranges
- Supports SNPs, insertions (`+`), deletions (`Δ`), and missing coverage regions
- Handles complex position formats (e.g., `1234-5678` for uncertain boundaries)

**3. Color-Coded Visualization**
- Mutations are color-coded by detection stage (evolution round)
- Uses HSV color gradient: earlier stages = more saturated magenta
- Enables quick visual identification of when mutations arose

**4. GenBank Feature Format**
```
     misc_feature    12345..12350
                     /label="[Mutation] Δ6 bp geneX"
                     /note="First detected in Round2"
                     /note="color: #FF00FF"
```

### Output

The script generates an annotated GenBank file (`.gbk`) that can be:
- Opened directly in SnapGene
- Imported into other genome viewers (Benchling, ApE, etc.)
- Used for publication-quality genome maps

### Example Workflow

```bash
# Step 1: Run breseq on each strain
breseq -r reference.gbk reads_WT.fastq -o breseq_WT/
breseq -r reference.gbk reads_R1.fastq -o breseq_R1/
breseq -r reference.gbk reads_R2.fastq -o breseq_R2/

# Step 2: Extract and compare in R
Rscript -e "
library(breseqConverter)
extract_named_tables('breseq_output/')
strains <- list(WT=index_WT, R1=index_R1, R2=index_R2)
results <- compare_multiple_tables(strains)
save_comparison_results_to_excel(results[[1]], 'merged_results.xlsx')
"

# Step 3: Generate SnapGene-compatible annotation
python breseq_for_snapgene.py \
  --dna_file reference.dna \
  --mutation_xlsx merged_results.xlsx \
  --output_dir snapgene_output/

# Step 4: Open in SnapGene
open snapgene_output/annotated_dna_with_mutations.gbk
```

---

## Typical Workflows

### Workflow 1: Strain Engineering Tracking

Track accumulated mutations across sequential engineering rounds.

```r
# Extract breseq outputs
extract_named_tables("./breseq_results/")

# Organize by engineering stage
strains <- list(
  WT = index_WT,
  KO_round1 = index_R1,
  KO_round2 = index_R2,
  final = index_final
)

# Compare sequentially
comparisons <- compare_multiple_tables(strains)

# Export for SnapGene visualization
save_comparison_results_to_excel(comparisons$Comparison_1_to_2, "R1_to_R2.xlsx")
```

### Workflow 2: Parallel Evolution Analysis

Compare independently evolved lineages to identify convergent mutations.

```r
ancestor <- index_ancestor

# Compare each lineage to ancestor
for (name in c("lineage_A", "lineage_B", "lineage_C")) {
  comparison <- compare_named_tables(ancestor, get(paste0("index_", name)))
  save_comparison_results_to_excel(comparison, paste0("vs_ancestor_", name, ".xlsx"))
}
```

---

## Output Interpretation

### Mutation Status Guide

| Status | Biological Interpretation |
|--------|---------------------------|
| `New in B` | Acquired mutation (likely beneficial or hitchhiker) |
| `Missing in B` | Reversion or secondary suppressor |
| `Unchanged` | Stable mutation maintained across rounds |

### SnapGene Color Coding

| Color Saturation | Meaning |
|------------------|---------|
| Bright magenta | Early-stage mutation (Round 1) |
| Pale magenta | Later-stage mutation (Round 3+) |
| Default magenta | Unknown detection stage |



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
