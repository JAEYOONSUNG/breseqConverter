# breseqConverter

An R package for parsing, comparing, and visualizing genomic mutation data from [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) output files. Designed for tracking genomic changes during bacterial strain engineering and adaptive laboratory evolution (ALE) experiments.

<img width="929" height="865" alt="image" src="https://github.com/user-attachments/assets/132092fe-64fe-4a5c-97d5-d675a02595dd" />


## Features

- **Extract** mutation tables from breseq HTML reports
- **Compare** genomic changes between strains or evolution time points
- **Track** sequential mutations across multiple rounds of strain improvement
- **Visualize** with Circos-style genome maps (COG annotations, GC content, RNA features)
- **Export** to SnapGene-compatible GenBank format with color-coded mutations

## Installation

```r
# From GitHub (recommended)
devtools::install_github("JAEYOONSUNG/breseqConverter")

# Or install locally
install.packages("path/to/breseqConverter", repos = NULL, type = "source")
```

## Dependencies

```r
# Required
install.packages(c("rvest", "dplyr", "openxlsx"))

# For Circos visualization
install.packages(c("circlize", "RColorBrewer", "gtools"))
BiocManager::install(c("ComplexHeatmap", "Biostrings"))
```

---

## Quick Start

### Full Pipeline with Circos Visualization

```r
library(breseqConverter)

# Define breseq output directories
breseq_dirs <- c(
  "WT" = "./breseq_WT",
  "45SJY3" = "./breseq_45SJY3",
  "45SJY5" = "./breseq_45SJY5",
  "45SJY7" = "./breseq_45SJY7"
)

# Run complete analysis
run_breseq_analysis(
  breseq_dirs = breseq_dirs,
  reference_gbk = "reference.gbk",
  eggnog_file = "eggnog_annotations.tsv",  # optional for COG colors
  output_dir = "breseq_analysis",
  add_genome_map = TRUE
)

# Export to SnapGene
export_mutations_to_genbank(
  reference_file = "reference.dna"
  # Auto-detects breseq_analysis/merged_mutation_results.xlsx
)
```

**Output files:**
- `breseq_analysis/merged_mutation_results.xlsx` — Combined mutation data
- `breseq_analysis/circos_plot.pdf` — Circos genome visualization
- `breseq_analysis/<reference>_with_mutations.gbk` — SnapGene-compatible GenBank

---

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

### `compare_named_tables()`

Compares mutation profiles between two strains.

```r
comparison <- compare_named_tables(strain_A, strain_B)
```

**Status values:**
- `New in B` — Mutations acquired in strain B
- `Missing in B` — Mutations lost in strain B
- `Unchanged or Modified` — Mutations present in both

### `compare_multiple_tables()`

Sequentially compares multiple strains for tracking mutations across evolution lineages.

```r
lineage <- list(WT = index_WT, R1 = index_R1, R2 = index_R2)
all_comparisons <- compare_multiple_tables(lineage)
```

---

## Circos Visualization

### `run_breseq_analysis()`

Complete pipeline for breseq analysis with Circos-style genome visualization.

```r
run_breseq_analysis(
  breseq_dirs,
  reference_gbk,
  eggnog_file = NULL,
  output_dir = "breseq_analysis",
  add_genome_map = TRUE
)
```

### Circos Plot Tracks (outer → inner)

1. **Chromosome ideogram** — Genome coordinates with Mb labels
2. **Gene annotations** — COG-colored genes (requires eggNOG file)
3. **RNA features** — tRNA (green), rRNA (red)
4. **GC skew** — Positive (blue) / Negative (yellow)
5. **GC ratio** — Above average (green) / Below average (red)
6. **Mutation tracks** — One track per sample

### Mutation Color Scheme

| Type | Color | Hex Code |
|------|-------|----------|
| SNP | Teal | `#139487` |
| INS (Insertion) | Pink | `#FF5C8D` |
| DEL (Deletion) | Dark Navy | `#03045E` |
| REP (Replacement) | Golden Yellow | `#E5C87B` |

### Missing Coverage (MC) Types

| Type | Description |
|------|-------------|
| `MC_repeat` | Repeat regions where reads map to multiple locations (IS elements, rRNA) |
| `MC_inversion` | Regions where reads map in opposite orientation (potential inversions) |

### COG Category Colors

Standard COG functional category colors are used when eggNOG annotations are provided.

---

## SnapGene Export

### `export_mutations_to_genbank()`

Exports breseq mutations to SnapGene-compatible GenBank format.

```r
export_mutations_to_genbank(
  reference_file,
  mutation_xlsx = NULL,      # Auto-detects from breseq_analysis/
  output_file = NULL,        # Auto-generates in breseq_analysis/
  sheets = c("Mutations", "Missing"),
  color_by = c("type", "stage")
)
```

### Features

- **Auto-detection** — Finds `breseq_analysis/merged_mutation_results.xlsx` automatically
- **SnapGene .dna support** — Converts .dna files via SnapGene CLI
- **Process check** — Prompts to close SnapGene if running (required for CLI conversion)
- **CONTIG removal** — Fixes SnapGene compatibility issues
- **ApE color format** — Colors visible in SnapGene

### Label Format

```
[SNP:45SJY3] A→G cyaA
[DEL:45SJY5] Δ13,543 bp hypothetical protein
[MC:45SJY7] (5,234 bp) intergenic region
```

- Mutation type and first detection stage in brackets
- Full mutation description with commas for readability
- Gene/feature annotation

### GenBank Feature Output

```
     misc_feature    12345..25678
                     /label=[DEL:45SJY3] Δ13,334 bp geneX
                     /ApEinfo_fwdcolor="#03045E"
                     /ApEinfo_revcolor="#03045E"
                     /note="color: #03045E"
```

### Color Options

- `color_by = "type"` — Color by mutation type (SNP, DEL, INS, REP)
- `color_by = "stage"` — Magenta gradient by detection stage (saturated → pale)

### SnapGene Visualization

<img width="792" height="768" alt="스크린샷 2026-02-06 21 02 41" src="https://github.com/user-attachments/assets/e094a335-ec2c-417b-b0e0-068eaeddd874" />

---

## Complete Workflow Example

### Tracking Mutations Across ALE Rounds

```r
library(breseqConverter)

# Step 1: Define breseq directories
breseq_dirs <- c(
  "WT" = "./breseq_output/WT",
  "45SJY3" = "./breseq_output/45SJY3",
  "45SJY5" = "./breseq_output/45SJY5",
  "45SJY7" = "./breseq_output/45SJY7"
)

# Step 2: Run full analysis with Circos plot
run_breseq_analysis(
  breseq_dirs = breseq_dirs,
  reference_gbk = "GCB.gbk",
  eggnog_file = "GCB_eggnog.tsv",
  output_dir = "breseq_analysis",
  add_genome_map = TRUE
)

# Step 3: Export to SnapGene
export_mutations_to_genbank(
  reference_file = "GCB.dna",
  color_by = "type"
)

# Generated files:
# - breseq_analysis/merged_mutation_results.xlsx
# - breseq_analysis/circos_plot.pdf
# - breseq_analysis/GCB_with_mutations.gbk
```

### Step-by-Step Usage

```r
# 1. Extract tables from breseq HTML output
extract_named_tables(directory = "path/to/breseq/output")

# 2. Compare two strains
comparison <- compare_named_tables(strain_A, strain_B)

# 3. Save results to Excel
save_comparison_results_to_excel(comparison, "comparison_results.xlsx")

# 4. Export to SnapGene
export_mutations_to_genbank(
  reference_file = "reference.gbk",
  mutation_xlsx = "comparison_results.xlsx"
)
```

---

## SnapGene Integration

### Automatic SnapGene Detection

```r
# Auto-detect SnapGene path
find_snapgene()

# Or set manually
options(breseqConverter.snapgene_path = "/custom/path/to/snapgene")
```

**Supported paths:**
- **macOS**: `/Applications/SnapGene.app/Contents/MacOS/SnapGene`
- **Windows**: `C:\Program Files\SnapGene\snapgene.exe`
- **Linux**: `/usr/bin/snapgene`, `/opt/snapgene/snapgene`

### SnapGene Process Check

When using `export_mutations_to_genbank()` with .dna files, the function checks if SnapGene is running:

```
⚠️  SnapGene is active. Close application to proceed? [Y/n]:
```

- Press Enter or `Y` to close SnapGene and proceed
- Press `n` to continue without closing

---

## Output Interpretation

### Mutation Status

| Status | Biological Interpretation |
|--------|---------------------------|
| `New in B` | Acquired mutation (beneficial or hitchhiker) |
| `Missing in B` | Reversion or secondary suppressor |
| `Unchanged` | Stable mutation maintained across rounds |

### SnapGene Label Guide

| Label Format | Meaning |
|--------------|---------|
| `[SNP:45SJY3]` | SNP first detected in sample 45SJY3 |
| `[DEL:45SJY5]` | Deletion first detected in sample 45SJY5 |
| `[MC:45SJY7]` | Missing coverage first detected in sample 45SJY7 |
| `Δ13,543 bp` | Deletion size with comma formatting |

---

## Utility Functions

### `save_comparison_results_to_excel()`

Exports comparison results to a multi-sheet Excel workbook.

```r
save_comparison_results_to_excel(comparison, "output.xlsx")
```

### `get_screen_size()`

Helper function for adaptive plot sizing (macOS).

```r
dims <- get_screen_size(default_width = 12, default_height = 7, dpi = 96)
```

---

## Notes

- Input files must be breseq HTML output (`index.html`)
- Comparison is based on genomic `position` column
- For Missing Coverage, comparison uses `start` column
- Excel sheet names are truncated to 31 characters (Excel limitation)
- SnapGene must be closed when converting .dna files
- Numbers with commas (e.g., `Δ13,543 bp`) are parsed correctly

---

## Citation

If you use this package, please cite:

```
[DNMB] DNMB: Deletion of non-canonical defense systems enable domestication of
       non-model Geobacillus for thermophilic platform engineering.
       Jae-Yoon Sung, Mun Hoe Lee, Jungsoo Park, Hyungbin Kim, Dariimaa Ganbat,
       Donggyu Kim, Hyun-Woo Cho, Sang Jae Lee, Seong Bo Kim, and Dong-Woo Lee. 2024.
       XXX, XXX, https://doi.org/XXX

[breseq] Identification of mutations in laboratory-evolved microbes from
         next-generation sequencing data using breseq.
         Deatherage DE, Barrick JE (2014) Methods Mol Biol 1151:165-188.
```

## License

MIT License
