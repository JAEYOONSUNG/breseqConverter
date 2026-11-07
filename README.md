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
# From GitHub
devtools::install_github("JAEYOONSUNG/breseqConverter")
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

Input files must be breseq HTML reports (`index.html`). Supply `sample_order`
when the intended chronological order differs from the file names.

```r
library(breseqConverter)

# Run complete analysis
run_full_analysis(
  input_dir = "./breseq_outputs",
  reference_file = "reference.gbff",      # Supports .gb, .gbff, .gbk
  eggnog_file = "eggnog_annotations.xlsx", # Optional, for COG colors
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

### `run_full_analysis()`

Complete pipeline for breseq analysis with Circos-style genome visualization.

```r
run_full_analysis(
  input_dir = "breseq_outputs",
  reference_file = "reference.gbk",
  eggnog_file = NULL,
  output_dir = "breseq_analysis",
  add_genome_map = TRUE
)
```

### Intergenic labels and publication Excel output

These options are available in `run_full_analysis()` and
`run_breseq_pipeline()`. Their defaults preserve the existing output.

```r
sample_ids <- c("WT", "sample1", "sample2", "sample3")

results <- run_full_analysis(
  input_dir = "breseq_outputs",
  output_dir = "publication_output",
  reference_file = "reference.gbk",
  add_genome_map = FALSE,
  create_genbank = FALSE,
  hide_intergenic_labels = TRUE,
  highlight_new = TRUE,
  clean_headers = TRUE,
  sample_order = sample_ids,
  sample_labels = c(
    WT = "Parent", sample1 = "Strain 1", sample2 = "Strain 2",
    sample3 = "Final strain"
  )
)
```

Use the exact extracted sample IDs in `sample_order` and the names of
`sample_labels`; the values of `sample_labels` are display names.

- `hide_intergenic_labels = TRUE` hides labels whose variant annotation starts
  with `intergenic`. Mutation marks and the underlying data remain present.
  The option is also available directly in `plot_circos_genome()`.
- `highlight_new = TRUE` shades every reported call with its stage's pale color,
  including retained calls, and uses bold text for new or changed calls.
  Evidence read counts and annotation changes alone do not create a new variant.
- `clean_headers = TRUE` replaces technical Excel column names with readable
  field and strain names. Canonical column names are retained in hidden workbook
  metadata so the package can reload the exported workbook.
- `sample_order` specifies the comparison order. `sample_labels` changes the
  displayed strain names without changing sample IDs in the returned R data.

The Excel export helpers accept the same four Excel options. Missing coverage
and junction evidence remain evidence: a missing call does not demonstrate
sequence restoration, and a shifted coverage boundary is not automatically a
new deletion. Intergenic annotations alone do not establish that a variant has
no functional effect.

Opt-in Excel exports use Arial 10 pt and restrained stage colors. Calls remain
shaded in every sample where they were reported. Bold text marks a new,
changed, or reappearing call; revised missing-coverage boundaries and evidence
inside a prior deletion stay in regular weight. An unreported call is unfilled,
and earlier calls are never carried forward into samples without evidence.
Variant classification uses the first reported allele at each locus; the sample
columns retain every reported allele.
Junction highlighting uses paired endpoint keys from native comparisons. To
highlight junctions after importing a plain Excel table, rerun the HTML
comparison; column-header metadata alone does not encode the endpoint pairs.

### Circos Plot Tracks (outer → inner)

1. **Chromosome ideogram** — Genome coordinates with Mb labels
2. **Gene annotations** — COG-colored genes (requires eggNOG file)
3. **RNA features** — tRNA (green), rRNA (red)
4. **GC skew** — At or above the window mean: green `#00FF00`; below: red `#FF0000`
5. **GC ratio** — At or above the window mean: blue `#0000FF`; below: yellow `#FFFF00`
6. **Mutation tracks** — One track per sample

GC skew uses `(C - G) / (C + G)`. The full workflow plots each statistic's
deviation from its window mean. Colors are assigned from the same values used
for each track, so cached color columns cannot swap the palettes. To change the
colors, pass a named vector to `run_full_analysis(gc_colors = ...)`:

```r
gc_colors <- c(
  skew_positive = "#00FF00", skew_negative = "#FF0000",
  ratio_positive = "#0000FF", ratio_negative = "#FFFF00"
)
```

The four values above match the legacy GenomeDrawer palette. Partial overrides
are accepted. The intergenic label setting is independent of these colors.

### Standalone script

Download [`one_shot.R`](https://raw.githubusercontent.com/JAEYOONSUNG/breseqConverter/master/inst/scripts/one_shot.R)
for a self-contained entry point. It includes the converter implementation and
automatically reads the inputs, prepares genome annotations, exports the
workbook, plots every replicon, and reconstructs changed GenBank sequences.
When SnapGene's CLI is available it also creates native `.dna` files. It requires the
listed R libraries, but no breseqConverter installation, external R scripts,
GenomeDrawer objects, or separate drawing call.

Edit the paths and settings at the top, then use Source in RStudio or run:

```sh
Rscript --vanilla one_shot.R
```

Set `INPUT_DIR` to your data folder. The defaults detect sample IDs, add an
empty `WT` reference baseline, and write a new output folder on Desktop. Supply
`SAMPLE_ORDER` when natural filename order differs from the intended chronology.
Two display settings are:

```r
HIDE_INTERGENIC_LABELS <- TRUE  # FALSE displays the labels
GC_COLORS <- c(
  skew_positive = "#00FF00", skew_negative = "#FF0000",
  ratio_positive = "#0000FF", ratio_negative = "#FFFF00"
)
```

Each ring contains its sample's observed calls, including retained mutations.
The configured sample order runs from outside to inside. Reference GenBank or
FASTA, GFF and eggNOG inputs can be detected automatically or supplied explicitly.
Canonical merged Excel inputs are also supported. Audited coverage display
evidence is loaded from `coverage_derived_DEL_rows.tsv` when available, kept
separate from the original calls, and copied into output provenance. The
converter does not infer missing coverage values from unavailable summary files.

For other datasets, change `INPUT_DIR` and set `SAMPLE_ORDER`, `SAMPLE_LABELS`,
and `SAMPLE_ALIASES` to `NULL` for automatic sample IDs. Set `BASELINE_ID` to
`NULL` when a reference baseline should not be added. `one_shot_legacy.R`
preserves the former script that requires prepared GenomeDrawer objects.

Sequence export settings are `CREATE_GENBANK = TRUE` and `CREATE_DNA = TRUE`.
Finalized SNPs, insertions, deletions, replacements and exact tandem-repeat
changes are applied to each sample's reference sequence. Features are lifted
to the edited coordinates, and affected coding annotations are reviewed or
retranslated. Invalid reference alleles, unsupported calls and overlapping
edits stop export with an edit report rather than silently making a sequence.

Each retained sample/replicon produces a `*_called_variants.gbk` under
`changed_sequences/`, plus `.dna` after successful SnapGene conversion.
`sequence_manifest.tsv` records exported files, lengths, losses and conversion
statuses. A finalized whole-replicon deletion produces a loss record instead
of an empty DNA file. Audited whole-replicon coverage loss is recorded
separately as `lost_by_coverage`. Unassigned MC boundaries are not exact
deletion breakpoints: these intervals go into `unresolved_coverage_regions.tsv`
and are not silently removed from the reconstructed sequence.

SnapGene must be installed and not already running for its documented CLI.
If it is open, unavailable or delayed, the GenBank files remain usable and the
manifest records the conversion status; the workflow never quits the app.

Developers can regenerate the entry point with `Rscript tools/build_standalone.R`.

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

### Changed-sequence export

`export_changed_sequences_to_genbank()` reconstructs each sample/replicon from
its finalized calls. The standalone workflow calls this automatically:

```r
export_changed_sequences_to_genbank(
  reference_file = "reference.gbk",
  mutation_list = raw_calls_by_sample,
  output_dir = "changed_sequences"
)
```

### Reference-coordinate annotation export (legacy)

`export_mutations_to_genbank()` and `create_annotated_genbank()` retain the
reference bases and add colored mutation annotations at reference coordinates.
They do not reconstruct the mutated sequence. Use the changed-sequence exporter
or the standalone script for genomes with the calls applied to their bases.

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

# Place each sample's breseq report under this input directory.
run_full_analysis(
  input_dir = "./breseq_output",
  reference_file = "GCB.gbk",
  eggnog_file = "GCB_eggnog.tsv",
  output_dir = "breseq_analysis",
  add_genome_map = TRUE,
  sample_order = c("WT", "45SJY3", "45SJY5", "45SJY7")
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
| `New in B` | Call reported in B and absent from A |
| `Missing in B` | Call reported in A and absent from B; does not by itself prove reversion |
| `Unchanged or Modified` | Locus reported in both samples; compare the alleles and evidence |

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
