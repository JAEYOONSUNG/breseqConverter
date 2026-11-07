# breseqConverter

Parse and compare [breseq](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) reports, draw per-sample genome maps, and reconstruct sequences with the reported mutations applied.

**Start with [one_shot.R](https://raw.githubusercontent.com/JAEYOONSUNG/breseqConverter/master/inst/scripts/one_shot.R).** One file creates **both Excel formats**, one PDF per replicon, and changed-sequence GenBank files. Native SnapGene `.dna` files are also created when its CLI is available. Installing the breseqConverter package is optional.

## Features

- Extract mutation, missing-coverage and junction tables from breseq HTML reports.
- Compare samples in an explicit chronological order while preserving replicon identities.
- Write a complete wide Excel workbook **and** a separate compact publication workbook.
- Show every observed call with its stage color; use bold for new, changed or reappearing calls. Both workbooks use Arial 10 pt.
- Plot each sample's observed variants, with optional intergenic-label hiding and configurable GC colors.
- Apply finalized calls to the reference bases and export annotated GenBank plus conditional native SnapGene DNA files.

## Quick Start: one file

Requires **R 4.1 or newer**. No external R scripts or prepared GenomeDrawer objects are needed.

### 1. Install the R libraries once

```r
install.packages(c(
  "dplyr", "rvest", "openxlsx", "circlize", "RColorBrewer",
  "gtools", "reshape2", "stringr"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")

# For native .dna export through an installed SnapGene CLI:
install.packages("processx")
```

### 2. Download and run

Change `input_dir` to your input folder, then run this block in R or RStudio:

```r
download.file(
  "https://raw.githubusercontent.com/JAEYOONSUNG/breseqConverter/master/inst/scripts/one_shot.R",
  destfile = "one_shot.R", mode = "wb"
)

options(breseqConverter.one_shot.config = list(
  input_dir = "/path/to/breseq_outputs",
  sample_order = NULL,       # Or c("WT", "sample1", "sample2")
  sample_labels = NULL,      # Optional named display labels
  hide_intergenic_labels = TRUE,
  create_genbank = TRUE,
  create_dna = TRUE
))
source("one_shot.R")
```

The script detects the reference, GFF and eggNOG files and saves a new timestamped folder on Desktop. Inputs are read only. For HTML reports, `sample_order = NULL` uses natural filename order and adds an empty WT baseline when no WT report exists. Canonical Excel input keeps its stored or detected sample order. Set an explicit order for the intended construction history. If several different references or annotation files are found, set their paths explicitly.

Alternatively, download the same file, edit the settings at its top, and use Source in RStudio or run `Rscript --vanilla one_shot.R`.

### Output files

| Output in the results folder | Contents |
| --- | --- |
| `merged_mutation_results.xlsx` | Complete wide comparison with all source fields and hidden canonical-column metadata. |
| `publication_mutation_results.xlsx` | Organized publication tables: mutations, missing coverage, junction evidence, and supplied plot-derived evidence. |
| `plots/circos_<replicon>.pdf` | One genome plot per reference replicon; sample rings run from outside to inside in the configured order. |
| `changed_sequences/*_called_variants.gbk` | Per-sample, per-replicon sequences with finalized mutations applied and annotations lifted to the edited coordinates. |
| `changed_sequences/*_called_variants.dna` | Native SnapGene files after successful CLI conversion. |
| `changed_sequences/sequence_manifest.tsv` | Sequence lengths, output paths, replicon losses and DNA-conversion status. |
| `sample_call_counts.tsv` | Original call counts and plotted call counts, including any audited display additions. |
| `provenance/` | Source hashes, evidence tables and the R result object. |

The two Excel workbooks are always written by the one-shot workflow. Sequence files depend on `create_genbank`; DNA also requires `create_dna`, `processx` and an available SnapGene CLI. Lost replicons receive manifest entries instead of empty sequence files. SnapGene sessions that are already open are left running, with conversion status recorded in the manifest.

## One-shot configuration

Settings can be edited at the top of `one_shot.R`, or passed as lower-case keys through `options(breseqConverter.one_shot.config = list(...))` before sourcing it.

| Setting in the file | R option-list key | Purpose / default |
| --- | --- | --- |
| `INPUT_DIR` | `input_dir` | Input folder; defaults to the current working directory. |
| `OUTPUT_DIR` | `output_dir` | Defaults to a new timestamped Desktop folder; must be outside the input folder. |
| `REFERENCE_FILE` | `reference_file` | GenBank or FASTA reference; `NULL` detects it automatically. |
| `EXCEL_FILE` | `excel_file` | Optional canonical merged workbook input instead of HTML reports. |
| `GFF_FILE`, `EGGNOG_FILE` | `gff_file`, `eggnog_file` | Optional annotations; `NULL` detects them automatically. |
| `SAMPLE_ORDER` | `sample_order` | Exact sample IDs in chronological order; `NULL` uses natural filename order. |
| `SAMPLE_LABELS` | `sample_labels` | Named ID-to-display-label map, used in both Excel files and plots. |
| `BASELINE_ID` | `baseline_id` | For HTML input, `"WT"` adds an empty baseline when no WT report exists; `NULL` adds none. For Excel input, include a desired empty baseline in `sample_order`. |
| `HIDE_INTERGENIC_LABELS` | `hide_intergenic_labels` | Default `FALSE`; `TRUE` hides labels only, preserving marks and data. |
| `GC_COLORS` | `gc_colors` | Skew: green `#00FF00` / red `#FF0000`; ratio: blue `#0000FF` / yellow `#FFFF00`. |
| `CREATE_GENBANK`, `CREATE_DNA` | `create_genbank`, `create_dna` | Both default `TRUE`; DNA conversion follows successful changed-sequence export. |
| `SNAPGENE_PATH` | `snapgene_path` | Optional CLI path; `NULL` detects it automatically. |
| `DISPLAY_DELETIONS_FILE` | `display_deletions_file` | Optional audited coverage-evidence TSV; otherwise auto-detect `coverage_derived_DEL_rows.tsv`. |
| `SAMPLE_ALIASES` | `sample_aliases` | Optional mapping from evidence-file IDs to exact sample IDs. |
| `PUBLICATION_ANNOTATIONS_FILE` | `publication_annotation_overrides` | Optional TSV of publication-only annotation labels. |

HTML reports may be named `index_<sample>.html`, or stored as `<sample>/output/index.html`. Use the extracted IDs in `sample_order` and the names of `sample_labels`; display-label values can be different. Canonical Excel files produced by the converter are also supported. FASTA supplies bases but no gene features; use GenBank or a GFF when gene annotations are needed.

## One-shot Excel and sequence details

Every one-shot run writes both Excel files to the output folder:

- `merged_mutation_results.xlsx` retains the complete wide comparison, all
  reported fields and hidden canonical column metadata for reloading.
- `publication_mutation_results.xlsx` presents one row per variant or evidence
  region with compact stage columns. It separates mutations, missing coverage,
  paired junction evidence and, when supplied, plot-derived evidence. Gene and
  product annotations come from the supplied reference or GFF. Missing coverage
  wholly inside a previously called deletion is omitted from the compact view;
  every source row remains in the wide workbook.

Both use Arial 10 pt. Reported calls remain shaded in each observed stage;
new or changed calls use bold text. Unreported calls are never filled forward.
Use `SAMPLE_LABELS` to set the displayed strain names in both files. The returned
R result includes `excel_files`, with `comparison` and `publication` paths;
`excel_file` continues to refer to the wide workbook.

Optional publication-only labels can be supplied in
`PUBLICATION_ANNOTATIONS_FILE`, a tab-separated file with columns `sheet`,
`seq_id`, `position`, `field` and `value`. Use the reference accession and the
mutation position, a source start interval for missing coverage, or the displayed
start coordinate for plot-derived evidence. Only
annotation fields can be overridden; sample calls remain the original data.
The source comparison, plots and reconstructed sequences are unaffected.

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

---

## Optional R package installation

The one-shot workflow above does not require installing breseqConverter. For individual R API calls in the following sections:

```r
install.packages("remotes")
remotes::install_github("JAEYOONSUNG/breseqConverter")
library(breseqConverter)
```

`run_full_analysis()` and `run_breseq_pipeline()` retain their existing API behavior. Use the one-shot workflow to create both Excel workbooks, per-replicon plots and changed sequences automatically, or call the export helpers explicitly.

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
- `New in B` — Call reported in B and absent from A
- `Missing in B` — Call absent from B; this alone does not establish reversion
- `Unchanged or Modified` — Mutations present in both

### `compare_multiple_tables()`

Sequentially compares multiple strains for tracking mutations across evolution lineages.

```r
lineage <- list(WT = index_WT, R1 = index_R1, R2 = index_R2)
all_comparisons <- compare_multiple_tables(lineage)
```

---

## Circos Visualization

<img width="929" height="865" alt="image" src="https://github.com/user-attachments/assets/132092fe-64fe-4a5c-97d5-d675a02595dd" />

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
- **Legacy process check** — The reference-annotation API can prompt to close SnapGene for CLI conversion; the one-shot workflow does not close it.
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

## One-shot example: explicit stage order

After downloading `one_shot.R`, this runs extraction, both Excel exports, all replicon plots and changed-sequence export in one call:

```r
options(breseqConverter.one_shot.config = list(
  input_dir = "/path/to/breseq_outputs",
  sample_order = c("WT", "sample1", "sample2", "sample3"),
  sample_labels = c(WT = "Parent", sample1 = "Strain 1",
                    sample2 = "Strain 2", sample3 = "Final strain"),
  hide_intergenic_labels = TRUE,
  create_genbank = TRUE,
  create_dna = TRUE
))
source("one_shot.R")

ONE_SHOT_RESULTS$excel_files  # comparison and publication paths
ONE_SHOT_RESULTS$plots        # one PDF per replicon
ONE_SHOT_RESULTS$genome_export$manifest
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

### SnapGene process check in the legacy annotation exporter

The one-shot workflow leaves running SnapGene sessions untouched and records a skipped or pending conversion. The older `export_mutations_to_genbank()` API has a separate interactive process check when using `.dna` input:

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

### `save_publication_results_to_excel()`

For an existing native merged comparison, the compact file is also available
without running the plots or sequence export:

```r
save_publication_results_to_excel(
  merged_results, "publication_mutation_results.xlsx",
  sample_order = sample_ids,
  sample_labels = c(WT = "Parent", sample1 = "Strain 1",
                    sample2 = "Strain 2", sample3 = "Final strain"),
  reference_file = "reference.gbk"
)
```

### `get_screen_size()`

Helper function for adaptive plot sizing (macOS).

```r
dims <- get_screen_size(default_width = 12, default_height = 7, dpi = 96)
```

---

## Notes

- The one-shot workflow accepts breseq HTML reports or canonical merged Excel inputs.
- Comparisons retain reference replicon identities as well as coordinates; junction events use both oriented endpoints when paired keys are available.
- A missing call is not evidence that a previously deleted sequence or plasmid has returned. No call is carried forward without a report in that sample.
- Unassigned missing coverage and junctions remain evidence. Audited display additions are listed separately from finalized calls and reconstructed sequence edits.
- Intergenic-label hiding changes only labels, and does not establish that those variants have no functional effect.
- Excel sheet names are limited to 31 characters. Comma-formatted mutation lengths such as `Δ13,543 bp` are supported.

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
