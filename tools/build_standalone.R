# Development-only builder. The generated file never sources package files.
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args)) normalizePath(args[1L]) else normalizePath(".")
target <- file.path(repo, "inst", "scripts", "one_shot.R")
legacy <- file.path(repo, "inst", "scripts", "one_shot_legacy.R")
if (!file.exists(legacy) && file.exists(target)) file.copy(target, legacy)
header <- c(
  "# breseqConverter: one file, automatic extraction and plotting (2026-11-07)",
  "# RStudio: edit the settings below, then Source. Terminal: Rscript --vanilla one_shot.R",
  "# No breseqConverter installation, external R scripts, or GenomeDrawer globals required.",
  "# Required R libraries are checked before reading data. Original inputs are read only.",
  "# Set INPUT_DIR to your data folder. New results are saved outside the inputs.",
  "# Excel outputs: merged_mutation_results.xlsx and publication_mutation_results.xlsx.",
  "INPUT_DIR <- getwd()  # Or set an absolute path to your breseq input folder",
  "OUTPUT_DIR <- file.path(path.expand('~/Desktop'), paste0('breseqConverter_run_', format(Sys.time(), '%Y%m%d_%H%M%S')))",
  "REFERENCE_FILE <- NULL  # NULL: find the reference GenBank/FASTA automatically",
  "EXCEL_FILE <- NULL      # NULL: read breseq HTML; or supply a canonical exported workbook",
  "GFF_FILE <- NULL        # NULL: auto-detect gene annotations",
  "EGGNOG_FILE <- NULL     # NULL: auto-detect COG annotations",
  "HIDE_INTERGENIC_LABELS <- FALSE  # TRUE hides labels only; marks remain",
  "SAMPLE_ORDER <- NULL  # Or c('WT', 'sample1', 'sample2') in chronological order",
  "SAMPLE_LABELS <- NULL # Optional named map, e.g. c(WT='Parent', sample1='Strain 1')",
  "BASELINE_ID <- 'WT'  # Empty reference baseline when WT was not resequenced",
  "SAMPLE_ALIASES <- NULL # Optional map from audited-evidence IDs to input sample IDs",
  "DISPLAY_DELETIONS_FILE <- NULL  # Auto-find the existing audited coverage evidence TSV",
  "PUBLICATION_ANNOTATIONS_FILE <- NULL # Optional TSV of publication-only annotation labels",
  "GC_COLORS <- c(skew_positive='#00FF00', skew_negative='#FF0000', ratio_positive='#0000FF', ratio_negative='#FFFF00')",
  "CREATE_GENBANK <- TRUE # Apply called mutations to the sequence, not just annotations",
  "CREATE_DNA <- TRUE     # Convert changed GenBank files when SnapGene CLI is available",
  "SNAPGENE_PATH <- NULL  # NULL: auto-detect; running SnapGene sessions are never closed",
  "# Everything below is bundled implementation; no other code file is loaded.",
  ".breseq_one_shot_config <- list(input_dir=INPUT_DIR, output_dir=OUTPUT_DIR, reference_file=REFERENCE_FILE, excel_file=EXCEL_FILE, gff_file=GFF_FILE, eggnog_file=EGGNOG_FILE, hide_intergenic_labels=HIDE_INTERGENIC_LABELS, sample_order=SAMPLE_ORDER, sample_labels=SAMPLE_LABELS, baseline_id=BASELINE_ID, sample_aliases=SAMPLE_ALIASES, display_deletions_file=DISPLAY_DELETIONS_FILE, publication_annotation_overrides=PUBLICATION_ANNOTATIONS_FILE, gc_colors=GC_COLORS, create_genbank=CREATE_GENBANK, create_dna=CREATE_DNA, snapgene_path=SNAPGENE_PATH)",
  "breseq_one_shot <- local({"
)
sources <- sort(list.files(file.path(repo, "R"), pattern = "[.]R$", full.names = TRUE))
body <- unlist(lapply(sources, function(file) c(paste0("# Bundled: ", basename(file)),
  sub("[[:blank:]]+$", "", readLines(file, warn = FALSE)), "")))
footer <- c(
  "  function(config) .breseq_run_standalone(config)",
  "})",
  "if (!identical(getOption('breseqConverter.autorun', TRUE), FALSE)) {",
  "  .breseq_overrides <- getOption('breseqConverter.one_shot.config', list())",
  "  ONE_SHOT_RESULTS <- breseq_one_shot(utils::modifyList(.breseq_one_shot_config, .breseq_overrides, keep.null=TRUE))",
  "}"
)
writeLines(c(header, body, footer), target, useBytes = TRUE)
invisible(parse(target))
message("Built self-contained single file: ", target)
