# Internal driver embedded into the generated single-file entry point.
.breseq_standalone_pick <- function(input_dir, pattern, explicit = NULL, optional = FALSE) {
  if (!is.null(explicit)) {
    if (!file.exists(explicit)) stop("Input file does not exist: ", explicit)
    return(normalizePath(explicit))
  }
  candidates <- list.files(input_dir, pattern = pattern, recursive = TRUE,
                           full.names = TRUE, ignore.case = TRUE)
  candidates <- candidates[!grepl("/provenance/|/\\.git/", candidates)]
  if (!length(candidates)) {
    if (optional) return(NULL)
    stop("No matching input file found in ", input_dir, ": ", pattern)
  }
  if (length(candidates) > 1L) {
    hashes <- unname(tools::md5sum(candidates))
    if (length(unique(hashes)) != 1L)
      stop("Multiple different input files match; set the corresponding path explicitly:\n",
           paste(candidates, collapse = "\n"))
    candidates <- candidates[order(nchar(candidates))]
  }
  normalizePath(candidates[1L])
}

.breseq_standalone_sample_id <- function(file) {
  stem <- tools::file_path_sans_ext(basename(file))
  if (tolower(stem) != "index")
    return(sub("_index$", "", sub("^index_", "", stem, ignore.case = TRUE), ignore.case = TRUE))
  parent <- dirname(file)
  if (basename(parent) == "output") parent <- dirname(parent)
  basename(parent)
}

.breseq_standalone_unresolved_coverage <- function(data, samples) {
  empty <- data.frame(sample = character(), seq_id = character(),
    source_start = character(), source_end = character(), gene = character(),
    interpretation = character(), stringsAsFactors = FALSE)
  if (is.null(data) || !nrow(data)) return(empty)
  columns <- .breseq_stage_columns(data, samples)
  seqid <- .breseq_wide_seq(data, columns)
  records <- list()
  for (sample in samples) {
    ends <- .breseq_field(data, columns, "end", sample)
    gene <- .breseq_field(data, columns, "gene", sample)
    present <- which(.breseq_present(ends))
    if (length(present)) records[[length(records) + 1L]] <- data.frame(
      sample = sample, seq_id = seqid[present], source_start = as.character(data$start[present]),
      source_end = ends[present], gene = gene[present],
      interpretation = "Unassigned missing coverage: exact deletion breakpoints are not established; not applied to reconstructed sequence.",
      stringsAsFactors = FALSE)
  }
  if (length(records)) do.call(rbind, records) else empty
}

.breseq_run_standalone <- function(config) {
  required <- c("dplyr", "rvest", "openxlsx", "circlize", "ComplexHeatmap",
                "RColorBrewer", "gtools", "reshape2", "stringr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". Install these packages, then run this same file again.")
  suppressPackageStartupMessages(library(dplyr))
  if (getRversion() < "4.1.0") stop("This script requires R 4.1 or newer.")
  input <- normalizePath(path.expand(config$input_dir), mustWork = TRUE)
  output <- path.expand(config$output_dir)
  if (!grepl("^/", output)) output <- file.path(getwd(), output)
  check_output <- file.path(normalizePath(dirname(output), mustWork = FALSE), basename(output))
  if (check_output == input || startsWith(check_output, paste0(input, "/")))
    stop("output_dir must be outside input_dir so original inputs remain unchanged.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  if (output == input || startsWith(output, paste0(input, "/")))
    stop("output_dir must be outside input_dir so original inputs remain unchanged.")
  message("Input: ", input, "\nOutput: ", output)
  work <- file.path(output, "provenance")
  dir.create(work, showWarnings = FALSE)

  reference <- .breseq_standalone_pick(input,
    "\\.(gbff|gbk|gb|genbank|fasta|fna|fa)$", config$reference_file)
  refs <- .breseq_prepare_standalone_reference(reference, file.path(work, "reference"))
  gff <- .breseq_standalone_pick(input, "\\.gff3?$", config$gff_file, optional = TRUE)
  eggnog <- .breseq_standalone_pick(input, "emapper[.]annotations([.](xlsx|tsv))?$",
                                   config$eggnog_file, optional = TRUE)
  gff_data <- if (is.null(gff)) NULL else parse_gff_file(gff)
  eggnog_data <- if (is.null(eggnog)) NULL else parse_eggnog_file(eggnog)
  html <- list.files(input, pattern = "index.*[.]html$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
  html <- html[!grepl("/provenance/|/\\.git/", html)]
  source_files <- c(reference, gff, eggnog)
  baseline <- config$baseline_id
  order <- config$sample_order
  if (is.null(config$excel_file) && length(html)) {
    ids <- vapply(html, .breseq_standalone_sample_id, character(1))
    if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("HTML sample IDs are empty or duplicated.")
    tables <- setNames(lapply(html, extract_single_html), ids)
    if (is.null(order)) {
      order <- ids[gtools::mixedorder(ids)]
      if (!is.null(baseline) && !baseline %in% order) order <- c(baseline, order)
    }
    if (anyNA(order) || anyDuplicated(order) || !all(ids %in% order) ||
        any(!order %in% c(ids, baseline)))
      stop("sample_order must contain each HTML sample once and, optionally, baseline_id.")
    if (!is.null(baseline) && baseline %in% order && !baseline %in% ids)
      tables[[baseline]] <- lapply(tables[[1L]], function(d) d[FALSE, , drop = FALSE])
    tables <- tables[order]
    if (length(tables) < 2L) stop("A single sample requires a distinct baseline_id.")
    merged <- merge_comparison_results(compare_multiple_tables(tables))
    mut_name <- names(merged)[vapply(names(merged), function(n) .breseq_table_kind(n) == "mutation", logical(1))]
    if (length(mut_name)) merged[[mut_name[1L]]] <- add_mutation_types(merged[[mut_name[1L]]])
    mutation_data <- setNames(merged, c("Mutations", "Missing", "New junction")[match(names(merged),
      c("Mutation predictions", "Unassigned missing coverage evidence", "Unassigned new junction evidence"))])
    source_files <- c(source_files, html)
  } else {
    excel <- .breseq_standalone_pick(input, "(merged_mutation_results|mutant_results)[.]xlsx$",
                                    config$excel_file)
    mutation_data <- load_mutation_excel(excel)
    if (is.null(mutation_data$Mutations)) stop("No canonical Mutations worksheet in ", excel)
    if (is.null(order)) {
      order <- attr(mutation_data$Mutations, "breseq_sample_order", exact = TRUE)
      if (is.null(order)) order <- attr(.breseq_stage_columns(mutation_data$Mutations), "sample_order")
    }
    merged <- setNames(mutation_data,
      c("Mutation predictions", "Unassigned missing coverage evidence", "Unassigned new junction evidence")[match(names(mutation_data), c("Mutations", "Missing", "New junction"))])
    source_files <- c(source_files, excel)
  }
  labels <- .breseq_validate_labels(config$sample_labels, order)
  excel_out <- file.path(output, "merged_mutation_results.xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
  sheet_map <- c("Mutation predictions" = "Mutations",
    "Unassigned missing coverage evidence" = "Missing",
    "Unassigned new junction evidence" = "New junction")
  for (name in names(merged)) {
    sheet <- unname(sheet_map[name])
    if (is.na(sheet)) sheet <- substr(name, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    highlight <- .breseq_table_kind(name) != "junction" ||
      !is.null(attr(merged[[name]], "breseq_keys")) || "comparison_key" %in% names(merged[[name]])
    if (!highlight) message("Junction pair identities are unavailable in the imported workbook; preserving its values without inferring stage event highlights.")
    .breseq_write_presentation(wb, sheet, merged[[name]], name,
      highlight_new = highlight, clean_headers = TRUE, sample_labels = labels,
      sample_order = order, mutation_data = mutation_data$Mutations)
  }
  .breseq_save_table_workbook(wb, excel_out)
  processed <- process_mutation_data_for_viz(mutation_data, sample_order = order)
  # A workbook may retain region annotations in a sample with no variant call.
  # Annotation/replicon metadata alone must not become an UNKNOWN mutation.
  calls <- lapply(processed$mutation_list, function(d) {
    d[.breseq_present(.breseq_text(d$mutation)), , drop = FALSE]
  })
  raw_calls <- calls
  raw_counts <- setNames(vapply(order, function(s) if (is.null(calls[[s]])) 0L else nrow(calls[[s]]), integer(1)), order)

  display_file <- config$display_deletions_file
  if (is.null(display_file)) {
    candidates <- list.files(input, pattern = "^coverage_derived_DEL_rows[.]tsv$",
                             recursive = TRUE, full.names = TRUE)
    if (length(candidates) > 1L && length(unique(unname(tools::md5sum(candidates)))) != 1L)
      stop("Several different coverage display tables exist; set display_deletions_file.")
    if (length(candidates)) display_file <- candidates[order(nchar(candidates))][1L]
  }
  display_rows <- NULL
  if (!is.null(display_file)) {
    aliases <- config$sample_aliases
    if (is.null(aliases)) {
      short_ids <- sub("^[0-9]+(SJY)?", "", names(calls), ignore.case = TRUE, perl = TRUE)
      if (anyDuplicated(short_ids))
        stop("Coverage evidence IDs are ambiguous; set SAMPLE_ALIASES explicitly.")
      aliases <- stats::setNames(names(calls), short_ids)
    }
    display <- .breseq_standalone_display_deletions(calls, display_file,
      aliases, refs$info, merged_missing = mutation_data$Missing)
    calls <- display$mutation_list
    display_rows <- display$display_rows
    utils::write.table(display_rows, file.path(work, "display_deletion_evidence.tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    source_files <- c(source_files, display_file)
  }
  plot_counts <- setNames(vapply(order, function(s) if (is.null(calls[[s]])) 0L else nrow(calls[[s]]), integer(1)), order)
  counts <- data.frame(sample = order, label = unname(labels[order]),
    raw_calls = unname(raw_counts), plot_calls = unname(plot_counts))
  utils::write.table(counts, file.path(output, "sample_call_counts.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  publication_excel_out <- file.path(output, "publication_mutation_results.xlsx")
  publication_overrides <- config$publication_annotation_overrides
  if (is.character(publication_overrides) && length(publication_overrides) == 1L) {
    if (!file.exists(publication_overrides)) stop("Publication annotation file does not exist: ", publication_overrides)
    source_files <- c(source_files, publication_overrides)
    publication_overrides <- utils::read.delim(publication_overrides, check.names = FALSE,
      stringsAsFactors = FALSE, comment.char = "")
  }
  if (!is.null(publication_overrides)) utils::write.table(publication_overrides,
    file.path(work, "publication_annotations.tsv"), sep = "\t", row.names = FALSE,
    quote = TRUE, qmethod = "double", na = "")
  save_publication_results_to_excel(merged, publication_excel_out,
    sample_order = order, sample_labels = labels, reference_file = reference,
    gff_file = gff, plot_evidence = display_rows,
    annotation_overrides = publication_overrides)
  genome_export <- NULL
  create_genbank <- !identical(config$create_genbank, FALSE)
  if (create_genbank) {
    losses <- NULL
    if (!is.null(display_rows)) {
      loss_at <- display_rows$evidence_origin == "COVERAGE_RATIO" & display_rows$appended
      losses <- display_rows[loss_at, , drop = FALSE]
      if (nrow(losses)) losses$sample <- losses$sample_id
    }
    genome_export <- export_changed_sequences_to_genbank(reference, raw_calls,
      file.path(output, "changed_sequences"), samples = order, sample_labels = labels,
      replicon_losses = losses)
    for (artifact in c("feature_audit", "edits", "call_report")) {
      table <- genome_export[[artifact]]
      if (is.data.frame(table) && ncol(table)) utils::write.table(table,
        file.path(output, "changed_sequences", paste0(artifact, ".tsv")),
        sep = "\t", row.names = FALSE, quote = TRUE, na = "NA")
    }
    if (any(genome_export$manifest$status == "unresolved_calls")) {
      utils::write.table(genome_export$manifest,
        file.path(output, "changed_sequences", "sequence_manifest.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
      stop("Some finalized calls could not be applied exactly. See changed_sequences/call_report.tsv and sequence_manifest.tsv; DNA conversion was not started.")
    }
    unresolved <- .breseq_standalone_unresolved_coverage(mutation_data$Missing, order)
    if (nrow(unresolved)) {
        utils::write.table(unresolved, file.path(output, "changed_sequences", "unresolved_coverage_regions.tsv"),
          sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
        key <- paste(unresolved$sample, .breseq_seq(unresolved$seq_id), sep = "\r")
        manifest_key <- paste(genome_export$manifest$sample, .breseq_seq(genome_export$manifest$seq_id), sep = "\r")
        genome_export$manifest$unresolved_coverage <- manifest_key %in% key
    }
    if (!"unresolved_coverage" %in% names(genome_export$manifest))
      genome_export$manifest$unresolved_coverage <- FALSE
    genome_export$manifest$sequence_basis <- ifelse(genome_export$manifest$unresolved_coverage,
      "Finalized calls only; unresolved missing coverage was not applied.",
      "Finalized calls relative to the supplied reference.")
    genome_export$manifest$dna_file <- NA_character_
    genome_export$manifest$dna_status <- "not_requested"
    if (!identical(config$create_dna, FALSE)) {
      for (i in which(!is.na(genome_export$manifest$genbank_file) & nzchar(genome_export$manifest$genbank_file))) {
        gbk <- genome_export$manifest$genbank_file[i]
        dna <- sub("[.]gb[k]?$", ".dna", gbk, ignore.case = TRUE)
        converted <- .breseq_convert_genbank_to_dna(gbk, dna, config$snapgene_path)
        genome_export$manifest$dna_status[i] <- converted$status
        if (isTRUE(converted$success)) genome_export$manifest$dna_file[i] <- converted$output_file
        else message("SnapGene export ", converted$status, ": ", converted$message)
      }
    }
    utils::write.table(genome_export$manifest,
      file.path(output, "changed_sequences", "sequence_manifest.tsv"), sep = "\t",
      row.names = FALSE, quote = FALSE, na = "NA")
  }
  plots <- character()
  plot_dir <- file.path(output, "plots")
  dir.create(plot_dir, showWarnings = FALSE)
  ref_ids <- .breseq_seq(refs$info$seq_id)
  known <- unique(unlist(lapply(calls, function(d) if (nrow(d)) .breseq_seq(d[["seq id"]]) else character())))
  if (any(!known %in% ref_ids)) stop("Mutation replicons missing from the reference: ", paste(setdiff(known, ref_ids), collapse = ", "))
  for (rep in refs$replicons) {
    seqid <- .breseq_seq(rep$genome_info$seq_id[1L])
    per_calls <- lapply(calls, function(d) d[.breseq_seq(d[["seq id"]]) == seqid, , drop = FALSE])
    per_gff <- if (is.null(gff_data)) NULL else gff_data[.breseq_seq(gff_data$seqid) == seqid, , drop = FALSE]
    gc <- calculate_gc_content(rep$reference_file,
      window_size = min(10000L, rep$genome_info$length[1L]), gc_colors = config$gc_colors)
    if (is.null(gc)) stop("Could not calculate GC statistics for ", seqid)
    filename <- file.path(plot_dir, paste0("circos_", gsub("[^[:alnum:]_-]", "_", seqid), ".pdf"))
    plot_calls <- setNames(per_calls, unname(labels[names(per_calls)]))
    create_full_circos_plot(mutation_list = plot_calls, missing_data = NULL,
      samples = unname(labels[order]), genome_info = rep$genome_info, gff_data = per_gff,
      eggnog_data = eggnog_data, gc_data = gc, output_file = filename,
      hide_intergenic_labels = config$hide_intergenic_labels, gc_colors = config$gc_colors,
      observed_tracks = TRUE, reverse_tracks = FALSE)
    if (!file.exists(filename) || file.info(filename)$size == 0) stop("Plot export failed: ", filename)
    plots[seqid] <- filename
    saveRDS(gc, file.path(work, paste0("gc_", seqid, ".rds")))
  }
  source_files <- unique(normalizePath(source_files, mustWork = TRUE))
  manifest <- data.frame(path = source_files, md5 = unname(tools::md5sum(source_files)))
  utils::write.table(manifest, file.path(work, "source_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  result <- list(output_dir = output, excel_file = excel_out,
    publication_excel_file = publication_excel_out,
    excel_files = c(comparison = excel_out, publication = publication_excel_out), plots = plots,
    samples = order, sample_labels = labels, counts = counts,
    mutation_list = calls, merged_results = merged, display_rows = display_rows,
    reference_info = refs$info, genome_export = genome_export)
  saveRDS(result, file.path(work, "one_shot_results.rds"))
  message("Complete: ", output, "\nCreated ", length(plots), " replicon plots and both Excel workbooks (comparison and publication).")
  if (!is.null(genome_export)) message("Changed-sequence files: ", file.path(output, "changed_sequences"))
  invisible(result)
}
