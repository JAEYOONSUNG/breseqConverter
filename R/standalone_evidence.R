# Add explicitly audited display evidence without changing canonical call rows.
# The sidecar is an input: MC cores are not inferred again from boundary depths.
.breseq_standalone_display_deletions <- function(mutation_list, evidence_file,
                                                sample_aliases, reference_info,
                                                merged_missing = NULL) {
  if (!is.character(evidence_file) || length(evidence_file) != 1L ||
      is.na(evidence_file) || !file.exists(evidence_file))
    stop("evidence_file must name an existing audited TSV file.")
  evidence <- utils::read.delim(evidence_file, check.names = FALSE,
                              stringsAsFactors = FALSE)
  required <- c("sample", "seq_id", "start", "end", "size", "evidence_origin",
                "source_start", "source_end", "depth_ratio", "criterion",
                "raw_final_DEL_call")
  if (!all(required %in% names(evidence)))
    stop("Audited evidence TSV is missing required columns: ",
         paste(setdiff(required, names(evidence)), collapse = ", "))
  evidence$sample_id <- character(nrow(evidence))
  evidence$appended <- rep(FALSE, nrow(evidence))
  evidence$skip_reason <- rep("", nrow(evidence))
  if (!nrow(evidence)) return(list(mutation_list = mutation_list, display_rows = evidence))

  samples <- names(mutation_list)
  if (!is.list(mutation_list) || is.null(samples) || anyNA(samples) ||
      any(!nzchar(samples)) || anyDuplicated(samples))
    stop("mutation_list must have unique, non-empty sample IDs in analysis order.")
  if (!is.character(sample_aliases) || is.null(names(sample_aliases)) ||
      anyNA(sample_aliases) || anyNA(names(sample_aliases)) ||
      any(!nzchar(sample_aliases)) || any(!nzchar(names(sample_aliases))) ||
      anyDuplicated(names(sample_aliases)))
    stop("sample_aliases must map unique source aliases to actual sample IDs.")
  source_samples <- .breseq_text(evidence$sample)
  mapped <- source_samples
  alias_at <- match(source_samples, names(sample_aliases))
  mapped[!is.na(alias_at)] <- unname(sample_aliases[alias_at[!is.na(alias_at)]])
  if (anyNA(mapped) || any(!mapped %in% samples))
    stop("Audited evidence contains an unknown sample or sample alias.")
  evidence$sample_id <- mapped

  if (!all(c("seq_id", "length") %in% names(reference_info)))
    stop("reference_info must contain seq_id and length.")
  reference_ids <- .breseq_seq(reference_info$seq_id)
  reference_lengths <- suppressWarnings(as.numeric(reference_info$length))
  if (anyNA(reference_ids) || any(!nzchar(reference_ids)) || anyDuplicated(reference_ids) ||
      any(!is.finite(reference_lengths) | reference_lengths < 1 |
          reference_lengths != floor(reference_lengths)))
    stop("reference_info must have unique replicon IDs and positive integer lengths.")
  seqids <- .breseq_seq(evidence$seq_id)
  reference_at <- match(seqids, reference_ids)
  if (anyNA(reference_at)) stop("Audited evidence contains an unknown reference replicon.")
  coordinates <- lapply(evidence[c("start", "end", "size")], function(x)
    suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE))))
  start <- coordinates$start; end <- coordinates$end; size <- coordinates$size
  if (any(!is.finite(c(start, end, size))) ||
      any(start != floor(start) | end != floor(end) | size != floor(size)) ||
      any(start < 1 | end < start | size < 1 | end > reference_lengths[reference_at]))
    stop("Audited evidence coordinates must be positive, within-reference integer intervals.")
  if (any(size != end - start + 1)) stop("Audited evidence size must equal end - start + 1.")
  origin <- .breseq_text(evidence$evidence_origin)
  if (anyNA(origin) || any(!origin %in% c("MC", "COVERAGE_RATIO")))
    stop("Audited evidence_origin must be MC or COVERAGE_RATIO.")
  final_call <- toupper(.breseq_text(evidence$raw_final_DEL_call))
  if (anyNA(final_call) || any(final_call != "FALSE"))
    stop("Audited display evidence must have raw_final_DEL_call = FALSE.")
  ratio <- suppressWarnings(as.numeric(evidence$depth_ratio))
  coverage <- origin == "COVERAGE_RATIO"
  if (any(!is.finite(ratio[coverage]) | ratio[coverage] < 0 | ratio[coverage] >= 0.1) ||
      any(start[coverage] != 1 | end[coverage] != reference_lengths[reference_at[coverage]]))
    stop("COVERAGE_RATIO evidence requires depth_ratio <0.1 and the whole replicon interval.")
  if (anyDuplicated(paste(mapped, seqids, start, end, sep = "\r")))
    stop("Audited evidence contains duplicate sample/replicon intervals.")

  missing_columns <- if (!is.null(merged_missing)) .breseq_stage_columns(merged_missing) else NULL
  missing_seqids <- if (!is.null(missing_columns)) .breseq_wide_seq(merged_missing, missing_columns) else NULL
  source_keys <- function(data) {
    if (!nrow(data)) return(character())
    seq_columns <- intersect(c("seq id", "seq_id", "seqid"), names(data))
    if (!length(seq_columns) || !all(c("position", "mutation") %in% names(data)))
      stop("Mutation rows require a replicon ID, position, and mutation.")
    paste(.breseq_seq(data[[seq_columns[1]]]), .breseq_coord(data$position),
          .breseq_allele(data$mutation), sep = "\r")
  }
  seen <- character()
  for (sample in samples) {
    data <- mutation_list[[sample]]
    if (!is.data.frame(data)) stop("Every mutation_list entry must be a data frame.")
    automatic_rownames <- .row_names_info(data, 1L) < 0L
    # Prior source calls and prior display rows determine first-appearance labels.
    prior_keys <- seen
    source_keys(data)
    for (i in which(mapped == sample)) {
      seq_columns <- intersect(c("seq id", "seq_id", "seqid"), names(data))
      needed <- c("position", "mutation", "annotation", "evidence", "gene",
                  "description", "new", "ALE_version", "version_color", "label")
      if (!length(seq_columns) || !all(needed %in% names(data)))
        stop("Display rows require the standard processed mutation columns.")
      if (any(vapply(data, is.factor, logical(1))))
        stop("Processed mutation columns must not be factors when appending display evidence.")
      same_replicon <- .breseq_seq(data[[seq_columns[1]]]) == seqids[i]
      deleted_size <- suppressWarnings(as.numeric(sub("^[\u0394\u2206]([0-9,]+).*", "\\1",
        .breseq_allele(data$mutation))))
      deleted_start <- suppressWarnings(as.numeric(.breseq_coord(data$position)))
      final_del <- same_replicon & grepl("^[\u0394\u2206][0-9,]+bp$", .breseq_allele(data$mutation)) &
        is.finite(deleted_start) & is.finite(deleted_size)
      deleted_end <- deleted_start + deleted_size - 1
      if (any(final_del & deleted_start <= start[i] & deleted_end >= end[i], na.rm = TRUE)) {
        evidence$skip_reason[i] <- "already covered by same-sample DEL"
        next
      }
      if (any(final_del & deleted_start <= end[i] & deleted_end >= start[i], na.rm = TRUE))
        stop("Audited display interval partially overlaps an existing same-sample DEL.")

      gene <- if (coverage[i]) paste0("p", seqids[i]) else ""
      if (!coverage[i] && !is.null(missing_columns)) {
        missing_end <- .breseq_field(merged_missing, missing_columns, "end", sample)
        at <- which(missing_seqids == seqids[i] &
          .breseq_coord(merged_missing$start) == .breseq_coord(evidence$source_start[i]) &
          .breseq_coord(missing_end) == .breseq_coord(evidence$source_end[i]))
        if (length(at) > 1L) stop("Ambiguous source MC row for an audited display interval.")
        if (length(at)) {
          source_gene <- .breseq_field(merged_missing, missing_columns, "gene", sample)[at]
          if (.breseq_present(source_gene)) gene <- source_gene
        }
      }
      row <- data[NA_integer_, , drop = FALSE]
      for (name in names(row)) {
        row[[name]] <- if (is.logical(data[[name]])) FALSE else
          if (is.integer(data[[name]])) NA_integer_ else
            if (is.numeric(data[[name]])) NA_real_ else ""
      }
      for (name in seq_columns) row[[name]] <- seqids[i]
      row$position <- if (is.integer(data$position)) as.integer(start[i]) else
        if (is.character(data$position)) as.character(start[i]) else start[i]
      row$mutation <- paste0("\u0394", format(size[i], big.mark = ",", scientific = FALSE, trim = TRUE), " bp")
      row$annotation <- "coverage-derived DEL display"
      row$evidence <- origin[i]
      row$gene <- gene
      row$description <- if (coverage[i]) "Whole-replicon depth / chromosome depth <0.1" else
        "Conservative internal span of source missing-coverage evidence"
      row$new <- !source_keys(row) %in% c(prior_keys, source_keys(data))
      row$ALE_version <- sample
      colors <- data$version_color[.breseq_present(data$version_color)]
      if (!length(colors)) stop("A display row requires the sample's existing version_color.")
      row$version_color <- colors[1]
      row$label <- trimws(paste(row$mutation, gene))
      if ("start" %in% names(row)) row$start <- if (is.integer(data$start)) as.integer(start[i]) else start[i]
      if ("end" %in% names(row)) row$end <- if (is.integer(data$end)) as.integer(end[i]) else end[i]
      if ("type" %in% names(row)) row$type <- "DEL"
      if ("mutation_type" %in% names(row)) row$mutation_type <- if (size[i] >= 50) "large_deletion" else "deletion"
      source_row_names <- attr(data, "row.names")
      rownames(row) <- make.unique(c(rownames(data), paste0("display_", i)))[nrow(data) + 1L]
      data[nrow(data) + 1L, ] <- row
      if (automatic_rownames) rownames(data) <- NULL else if (is.integer(source_row_names))
        attr(data, "row.names") <- c(source_row_names, max(c(0L, source_row_names)) + 1L)
      evidence$appended[i] <- TRUE
    }
    mutation_list[[sample]] <- data
    seen <- union(seen, source_keys(data))
  }
  list(mutation_list = mutation_list, display_rows = evidence)
}
