# Shared presentation helpers. Canonical data and sample IDs are never renamed.
.breseq_save_table_workbook <- function(wb, file) {
  # openxlsx initializes unused drawing relationships even on text-only sheets.
  # Some readers follow those links and fail when their target parts do not exist.
  # Retain relationships when actual drawings, VML or comments are present.
  for (i in seq_along(wb$worksheets_rels)) {
    rels <- wb$worksheets_rels[[i]]
    if (!length(wb$drawings[[i]]))
      rels <- rels[!grepl('/drawing"', rels, fixed = TRUE)]
    if (!length(wb$vml[[i]]) && !length(wb$comments[[i]]))
      rels <- rels[!grepl('/vmlDrawing"', rels, fixed = TRUE)]
    wb$worksheets_rels[[i]] <- rels
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

.breseq_text <- function(x) trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))

.breseq_stage_columns <- function(data, sample_order = NULL) {
  fields <- c("Forward aligned reads", "Reverse aligned reads", "Read Coverage Depth",
              "aligned reads1", "aligned reads2", "reads (cov).1", "reads (cov)",
              "description", "annotation", "evidence", "mutation", "product",
              "seq id", "seq_id", "gene", "end", "size", "\u2190reads", "reads\u2192",
              "score", "skew", "freq")
  fields <- fields[order(nchar(fields), decreasing = TRUE)]
  found <- list()
  for (i in seq_along(names(data))) {
    nm <- .breseq_text(names(data)[i])
    if (nm %in% fields) next
    for (field in fields) {
      prefix <- paste0(field, ".")
      if (!startsWith(nm, prefix)) next
      sample <- substring(nm, nchar(prefix) + 1L)
      sample <- sub("/output(?:/index)?$", "", sample, perl = TRUE)
      if (nzchar(sample)) found[[length(found) + 1L]] <-
        data.frame(index = i, field = if (field == "seq_id") "seq id" else field,
                   sample = sample, stringsAsFactors = FALSE)
      break
    }
  }
  result <- if (length(found)) do.call(rbind, found) else
    data.frame(index = integer(), field = character(), sample = character())
  detected <- unique(result$sample)
  if (!is.null(sample_order)) {
    if (!is.character(sample_order) || anyNA(sample_order) || any(!nzchar(sample_order)) ||
        anyDuplicated(sample_order) || !all(detected %in% sample_order))
      stop("sample_order must contain each detected sample exactly once (additional empty baselines are allowed).")
  }
  if (anyDuplicated(paste(result$field, result$sample, sep = "\r")))
    stop("Ambiguous duplicate field/sample columns in the wide table.")
  attr(result, "sample_order") <- if (is.null(sample_order)) detected else sample_order
  result
}

.breseq_table_kind <- function(name) {
  name <- tolower(name)
  if (grepl("missing|coverage|^source mc$", name)) return("missing")
  if (grepl("junction", name)) return("junction")
  if (grepl("mutation", name)) return("mutation")
  "other"
}

.breseq_seq <- function(x) sub("\\.[0-9]+$", "", .breseq_text(x))
.breseq_coord <- function(x) gsub("[,[:space:]]", "", .breseq_text(x))
.breseq_allele <- function(x) gsub("\u2206", "\u0394", .breseq_coord(x), fixed = TRUE)

# Keys are kept as an attribute, so legacy exported columns remain unchanged.
.breseq_source_keys <- function(data, kind, default_seq = "") {
  n <- nrow(data)
  if (!n) return(character())
  if ("comparison_key" %in% names(data)) return(as.character(data$comparison_key))
  seq_col <- intersect(c("seq id", "seq_id", "seqid"), names(data))
  seqid <- if (length(seq_col)) .breseq_seq(data[[seq_col[1]]]) else rep(default_seq, n)
  seqid[is.na(seqid)] <- default_seq
  coordinate <- if (kind == "missing") "start" else "position"
  endpoint <- paste(seqid, .breseq_coord(data[[coordinate]]), sep = "::")
  if (kind != "junction") return(endpoint)
  if (n %% 2L) stop("Junction tables must contain two endpoint rows per event.")
  keys <- character(n)
  for (i in seq.int(1L, n, by = 2L)) {
    pair <- endpoint[i:(i + 1L)]
    ordering <- order(pair)
    roles <- order(ordering)
    keys[i:(i + 1L)] <- paste(paste(pair[ordering], collapse = "||"),
                             paste0("end", roles), sep = "||")
  }
  keys
}

.breseq_wide_seq <- function(data, columns) {
  result <- rep(NA_character_, nrow(data))
  globals <- intersect(c("seq_id", "seq id", "seqid"), names(data))
  indices <- unique(c(match(globals, names(data)), columns$index[columns$field == "seq id"]))
  for (i in indices) {
    x <- .breseq_seq(data[[i]])
    ok <- !is.na(x) & nzchar(x)
    conflict <- ok & !is.na(result) & result != x
    if (any(conflict)) stop("A merged row contains conflicting replicon IDs; rebuild comparisons with composite keys.")
    result[ok & is.na(result)] <- x[ok & is.na(result)]
  }
  result[is.na(result)] <- ""
  result
}

.breseq_field <- function(data, columns, field, sample) {
  at <- columns$index[columns$field == field & columns$sample == sample]
  if (!length(at)) return(rep(NA_character_, nrow(data)))
  .breseq_text(data[[at[1]]])
}
.breseq_present <- function(x) !is.na(x) & nzchar(x)
.breseq_bounds <- function(x) {
  if (!.breseq_present(x)) return(c(NA_real_, NA_real_))
  values <- as.numeric(regmatches(gsub(",", "", x, fixed = TRUE),
                                 gregexpr("[0-9]+", gsub(",", "", x, fixed = TRUE)))[[1]])
  if (!length(values) || anyNA(values)) c(NA_real_, NA_real_) else range(values)
}

# Return row-by-stage states without modifying data or using changing read metrics
# as event identity. MC states describe evidence, never newly acquired deletions.
.breseq_stage_status <- function(data, table_name, sample_order = NULL, mutation_data = NULL) {
  columns <- .breseq_stage_columns(data, sample_order)
  samples <- attr(columns, "sample_order")
  states <- matrix("not_reported", nrow(data), length(samples), dimnames = list(NULL, samples))
  kind <- .breseq_table_kind(table_name)
  if (!nrow(data) || !length(samples) || kind == "other") return(list(samples = samples, states = states))
  seqid <- .breseq_wide_seq(data, columns)
  if (kind == "mutation") {
    if (!"position" %in% names(data)) stop("Mutation tables require position for stage highlighting.")
    locus <- paste(seqid, .breseq_coord(data$position), sep = "::")
    if (anyDuplicated(locus)) stop("Duplicate mutation loci cannot be highlighted safely; rebuild the merged table.")
    history <- vector("list", nrow(data))
    previous <- rep(NA_character_, nrow(data))
    for (sample in samples) {
      allele <- .breseq_allele(.breseq_field(data, columns, "mutation", sample))
      for (i in which(.breseq_present(allele))) {
        states[i, sample] <- if (.breseq_present(previous[i])) {
          if (allele[i] == previous[i]) "retained" else "changed"
        } else if (allele[i] %in% history[[i]]) "reappeared" else "new"
        history[[i]] <- unique(c(history[[i]], allele[i]))
      }
      previous <- allele
    }
  } else if (kind == "junction") {
    keys <- if ("comparison_key" %in% names(data)) as.character(data$comparison_key) else attr(data, "breseq_keys")
    if (is.null(keys) || length(keys) != nrow(data) || anyNA(keys) || any(!grepl("\\|\\|", keys)))
      stop("Junction highlighting requires paired comparison keys from compare_multiple_tables().")
    event <- sub("\\|\\|end[12]$", "", keys)
    groups <- split(seq_len(nrow(data)), event)
    if (any(lengths(groups) != 2L)) stop("Every junction event must retain both endpoint rows.")
    previous <- character(); seen <- character()
    for (sample in samples) {
      present <- .breseq_present(.breseq_field(data, columns, "seq id", sample))
      if (!any(columns$field == "seq id" & columns$sample == sample))
        present <- .breseq_present(.breseq_field(data, columns, "reads (cov)", sample))
      current <- character()
      for (e in names(groups)) {
        at <- groups[[e]]
        if (any(present[at]) && !all(present[at])) stop("Only one endpoint is present for a junction sample/event.")
        if (!all(present[at])) next
        states[at, sample] <- if (e %in% previous) "retained" else if (e %in% seen) "reappeared" else "new"
        current <- c(current, e)
      }
      seen <- union(seen, current); previous <- current
    }
  } else if (kind == "missing") {
    if (!"start" %in% names(data)) stop("Missing-coverage tables require start for highlighting.")
    previous <- list(); history <- list()
    mcols <- if (!is.null(mutation_data)) .breseq_stage_columns(mutation_data, sample_order) else NULL
    mseq <- if (!is.null(mcols)) .breseq_wide_seq(mutation_data, mcols) else NULL
    for (j in seq_along(samples)) {
      sample <- samples[j]
      ends <- .breseq_field(data, columns, "end", sample)
      genes <- .breseq_field(data, columns, "gene", sample)
      current <- list()
      for (i in which(.breseq_present(ends))) {
        a <- .breseq_bounds(as.character(data$start[i])); b <- .breseq_bounds(ends[i])
        record <- list(seq = seqid[i], start = a[2], end = b[1],
                       raw = paste(.breseq_coord(data$start[i]), .breseq_coord(ends[i])), gene = genes[i])
        if (anyNA(c(record$start, record$end)) || record$start > record$end)
          stop("Cannot classify an invalid missing-coverage interval.")
        overlap <- function(other) {
          common_locus <- TRUE
          if (.breseq_present(record$gene) && .breseq_present(other$gene)) {
            loci <- function(x) regmatches(x, gregexpr("[A-Za-z0-9]+_RS[0-9]+", x))[[1]]
            left <- loci(record$gene); right <- loci(other$gene)
            common_locus <- if (length(left) && length(right)) length(intersect(left, right)) > 0L else
              identical(record$gene, other$gene)
          }
          identical(record$seq, other$seq) && common_locus &&
            max(record$start, other$start) <= min(record$end, other$end)
        }
        old <- Filter(overlap, previous)
        state <- if (length(old)) {
          if (any(vapply(old, function(o) identical(record$raw, o$raw), logical(1)))) "retained" else "updated"
        } else if (length(Filter(overlap, history))) "reappeared" else "new"
        if (state == "new" && j > 1L && !is.null(mcols)) {
          alleles <- .breseq_field(mutation_data, mcols, "mutation", samples[j - 1L])
          is_del <- !is.na(alleles) & grepl("^[\u0394\u2206][0-9,]+[[:space:]]*bp$", alleles)
          size <- suppressWarnings(as.numeric(gsub("[^0-9]", "", alleles)))
          pos <- suppressWarnings(as.numeric(gsub(",", "", mutation_data$position, fixed = TRUE)))
          covered <- is_del & mseq == record$seq & pos <= record$start & pos + size - 1 >= record$end
          if (any(covered, na.rm = TRUE)) state <- "prior_del_region"
        }
        states[i, sample] <- state
        current[[length(current) + 1L]] <- record
      }
      history <- c(history, current); previous <- current
    }
  }
  list(samples = samples, states = states)
}

.breseq_validate_labels <- function(labels, samples) {
  if (is.null(labels)) return(stats::setNames(samples, samples))
  if (!is.character(labels) || is.null(names(labels)) || anyNA(labels) || any(!nzchar(labels)) ||
      anyNA(names(labels)) || any(!nzchar(names(labels))) || anyDuplicated(names(labels)) ||
      !all(names(labels) %in% samples)) stop("sample_labels must be a named character map of known sample IDs.")
  result <- stats::setNames(samples, samples); result[names(labels)] <- labels
  if (anyDuplicated(unname(result))) stop("sample_labels must produce distinct display labels.")
  result
}

.breseq_write_presentation <- function(wb, sheet, data, table_name, highlight_new = FALSE,
                                     clean_headers = FALSE, sample_labels = NULL,
                                     sample_order = NULL, mutation_data = NULL) {
  sheet <- unname(as.character(sheet))
  if (!is.logical(highlight_new) || length(highlight_new) != 1L || is.na(highlight_new) ||
      !is.logical(clean_headers) || length(clean_headers) != 1L || is.na(clean_headers))
    stop("highlight_new and clean_headers must be single TRUE/FALSE values.")
  if (!highlight_new && !clean_headers) {
    openxlsx::writeData(wb, sheet, data)
    return(invisible(NULL))
  }
  openxlsx::modifyBaseFont(wb, fontSize = 10, fontColour = "#333333", fontName = "Arial")
  body_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10,
    fontColour = "#333333", wrapText = TRUE, valign = "center")
  header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10,
    fontColour = "#333333", textDecoration = "bold", fgFill = "#FFFFFF",
    wrapText = TRUE, valign = "center")
  columns <- .breseq_stage_columns(data, sample_order)
  samples <- attr(columns, "sample_order")
  labels <- .breseq_validate_labels(sample_labels, samples)
  displayed <- data
  if (clean_headers) {
    base <- c(seq_id = "Replicon", position = "Position (bp)", start = "Start (bp)", comparison_key = "Comparison key")
    at <- match(names(displayed), names(base)); use <- !is.na(at)
    names(displayed)[use] <- unname(base[at[use]])
    for (k in seq_len(nrow(columns))) {
      field <- columns$field[k]
      field <- paste0(toupper(substring(field, 1, 1)), substring(field, 2))
      names(displayed)[columns$index[k]] <- paste(unname(labels[columns$sample[k]]), field, sep = " | ")
    }
    if (anyDuplicated(names(displayed))) stop("Clean headers would be ambiguous.")
    metadata <- data.frame(sheet = sheet, column = seq_along(data), canonical = names(data),
                           display = names(displayed), sample_order = paste(samples, collapse = "\t"),
                           stringsAsFactors = FALSE)
    meta_sheet <- "_breseq_columns"
    if (!meta_sheet %in% names(wb)) {
      openxlsx::addWorksheet(wb, meta_sheet, visible = FALSE)
      openxlsx::writeData(wb, meta_sheet, metadata)
      metadata_start <- 2L
    } else {
      old <- openxlsx::readWorkbook(wb, sheet = meta_sheet, check.names = FALSE)
      metadata_start <- nrow(old) + 2L
      openxlsx::writeData(wb, meta_sheet, metadata, startRow = metadata_start, colNames = FALSE)
    }
    openxlsx::addStyle(wb, meta_sheet, header_style, rows = 1L,
      cols = seq_len(ncol(metadata)), gridExpand = TRUE)
    if (nrow(metadata)) openxlsx::addStyle(wb, meta_sheet, body_style,
      rows = seq.int(metadata_start, length.out = nrow(metadata)),
      cols = seq_len(ncol(metadata)), gridExpand = TRUE)
  }
  openxlsx::writeData(wb, sheet, displayed)
  if (ncol(displayed)) {
    widths <- rep(20, ncol(displayed))
    for (k in seq_len(nrow(columns))) {
      field <- columns$field[k]
      widths[columns$index[k]] <- if (field %in% c("description", "product")) 60 else
        if (field == "annotation") 42 else if (field == "gene") 36 else
        if (field == "mutation") 24 else 22
    }
    widths[names(data) %in% c("position", "start")] <- 17
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(displayed)), widths = widths)
    openxlsx::setRowHeights(wb, sheet, rows = 1L, heights = 44)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 2L,
      firstActiveCol = min(3L, ncol(displayed) + 1L))
    openxlsx::addStyle(wb, sheet, header_style, rows = 1L,
      cols = seq_len(ncol(displayed)), gridExpand = TRUE)
    if (nrow(displayed)) {
      openxlsx::addStyle(wb, sheet, body_style, rows = seq_len(nrow(displayed)) + 1L,
        cols = seq_len(ncol(displayed)), gridExpand = TRUE)
      heights <- vapply(seq_len(nrow(displayed)), function(i) {
        lines <- vapply(seq_len(ncol(displayed)), function(j) {
          text <- as.character(displayed[[j]][i]); if (is.na(text)) return(1)
          sum(pmax(1, ceiling(nchar(strsplit(text, "\n", fixed = TRUE)[[1]]) / (widths[j] - 3))))
        }, 1)
        min(409.5, max(20, max(lines) * 14 + 6))
      }, 1)
      openxlsx::setRowHeights(wb, sheet, rows = seq_len(nrow(displayed)) + 1L, heights = heights)
    }
  }
  if (highlight_new && nrow(data)) {
    status <- .breseq_stage_status(data, table_name, sample_order, mutation_data)
    # Color means observed in this sample, including retained calls. Bold marks
    # a first report, changed allele or reappearance; it is not a de novo claim.
    palette <- c("#EEF3F8", "#F0F5EF", "#F7F2EA", "#F3F0F6")
    is_wt <- toupper(status$samples) == "WT"
    fills <- stats::setNames(rep("#F3F3F3", length(status$samples)), status$samples)
    fills[!is_wt] <- rep(palette, length.out = sum(!is_wt))
    bold_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10,
      fontColour = "#333333", textDecoration = "bold")
    for (sample in status$samples) {
      rows <- which(status$states[, sample] != "not_reported")
      cols <- columns$index[columns$sample == sample]
      if (length(rows) && length(cols)) openxlsx::addStyle(wb, sheet,
        openxlsx::createStyle(fgFill = fills[[sample]]), rows = rows + 1L, cols = cols,
        gridExpand = TRUE, stack = TRUE)
      bold_rows <- which(status$states[, sample] %in% c("new", "changed", "reappeared"))
      if (length(bold_rows) && length(cols)) openxlsx::addStyle(wb, sheet,
        bold_style, rows = bold_rows + 1L, cols = cols, gridExpand = TRUE, stack = TRUE)
    }
  }
  invisible(NULL)
}

#' Restore Canonical Headers from a Presentation Export
#'
#' Internal loader used for clean-header workbook round trips. Files without
#' breseqConverter metadata are returned unchanged.
#' @keywords internal
restore_canonical_export_headers <- function(file, sheet, data) {
  metadata_sheet <- "_breseq_columns"
  if (!metadata_sheet %in% openxlsx::getSheetNames(file)) return(data)
  meta <- openxlsx::read.xlsx(file, sheet = metadata_sheet, check.names = FALSE)
  if (!all(c("sheet", "column", "canonical", "display", "sample_order") %in% names(meta)))
    stop("Invalid breseqConverter canonical-column metadata.")
  if (is.numeric(sheet)) sheet <- openxlsx::getSheetNames(file)[sheet]
  meta <- meta[meta$sheet == sheet, , drop = FALSE]
  if (!nrow(meta)) return(data)
  meta <- meta[order(as.integer(meta$column)), , drop = FALSE]
  if (nrow(meta) != ncol(data) || !identical(as.integer(meta$column), seq_len(ncol(data))))
    stop("Canonical-column metadata does not match the loaded worksheet; preserve empty columns when reading.")
  displayed_names <- as.character(meta$display)
  # read.xlsx(check.names=FALSE) still uses sep.names='.' by default.
  if (!identical(names(data), displayed_names) &&
      !identical(names(data), gsub(" ", ".", displayed_names, fixed = TRUE)))
    stop("Worksheet headers differ from the presentation metadata; refusing to assign incorrect canonical names.")
  names(data) <- as.character(meta$canonical)
  order <- unique(as.character(meta$sample_order))
  order <- order[!is.na(order) & nzchar(order)]
  if (length(order) == 1L) attr(data, "breseq_sample_order") <- strsplit(order, "\t", fixed = TRUE)[[1]]
  data
}
