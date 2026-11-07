# Reconstruct one sample/replicon from finalized breseq mutation calls.
# Coordinate rules: https://breseq.barricklab.org/latest/genomediff-file-format/
# HTML repeat notation retains the underlying INS/DEL position, rather than
# changing it to the start of the repeat tract (breseq src/breseq/output.cpp;
# MutationPredictor::normalize_and_annotate_tandem_repeat_mutations).
# No unassigned MC evidence or coverage-derived plotting intervals belong here.
.breseq_apply_sequence_edits <- function(sequence, calls, strict = TRUE) {
  if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence) ||
      !grepl("^[ACGTRYSWKMBDHVN]*$", toupper(sequence))) {
    stop("sequence must be one DNA string containing only IUPAC nucleotide symbols.")
  }
  if (!is.data.frame(calls)) stop("calls must be a data frame.")
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be a single TRUE/FALSE value.")
  }
  sequence <- toupper(sequence)
  reference_length <- nchar(sequence)
  empty_edits <- data.frame(position = numeric(), mutation = character(), type = character(),
    ref_start = numeric(), ref_end = numeric(), ref_length = numeric(), ref_sequence = character(),
    alt_length = numeric(), alt_sequence = character(), new_start = numeric(), new_end = numeric(),
    source_rows = character(), stringsAsFactors = FALSE)
  report <- data.frame(row = seq_len(nrow(calls)), position = rep(NA_character_, nrow(calls)),
    mutation = rep(NA_character_, nrow(calls)), status = rep("pending", nrow(calls)),
    message = rep("", nrow(calls)), edit_index = rep(NA_integer_, nrow(calls)),
    stringsAsFactors = FALSE)
  if (!nrow(calls)) return(list(sequence = sequence, edits = empty_edits, report = report,
                              complete = TRUE, reference_length = reference_length))
  if (!all(c("position", "mutation") %in% names(calls))) {
    stop("calls must contain position and mutation columns.")
  }
  seq_columns <- which(gsub("[[:space:]\u00a0_]", "", names(calls), perl = TRUE) == "seqid")
  if (length(seq_columns)) {
    ids <- unlist(lapply(calls[seq_columns], as.character), use.names = FALSE)
    ids <- unique(sub("\\.[0-9]+$", "", ids[!is.na(ids) & nzchar(ids)]))
    if (length(ids) > 1L) stop("calls contain multiple replicons; apply edits to one replicon at a time.")
  }
  report$position <- as.character(calls$position)
  report$mutation <- as.character(calls$mutation)
  # Nonbreaking/zero-width spaces come from breseq's HTML display.
  normalized <- gsub("[[:space:]\u00a0\u200b\ufeff]", "", report$mutation, perl = TRUE)
  normalized <- gsub("\u2192", "->", normalized, fixed = TRUE)
  normalized <- gsub("\u2212", "-", normalized, fixed = TRUE)
  normalized <- gsub("\u2206", "\u0394", normalized, fixed = TRUE)
  positions <- gsub(",", "", trimws(report$position), fixed = TRUE)
  matches <- function(pattern, text) regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
  fail <- function(status, message) {
    stop(structure(list(message = message, call = NULL, status = status),
                   class = c("breseq_edit_parse_error", "error", "condition")))
  }
  parsed <- vector("list", nrow(calls))
  for (i in seq_len(nrow(calls))) {
    parsed[i] <- list(tryCatch({
      if (is.na(positions[i]) || !grepl("^[0-9]+$", positions[i])) {
        fail("invalid_position", "Expected an integer reference coordinate, without ranges or junction symbols.")
      }
      position <- as.numeric(positions[i])
      if (!is.finite(position) || position < 1 || position > reference_length) {
        fail("invalid_position", "Position is outside the reference sequence.")
      }
      mutation <- normalized[i]
      if (is.na(mutation) || !nzchar(mutation)) fail("unsupported", "Mutation is missing.")
      mutation <- toupper(mutation)
      ref_start <- position
      ref_length <- 0
      alt <- ""
      expected <- NULL
      type <- NULL
      repeat_unit <- NULL
      old_copies <- NULL
      m <- matches("^([ACGT])->([ACGT])$", mutation)
      if (length(m)) {
        type <- "SNP"; ref_length <- 1; expected <- m[2]; alt <- m[3]
      } else if (grepl("^\\+[ACGT]+$", mutation)) {
        type <- "INS"; ref_start <- position + 1; alt <- substring(mutation, 2)
      } else if (length(m <- matches("^[-\u0394]([0-9][0-9,]*)(?:BP)?$", mutation))) {
        type <- "DEL"; ref_length <- as.numeric(gsub(",", "", m[2], fixed = TRUE))
      } else if (length(m <- matches("^-([ACGT]+)$", mutation))) {
        type <- "DEL"; expected <- m[2]; ref_length <- nchar(expected)
      } else if (length(m <- matches("^([0-9][0-9,]*)BP->([ACGT]+)$", mutation))) {
        type <- "REP"; ref_length <- as.numeric(gsub(",", "", m[2], fixed = TRUE)); alt <- m[3]
      } else if (length(m <- matches("^\\(([ACGT]+)\\)([0-9]+)->([0-9]+)$", mutation))) {
        repeat_unit <- m[2]; old_copies <- as.numeric(m[3]); new_copies <- as.numeric(m[4])
        difference <- new_copies - old_copies
        if (!is.finite(old_copies) || old_copies < 1 || !is.finite(new_copies) || difference == 0) {
          fail("invalid_repeat", "Repeat copy numbers must specify a change from a positive reference count.")
        }
        if (difference > 0) {
          type <- "INS"; ref_start <- position + 1
          alt <- strrep(repeat_unit, difference)
          repeat_end <- position
        } else {
          type <- "DEL"; ref_length <- -difference * nchar(repeat_unit)
          repeat_end <- position + ref_length - 1
          expected <- strrep(repeat_unit, -difference)
        }
        repeat_start <- repeat_end - old_copies * nchar(repeat_unit) + 1
        if (repeat_start < 1 || repeat_end > reference_length ||
            substr(sequence, repeat_start, repeat_end) != strrep(repeat_unit, old_copies)) {
          fail("reference_mismatch", "The reported repeat sequence/count does not match the reference at the breseq INS/DEL coordinate.")
        }
        # breseq counts complete units backwards at this phase, stopping before
        # reference position 1. Preserve that convention for origin-edge calls.
        preceding_start <- repeat_start - nchar(repeat_unit)
        if (preceding_start > 1 && substr(sequence, preceding_start, repeat_start - 1) == repeat_unit) {
          fail("reference_mismatch", "The reference contains more preceding repeat copies than reported.")
        }
      } else {
        fail("unsupported", "Unsupported or sequence-abbreviated mutation; supply an exact finalized call/GenomeDiff record.")
      }
      if (!is.finite(ref_length) || ref_length < 0 ||
          (type %in% c("DEL", "REP") && ref_length == 0)) {
        fail("invalid_length", "Deletion/replacement length must be a positive integer.")
      }
      ref_end <- ref_start + ref_length - 1
      if (ref_end > reference_length) fail("invalid_position", "Edit crosses the reference end; an origin-spanning edit needs explicit resolution.")
      reference <- if (ref_length) substr(sequence, ref_start, ref_end) else ""
      if (!is.null(expected) && !identical(reference, expected)) {
        fail("reference_mismatch", paste0("Reference allele mismatch: expected ", expected, ", found ", reference, "."))
      }
      if (ref_length && identical(reference, alt)) fail("no_change", "The call does not change the reference sequence.")
      report$status[i] <- "validated"
      data.frame(position = position, mutation = report$mutation[i], type = type,
        ref_start = ref_start, ref_end = ref_end, ref_length = ref_length, ref_sequence = reference,
        alt_length = nchar(alt), alt_sequence = alt, new_start = NA_real_, new_end = NA_real_,
        source_rows = as.character(i), stringsAsFactors = FALSE)
    }, breseq_edit_parse_error = function(e) {
      report$status[i] <<- e$status
      report$message[i] <<- conditionMessage(e)
      NULL
    }))
  }
  edits <- if (any(lengths(parsed))) do.call(rbind, parsed) else empty_edits
  # Equivalent exact edits are safe to collapse; every source row remains in the
  # report, including duplicate evidence rows and equivalent display spellings.
  if (nrow(edits)) {
    key <- paste(edits$ref_start, edits$ref_end, edits$alt_sequence, sep = "\r")
    unique_key <- unique(key)
    kept <- match(unique_key, key)
    sources <- lapply(unique_key, function(k) as.integer(edits$source_rows[key == k]))
    edits <- edits[kept, , drop = FALSE]
    edits$source_rows <- vapply(sources, paste, collapse = ",", character(1))
    edits <- edits[order(edits$ref_start, edits$ref_end), , drop = FALSE]
    rownames(edits) <- NULL
    sources <- lapply(strsplit(edits$source_rows, ",", fixed = TRUE), as.integer)
    for (j in seq_len(nrow(edits))) {
      report$edit_index[sources[[j]]] <- j
      if (length(sources[[j]]) > 1L) {
        duplicates <- sources[[j]][-1]
        report$status[duplicates] <- "duplicate"
        report$message[duplicates] <- paste0("Same exact edit as source row ", sources[[j]][1], ".")
      }
    }
    conflict <- logical(nrow(edits))
    if (nrow(edits) > 1L) for (a in seq_len(nrow(edits) - 1L)) for (b in seq.int(a + 1L, nrow(edits))) {
      ia <- edits$ref_length[a] == 0; ib <- edits$ref_length[b] == 0
      overlap <- if (ia && ib) {
        edits$ref_end[a] == edits$ref_end[b]
      } else if (!ia && !ib) {
        max(edits$ref_start[a], edits$ref_start[b]) <= min(edits$ref_end[a], edits$ref_end[b])
      } else {
        ins <- if (ia) a else b; other <- if (ia) b else a
        edits$ref_end[ins] >= edits$ref_start[other] && edits$ref_end[ins] < edits$ref_end[other]
      }
      if (overlap) conflict[c(a, b)] <- TRUE
    }
    for (j in which(conflict)) {
      report$status[sources[[j]]] <- "overlap"
      report$message[sources[[j]]] <- "Overlapping edits or differently ordered insertions require explicit resolution; no sequence was generated."
    }
  }
  invalid <- !report$status %in% c("validated", "duplicate")
  if (any(invalid)) {
    result <- list(sequence = NULL, edits = edits, report = report, complete = FALSE,
                   reference_length = reference_length)
    if (!strict) return(result)
    details <- paste0("row ", report$row[invalid], " (", report$position[invalid], ", ",
                      report$mutation[invalid], "): ", report$message[invalid])
    message <- paste("Cannot reconstruct sequence:", paste(utils::head(details, 8), collapse = "; "))
    if (sum(invalid) > 8L) message <- paste0(message, "; see $report for all ", sum(invalid), " unresolved rows.")
    stop(structure(c(list(message = message, call = NULL), result),
                   class = c("breseq_sequence_edit_error", "error", "condition")))
  }
  pieces <- character(2 * nrow(edits) + 1L)
  cursor <- 1
  shift <- 0
  for (j in seq_len(nrow(edits))) {
    pieces[2 * j - 1L] <- if (edits$ref_start[j] > cursor) substr(sequence, cursor, edits$ref_start[j] - 1) else ""
    pieces[2 * j] <- edits$alt_sequence[j]
    edits$new_start[j] <- edits$ref_start[j] + shift
    edits$new_end[j] <- edits$new_start[j] + edits$alt_length[j] - 1
    cursor <- edits$ref_end[j] + 1
    shift <- shift + edits$alt_length[j] - edits$ref_length[j]
  }
  pieces[length(pieces)] <- if (cursor <= reference_length) substr(sequence, cursor, reference_length) else ""
  mutated <- paste0(pieces, collapse = "")
  stopifnot(nchar(mutated) == reference_length + sum(edits$alt_length - edits$ref_length))
  report$status[report$status == "validated"] <- "applied"
  list(sequence = mutated, edits = edits, report = report, complete = TRUE,
       reference_length = reference_length)
}
