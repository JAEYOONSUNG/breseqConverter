# GenBank reconstruction uses edited bases, independently of display-only marks.
.breseq_gb_split_location <- function(text) {
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  depth <- 0L; starts <- 1L; result <- character()
  for (i in seq_along(chars)) {
    if (chars[i] == "(") depth <- depth + 1L
    if (chars[i] == ")") depth <- depth - 1L
    if (chars[i] == "," && depth == 0L) {
      result <- c(result, substr(text, starts, i - 1L)); starts <- i + 1L
    }
  }
  c(result, substr(text, starts, nchar(text)))
}

.breseq_gb_location <- function(text) {
  text <- gsub("[[:space:]]", "", text)
  if (grepl("^(complement|join|order)\\(", text)) {
    operator <- sub("\\(.*", "", text)
    if (!endsWith(text, ")")) stop("Unbalanced feature location")
    inner <- substr(text, nchar(operator) + 2L, nchar(text) - 1L)
    parts <- .breseq_gb_split_location(inner)
    if (operator == "complement" && length(parts) != 1L) stop("Invalid complement")
    return(list(type = operator, children = lapply(parts, .breseq_gb_location)))
  }
  if (grepl("^[<>]?[0-9]+(\\.\\.[<>]?[0-9]+)?$", text)) {
    parts <- strsplit(text, "..", fixed = TRUE)[[1]]
    if (length(parts) == 1L) parts <- rep(parts, 2L)
    return(list(type = "range", start = as.numeric(gsub("[<>]", "", parts[1])),
      end = as.numeric(gsub("[<>]", "", parts[2])),
      left = sub("[0-9].*", "", parts[1]), right = sub("[0-9].*", "", parts[2]),
      point = !grepl("..", text, fixed = TRUE)))
  }
  if (grepl("^[0-9]+\\^[0-9]+$", text)) {
    parts <- as.numeric(strsplit(text, "^", fixed = TRUE)[[1]])
    return(list(type = "between", start = parts[1], end = parts[2]))
  }
  stop("Unsupported feature location: ", text)
}

.breseq_gb_location_text <- function(node) {
  if (node$type == "range") {
    left <- paste0(node$left, node$start)
    if (isTRUE(node$point) && node$start == node$end) return(left)
    return(paste0(left, "..", node$right, node$end))
  }
  if (node$type == "between") return(paste0(node$start, "^", node$end))
  paste0(node$type, "(", paste(vapply(node$children, .breseq_gb_location_text, ""), collapse = ","), ")")
}

.breseq_gb_map_base <- function(position, edits) {
  position + sum((edits$alt_length - edits$ref_length)[edits$ref_end < position])
}

.breseq_gb_lift_location <- function(node, edits, reference_length, new_length) {
  if (!node$type %in% c("range", "between")) {
    lifted <- lapply(node$children, .breseq_gb_lift_location, edits = edits,
      reference_length = reference_length, new_length = new_length)
    children <- lapply(lifted, `[[`, "node")
    children <- Filter(Negate(is.null), children)
    changed <- any(vapply(lifted, `[[`, logical(1), "changed"))
    partial <- length(children) != length(node$children) || any(vapply(lifted, `[[`, logical(1), "partial"))
    if (!length(children)) return(list(node = NULL, changed = TRUE, partial = TRUE))
    node$children <- children
    if (node$type %in% c("join", "order") && length(children) == 1L) node <- children[[1]]
    return(list(node = node, changed = changed, partial = partial))
  }
  a <- node$start; b <- node$end
  if (node$type == "between") {
    if (!(b == a + 1 || (a == reference_length && b == 1))) stop("Invalid between-base location")
    inside <- edits$ref_length > 0 & edits$ref_start <= a & edits$ref_end > a
    if (any(inside)) return(list(node = NULL, changed = TRUE, partial = TRUE))
    boundary <- a + sum((edits$alt_length - edits$ref_length)[edits$ref_end <= a])
    if (a == reference_length && b == 1) {
      node$start <- new_length; node$end <- 1
    } else {
      if (boundary < 1 || boundary >= new_length) return(list(node = NULL, changed = TRUE, partial = TRUE))
      node$start <- boundary; node$end <- boundary + 1
    }
    return(list(node = node, changed = any(edits$ref_length == 0 & edits$ref_end == a), partial = FALSE))
  }
  if (a < 1 || b < a || b > reference_length) stop("Out-of-reference feature location")
  relevant <- which((edits$ref_length > 0 & edits$ref_start <= b & edits$ref_end >= a) |
    (edits$ref_length == 0 & edits$ref_end >= a & edits$ref_end < b))
  if (!length(relevant)) {
    node$start <- .breseq_gb_map_base(a, edits); node$end <- .breseq_gb_map_base(b, edits)
    return(list(node = node, changed = FALSE, partial = FALSE))
  }
  segments <- list(); cursor <- a; partial_left <- FALSE; partial_right <- FALSE
  add <- function(left, right) { if (left <= right) segments[[length(segments) + 1L]] <<- c(left, right) }
  for (i in relevant) {
    e <- edits[i, ]
    if (e$ref_length == 0) {
      add(e$new_start, e$new_end)
      next
    }
    before <- min(b, e$ref_start - 1)
    if (cursor <= before) add(.breseq_gb_map_base(cursor, edits), .breseq_gb_map_base(before, edits))
    if (e$alt_length > 0) {
      if (e$ref_length == e$alt_length) {
        add(e$new_start + max(a, e$ref_start) - e$ref_start,
            e$new_start + min(b, e$ref_end) - e$ref_start)
      } else if (a <= e$ref_start && b >= e$ref_end) add(e$new_start, e$new_end)
    }
    if (e$ref_length != e$alt_length && e$ref_start <= a && e$ref_end >= a) partial_left <- TRUE
    if (e$ref_length != e$alt_length && e$ref_start <= b && e$ref_end >= b) partial_right <- TRUE
    cursor <- max(cursor, e$ref_end + 1)
  }
  if (cursor <= b) add(.breseq_gb_map_base(cursor, edits), .breseq_gb_map_base(b, edits))
  if (!length(segments)) return(list(node = NULL, changed = TRUE, partial = TRUE))
  node$start <- min(vapply(segments, `[`, 0, 1)); node$end <- max(vapply(segments, `[`, 0, 2))
  if (partial_left) node$left <- "<"
  if (partial_right) node$right <- ">"
  if (node$start != node$end) node$point <- FALSE
  list(node = node, changed = TRUE, partial = partial_left || partial_right)
}

.breseq_gb_reverse_complement <- function(sequence) {
  paste(rev(strsplit(chartr("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN", sequence), "", fixed = TRUE)[[1]]), collapse = "")
}

.breseq_gb_feature_sequence <- function(node, sequence) {
  if (node$type == "range") return(substr(sequence, node$start, node$end))
  if (node$type == "between") return("")
  if (node$type == "order") stop("Ordered segments cannot be translated as a joined CDS")
  joined <- paste(vapply(node$children, .breseq_gb_feature_sequence, "", sequence = sequence), collapse = "")
  if (node$type == "complement") .breseq_gb_reverse_complement(joined) else joined
}

.breseq_gb_translation <- function(node, sequence, code = 1L, codon_start = 1L) {
  if (!code %in% c(1L, 11L)) stop("unsupported translation table")
  if (grepl("[<>]", .breseq_gb_location_text(node)) || codon_start != 1L) stop("partial CDS")
  dna <- .breseq_gb_feature_sequence(node, sequence)
  if (!nchar(dna) || nchar(dna) %% 3L != 0L) stop("CDS length is not a multiple of three")
  if (grepl("[^ACGT]", dna)) stop("ambiguous coding bases")
  starts <- seq.int(1L, nchar(dna), by = 3L)
  codons <- substring(dna, starts, starts + 2L)
  bases <- c("T", "C", "A", "G")
  # Explicit NCBI standard-code order: first base slowest, third base fastest.
  keys <- unlist(lapply(bases, function(a) unlist(lapply(bases, function(b) paste0(a, b, bases)))))
  amino <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "", fixed = TRUE)[[1]]
  protein <- unname(stats::setNames(amino, keys)[codons])
  allowed_starts <- if (code == 11L) c("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG") else c("TTG", "CTG", "ATG")
  if (!codons[1] %in% allowed_starts) stop("start codon is not valid for this translation table")
  if (tail(protein, 1) != "*" || any(head(protein, -1L) == "*")) stop("missing terminal stop or internal stop codon")
  protein[1] <- "M"
  paste(head(protein, -1L), collapse = "")
}

.breseq_gb_qualifier <- function(key, value) {
  value <- gsub('"', "'", as.character(value), fixed = TRUE)
  chunks <- character(); remaining <- value
  repeat {
    width <- if (!length(chunks)) max(1L, 76L - 21L - nchar(key)) else 58L
    cut <- min(width, nchar(remaining))
    if (key != "translation" && nchar(remaining) > cut) {
      spaces <- gregexpr(" ", substr(remaining, 1L, cut), fixed = TRUE)[[1]]
      if (any(spaces > 0)) cut <- max(spaces)
    }
    chunks <- c(chunks, substr(remaining, 1L, cut))
    remaining <- substring(remaining, cut + 1L)
    if (key != "translation") remaining <- sub("^ +", "", remaining)
    if (!nzchar(remaining)) break
  }
  lines <- c(paste0(strrep(" ", 21), "/", key, '="', chunks[1]),
             if (length(chunks) > 1L) paste0(strrep(" ", 21), chunks[-1]) else character())
  lines[length(lines)] <- paste0(lines[length(lines)], '"')
  lines
}

.breseq_gb_feature_lines <- function(key, location) {
  starts <- seq.int(1L, nchar(location), by = 59L)
  chunks <- substring(location, starts, starts + 58L)
  c(paste0("     ", sprintf("%-16s", key), chunks[1]),
    if (length(chunks) > 1L) paste0(strrep(" ", 21), chunks[-1]) else character())
}

.breseq_gb_features <- function(lines) {
  # Some complete reference records have CONTIG/BASE COUNT footer fields even
  # when ORIGIN bases are present. They are not feature qualifier continuations.
  footer <- grep("^[^ ]", lines)
  if (length(footer)) lines <- lines[seq_len(footer[1] - 1L)]
  starts <- grep("^ {5}[^ ]", lines)
  if (!length(starts)) return(list())
  ends <- c(starts[-1] - 1L, length(lines))
  lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    q <- grep("^ {21}/", block)
    location_end <- if (length(q)) q[1] - 1L else length(block)
    qualifiers <- if (length(q)) lapply(seq_along(q), function(j) {
      last <- if (j < length(q)) q[j + 1L] - 1L else length(block)
      block[q[j]:last]
    }) else list()
    list(key = trimws(substr(block[1], 6L, 20L)),
      location = paste(trimws(substring(block[seq_len(location_end)], 22L)), collapse = ""),
      qualifiers = qualifiers)
  })
}

.breseq_gb_qualifier_key <- function(block) sub("[=[:space:]].*$", "", sub("^ */", "", block[1]))
.breseq_gb_qualifier_value <- function(block) {
  value <- paste(trimws(block), collapse = "")
  sub('"$', "", sub('^[^=]*="?', "", value))
}

.breseq_gb_read_records <- function(reference_file) {
  lines <- readLines(reference_file, warn = FALSE)
  if (any(grepl("^LOCUS", lines))) {
    starts <- grep("^LOCUS", lines); ends <- grep("^//", lines)
    if (length(starts) != length(ends)) stop("Incomplete GenBank reference records")
    return(lapply(seq_along(starts), function(i) {
      record <- lines[starts[i]:ends[i]]
      origin <- grep("^ORIGIN", record)[1]; features <- grep("^FEATURES", record)[1]
      if (is.na(origin)) stop("Reference record has no ORIGIN sequence")
      id <- strsplit(trimws(record[1]), "[[:space:]]+")[[1]][2]
      version <- grep("^VERSION[[:space:]]", record, value = TRUE)
      accession <- grep("^ACCESSION[[:space:]]", record, value = TRUE)
      if (length(version)) id <- strsplit(trimws(version[1]), "[[:space:]]+")[[1]][2] else
        if (length(accession)) id <- strsplit(trimws(accession[1]), "[[:space:]]+")[[1]][2]
      sequence <- toupper(gsub("[^A-Za-z]", "", paste(record[(origin + 1L):(length(record) - 1L)], collapse = "")))
      if (!nzchar(sequence) || grepl("[^ACGTRYMKSWBDHVN]", sequence)) stop("Invalid reference DNA bases")
      list(seq_id = .breseq_seq(id), sequence = sequence,
        topology = if (grepl("circular", record[1], fixed = TRUE)) "circular" else "linear",
        header = record[seq_len(if (is.na(features)) origin - 1L else features - 1L)],
        features = if (is.na(features)) list() else .breseq_gb_features(record[(features + 1L):(origin - 1L)]))
    }))
  }
  starts <- grep("^>", lines)
  if (!length(starts)) stop("Reference must be a GenBank or FASTA file")
  ends <- c(starts[-1] - 1L, length(lines))
  lapply(seq_along(starts), function(i) {
    id <- strsplit(sub("^>", "", lines[starts[i]]), "[[:space:]]+")[[1]][1]
    sequence <- toupper(gsub("[[:space:]]", "", paste(lines[(starts[i] + 1L):ends[i]], collapse = "")))
    if (!nzchar(sequence) || grepl("[^ACGTRYMKSWBDHVN]", sequence)) stop("Invalid reference DNA bases")
    topology <- if (grepl("linear", lines[starts[i]], ignore.case = TRUE)) "linear" else "circular"
    list(seq_id = .breseq_seq(id), sequence = sequence, topology = topology,
         header = c(paste0("DEFINITION  ", sub("^>", "", lines[starts[i]]))), features = list())
  })
}

.breseq_gb_lift_features <- function(record, edited_sequence, edits) {
  output <- character(); audit <- list(); length_new <- nchar(edited_sequence)
  features <- record$features
  if (!any(vapply(features, function(x) x$key == "source", logical(1))))
    features <- c(list(list(key = "source", location = paste0("1..", nchar(record$sequence)), qualifiers = list())), features)
  for (i in seq_along(features)) {
    f <- features[[i]]; q <- f$qualifiers
    keys <- vapply(q, .breseq_gb_qualifier_key, "")
    locus <- which(keys %in% c("locus_tag", "gene"))
    label <- if (length(locus)) .breseq_gb_qualifier_value(q[[locus[1]]]) else ""
    old_location <- f$location
    result <- tryCatch({
      node <- .breseq_gb_location(old_location)
      if (f$key == "source") {
        node <- .breseq_gb_location(paste0("1..", length_new))
        list(node = node, changed = FALSE, partial = FALSE)
      } else .breseq_gb_lift_location(node, edits, nchar(record$sequence), length_new)
    }, error = function(e) list(node = NULL, changed = FALSE, partial = FALSE, error = conditionMessage(e)))
    status <- "preserved"; note <- ""; location <- ""
    if (is.null(result$node)) {
      status <- if (!is.null(result$error)) "dropped_unsupported_location" else "dropped_fully_deleted"
      note <- if (!is.null(result$error)) result$error else "All supported location segments were deleted"
    } else {
      location <- .breseq_gb_location_text(result$node)
      if (location != old_location) status <- "lifted"
      if (result$changed) status <- if (result$partial) "partial_sequence_change" else "sequence_changed"
      # Positional qualifiers must not retain old coordinates after liftover.
      positional_failure <- FALSE
      for (j in which(keys %in% c("transl_except", "anticodon", "rpt_unit_range"))) {
        value <- .breseq_gb_qualifier_value(q[[j]])
        old <- if (keys[j] == "rpt_unit_range") value else sub(",aa:.*$", "", sub("^\\(pos:", "", value))
        lifted <- tryCatch(.breseq_gb_lift_location(.breseq_gb_location(old), edits,
          nchar(record$sequence), length_new), error = function(e) NULL)
        if (is.null(lifted) || is.null(lifted$node) || lifted$changed) {
          q[[j]] <- character(); positional_failure <- TRUE
        } else {
          replacement <- .breseq_gb_location_text(lifted$node)
          q[[j]] <- .breseq_gb_qualifier(keys[j], sub(old, replacement, value, fixed = TRUE))
        }
      }
      if (f$key == "CDS" && (result$changed || positional_failure)) {
        code_at <- which(keys == "transl_table"); start_at <- which(keys == "codon_start")
        code <- if (length(code_at)) as.integer(.breseq_gb_qualifier_value(f$qualifiers[[code_at[1]]])) else 1L
        offset <- if (length(start_at)) as.integer(.breseq_gb_qualifier_value(f$qualifiers[[start_at[1]]])) else 1L
        translated <- tryCatch({
          if (result$partial || positional_failure || any(keys %in% c("pseudo", "pseudogene", "transl_except"))) stop("partial or exceptional CDS requires reannotation")
          .breseq_gb_translation(result$node, edited_sequence, code, offset)
        }, error = function(e) e)
        q <- q[!keys %in% c("translation", "protein_id", "inference")]
        if (inherits(translated, "error")) {
          if (!any(keys %in% c("pseudo", "pseudogene"))) q <- c(q, list(paste0(strrep(" ", 21), "/pseudo")))
          note <- paste0("Altered CDS; translation not asserted: ", conditionMessage(translated),
                         ". Reference annotation retained conservatively; functional status unverified.")
          status <- "CDS_requires_reannotation"
        } else {
          q <- c(q, list(.breseq_gb_qualifier("translation", translated)))
          note <- "CDS translation recomputed from the called-variant sequence; reference protein_id removed."
          status <- "CDS_translation_recomputed"
        }
      } else if (result$changed || positional_failure) {
        note <- "Reference feature lifted across called sequence changes; altered or truncated annotation requires review."
      }
      if (nzchar(note)) q <- c(q, list(.breseq_gb_qualifier("note", note)))
      output <- c(output, .breseq_gb_feature_lines(f$key, location), unlist(q, use.names = FALSE))
    }
    audit[[i]] <- data.frame(feature_key = f$key, locus = label, original_location = old_location,
      new_location = location, status = status, note = note, stringsAsFactors = FALSE)
  }
  list(lines = output, audit = do.call(rbind, audit))
}

.breseq_gb_header <- function(record, id, sample, length_new) {
  header <- record$header
  # Reference accession IDs must not identify a newly reconstructed sequence.
  starts <- grep("^[A-Z]", header)
  ends <- c(starts[-1] - 1L, length(header))
  keep <- rep(TRUE, length(header))
  for (i in seq_along(starts)) if (trimws(substr(header[starts[i]], 1, 12)) %in%
      c("LOCUS", "ACCESSION", "VERSION", "BASE COUNT", "CONTIG")) keep[starts[i]:ends[i]] <- FALSE
  header <- header[keep]
  date <- Sys.Date(); months <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
  date <- paste(format(date, "%d"), months[as.integer(format(date, "%m"))], format(date, "%Y"), sep = "-")
  locus <- sprintf("LOCUS       %-16s %11d bp    DNA     %-8s BCT %s", id, length_new, record$topology, date)
  # Insert local identifiers after DEFINITION, ahead of source/reference metadata.
  at <- grep("^[A-Z]", header)
  after_definition <- if (length(at) > 1L) at[2] - 1L else length(header)
  header <- append(header, c(paste0("ACCESSION   ", id), paste0("VERSION     ", id)), after_definition)
  header <- sub(paste0("(bases 1 to ", nchar(record$sequence), ")"),
                paste0("(bases 1 to ", length_new, ")"), header, fixed = TRUE)
  c(locus, header,
    paste0("COMMENT     Called-variant sequence for sample ", sample, "; reference ", record$seq_id, "."),
    "            Uncertain missing-coverage intervals are not applied as exact edits.")
}

.breseq_gb_origin <- function(sequence) {
  starts <- seq.int(1L, nchar(sequence), by = 60L)
  vapply(starts, function(start) {
    part <- tolower(substr(sequence, start, start + 59L))
    chunks <- substring(part, seq.int(1L, nchar(part), by = 10L), seq.int(1L, nchar(part), by = 10L) + 9L)
    paste0(sprintf("%9d", start), " ", paste(chunks, collapse = " "))
  }, "")
}

.breseq_gb_edit_features <- function(edits, length_new, sample) {
  output <- character()
  for (i in seq_len(nrow(edits))) {
    edit <- edits[i, ]
    if (edit$alt_length > 0) {
      location <- if (edit$new_start == edit$new_end) as.character(edit$new_start) else
        paste0(edit$new_start, "..", edit$new_end)
      note <- if (edit$ref_length == 0) paste0("Insertion after reference position ", edit$ref_end, ".") else
        paste0("Reference interval ", edit$ref_start, "..", edit$ref_end, ".")
    } else {
      location <- as.character(max(1, min(length_new, edit$new_start - 1)))
      note <- paste0("Deletion junction marker on an adjacent remaining base; deleted reference interval ",
        edit$ref_start, "..", edit$ref_end, ". No deleted bases remain at this marker.")
    }
    note <- paste0(note, " Sample ", sample, "; applied ", edit$mutation,
                   "; reference length ", edit$ref_length, "; alternate length ", edit$alt_length, ".")
    color <- unname(.get_genbank_mutation_color(edit$type))
    output <- c(output, .breseq_gb_feature_lines("misc_feature", location),
      .breseq_gb_qualifier("label", paste0("[", edit$type, ":", sample, "] ", edit$mutation)),
      .breseq_gb_qualifier("note", note), .breseq_gb_qualifier("ApEinfo_fwdcolor", color),
      .breseq_gb_qualifier("ApEinfo_revcolor", color))
  }
  output
}

#' Export GenBank Files Containing Actually Reconstructed Sample Sequences
#'
#' Applies only final mutation calls, lifts reference features through sequence
#' edits, and recomputes safe table 1/11 CDS translations. Display-only coverage
#' evidence is excluded. Full plasmid losses have manifest entries without files.
#' @param reference_file GenBank or FASTA reference, including multiple records.
#' @param mutation_list Named processed per-sample final mutation tables.
#' @param output_dir Directory for per-sample, per-replicon GenBank and FASTA files.
#' @param samples Sample IDs in output order; a missing WT entry is an empty baseline.
#' @param sample_labels Optional named display labels for sample IDs.
#' @param replicon_losses Optional validated whole-replicon coverage-loss records
#'   with sample, seq_id, evidence_origin, and depth_ratio columns.
#' @return List with manifest, feature_audit, edits, and call_report tables.
#' @export
export_changed_sequences_to_genbank <- function(reference_file, mutation_list, output_dir,
                                                samples = names(mutation_list), sample_labels = NULL,
                                                replicon_losses = NULL) {
  records <- .breseq_gb_read_records(reference_file)
  ids <- vapply(records, `[[`, "", "seq_id")
  if (anyDuplicated(ids)) stop("Reference replicon IDs must be unique")
  if (!is.list(mutation_list) || !is.character(samples) || !length(samples) || anyNA(samples) ||
      any(!nzchar(samples)) || anyDuplicated(samples) || any(!samples %in% c(names(mutation_list), "WT")))
    stop("samples must identify unique mutation_list stages (or an empty WT baseline)")
  if (!is.null(replicon_losses)) {
    required <- c("sample", "seq_id", "evidence_origin", "depth_ratio")
    if (!all(required %in% names(replicon_losses))) stop("replicon_losses is missing required evidence columns")
    ratio <- suppressWarnings(as.numeric(replicon_losses$depth_ratio))
    if (any(!replicon_losses$sample %in% samples) || any(!.breseq_seq(replicon_losses$seq_id) %in% ids) ||
        anyNA(replicon_losses$evidence_origin) || any(replicon_losses$evidence_origin != "COVERAGE_RATIO") ||
        any(!is.finite(ratio) | ratio < 0 | ratio >= 0.1) ||
        anyDuplicated(paste(replicon_losses$sample, .breseq_seq(replicon_losses$seq_id))))
      stop("replicon_losses must contain unique, valid whole-replicon COVERAGE_RATIO evidence below 0.1")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  manifests <- list(); audits <- list(); all_edits <- list(); reports <- list()
  for (s in seq_along(samples)) for (r in seq_along(records)) {
    sample <- samples[s]; record <- records[[r]]; data <- mutation_list[[sample]]
    if (is.null(data)) data <- data.frame()
    calls <- data
    if (nrow(calls)) {
      names(calls) <- .breseq_text(names(calls))
      seq_col <- intersect(c("seq id", "seq_id", "seqid"), names(calls))
      if (!length(seq_col)) stop("Final mutation rows need a replicon ID column")
      if (any(!.breseq_seq(calls[[seq_col[1]]]) %in% ids)) stop("Mutation rows contain an unknown replicon ID")
      keep <- .breseq_seq(calls[[seq_col[1]]]) == record$seq_id
      if ("annotation" %in% names(calls)) keep <- keep &
        !(tolower(.breseq_text(calls$annotation)) %in% "coverage-derived del display")
      if ("evidence" %in% names(calls)) keep <- keep & !(calls$evidence %in% "COVERAGE_RATIO")
      calls <- calls[keep, , drop = FALSE]
    }
    loss <- !is.null(replicon_losses) && any(replicon_losses$sample == sample &
      .breseq_seq(replicon_losses$seq_id) == record$seq_id)
    result <- if (loss) NULL else tryCatch(.breseq_apply_sequence_edits(record$sequence, calls), error = function(e) e)
    status <- if (loss) "lost_by_coverage" else if (inherits(result, "error")) "unresolved_calls" else
      if (!nchar(result$sequence)) "lost_confirmed" else if (!nrow(calls)) "reference_baseline" else "exported_called_variants"
    filename <- fasta <- checksum <- ""; new_length <- if (loss) 0 else if (inherits(result, "error")) NA_real_ else nchar(result$sequence)
    label <- if (!is.null(sample_labels) && sample %in% names(sample_labels)) unname(sample_labels[sample]) else sample
    if (!loss && !inherits(result, "error")) {
      edits <- result$edits
      if (nrow(edits)) {
        if ("source_rows" %in% names(edits) && is.list(edits$source_rows)) edits$source_rows <- vapply(edits$source_rows, paste, "", collapse = ",")
        all_edits[[length(all_edits) + 1L]] <- cbind(sample = sample, seq_id = record$seq_id, edits, stringsAsFactors = FALSE)
      }
      if (!is.null(result$report) && nrow(result$report)) reports[[length(reports) + 1L]] <-
        cbind(sample = sample, seq_id = record$seq_id, result$report, stringsAsFactors = FALSE)
      if (new_length > 0) {
        lifted <- .breseq_gb_lift_features(record, result$sequence, edits)
        audits[[length(audits) + 1L]] <- cbind(sample = sample, seq_id = record$seq_id, lifted$audit, stringsAsFactors = FALSE)
        stem <- paste0(gsub("[^A-Za-z0-9_.-]", "_", sample), "_", record$seq_id, "_called_variants")
        filename <- file.path(output_dir, paste0(stem, ".gbk")); fasta <- file.path(output_dir, paste0(stem, ".fasta"))
        local_id <- sprintf("BRS_%03d_%03d", s, r)
        writeLines(c(.breseq_gb_header(record, local_id, label, new_length), "FEATURES             Location/Qualifiers",
          lifted$lines, .breseq_gb_edit_features(edits, new_length, sample),
          "ORIGIN", .breseq_gb_origin(result$sequence), "//"), filename, useBytes = TRUE)
        starts <- seq.int(1L, new_length, by = 80L)
        writeLines(c(paste0(">", local_id, " ", sample, " ", record$seq_id, " called variants"),
          substring(result$sequence, starts, starts + 79L)), fasta, useBytes = TRUE)
        hash_file <- tempfile(); writeChar(result$sequence, hash_file, eos = NULL, useBytes = TRUE)
        checksum <- unname(tools::md5sum(hash_file)); unlink(hash_file)
      }
    } else if (inherits(result, "error") && !is.null(result$report) && nrow(result$report))
      reports[[length(reports) + 1L]] <- cbind(sample = sample, seq_id = record$seq_id, result$report, stringsAsFactors = FALSE)
    manifests[[length(manifests) + 1L]] <- data.frame(sample = sample, sample_label = label, seq_id = record$seq_id,
      reference_length = nchar(record$sequence), new_length = new_length, sequence_length = new_length,
      final_call_rows = nrow(calls), status = status, file = filename, genbank_file = filename, fasta_file = fasta,
      sequence_md5 = checksum, message = if (inherits(result, "error")) conditionMessage(result) else "", stringsAsFactors = FALSE)
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
  list(manifest = bind(manifests), feature_audit = bind(audits), edits = bind(all_edits), call_report = bind(reports))
}
