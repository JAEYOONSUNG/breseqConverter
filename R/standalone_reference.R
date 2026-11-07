# Prepare one reference per replicon for the standalone driver.  The GC reader
# consumes one record at a time, so the original multi-record file must not be
# reused for every replicon.  This helper deliberately does not parse features;
# the driver can supply separately parsed GFF annotations through gff_data.
#
# GenBank records retain their sequence and feature text.  Only the temporary
# LOCUS identifier is aligned to VERSION (ACCESSION/LOCUS fallback), with a final
# numeric accession version removed to match breseq/GFF sequence identifiers.
# FASTA records retain their sequence and description with the same ID rule.
# Explicit topology is retained.  FASTA without topology defaults to circular,
# as in parse_fasta_info(); GenBank without topology remains unknown.
.breseq_prepare_standalone_reference <- function(reference_file, work_dir) {
  scalar_path <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!scalar_path(reference_file) || !file.exists(reference_file) ||
      dir.exists(reference_file)) {
    stop("Reference file not found: ", paste(reference_file, collapse = ", "))
  }
  if (!scalar_path(work_dir)) stop("work_dir must be a nonempty directory path")

  extension <- tolower(tools::file_ext(reference_file))
  genbank <- extension %in% c("gb", "gbk", "gbff", "genbank")
  fasta <- extension %in% c("fa", "fna", "fasta")
  if (!genbank && !fasta) {
    stop("Standalone references must be GenBank (.gb, .gbk, .gbff, .genbank) ",
         "or nucleotide FASTA (.fa, .fna, .fasta).")
  }

  lines <- readLines(reference_file, warn = FALSE)
  starts <- grep(if (genbank) "^LOCUS[[:space:]]" else "^>", lines)
  if (!length(starts)) {
    stop(if (genbank) "No LOCUS records found in GenBank reference" else
           "No FASTA records found in reference")
  }
  ends <- c(starts[-1L] - 1L, length(lines))
  normalize_id <- function(id) sub("\\.[0-9]+$", "", id)
  first_field <- function(text) {
    fields <- strsplit(trimws(text), "[[:space:]]+")[[1L]]
    if (length(fields) && nzchar(fields[1L])) fields[1L] else ""
  }
  topology_from <- function(text, default) {
    if (grepl("\\bcircular\\b", text, ignore.case = TRUE)) "circular" else
      if (grepl("\\blinear\\b", text, ignore.case = TRUE)) "linear" else default
  }

  records <- lapply(seq_along(starts), function(i) {
    record <- lines[seq.int(starts[i], ends[i])]
    if (genbank) {
      terminator <- grep("^//[[:space:]]*$", record)
      if (length(terminator) != 1L) {
        stop("GenBank record ", i, " must have exactly one // terminator")
      }
      record <- record[seq_len(terminator)]
      origin <- grep("^ORIGIN([[:space:]]|$)", record)
      if (length(origin) != 1L || origin >= length(record) - 1L) {
        stop("GenBank record ", i, " has no usable ORIGIN sequence")
      }
      sequence <- toupper(gsub("[[:space:][:digit:]]", "",
        paste(record[seq.int(origin + 1L, length(record) - 1L)], collapse = "")))
      locus_id <- first_field(sub("^LOCUS[[:space:]]+", "", record[1L]))
      field_value <- function(field) {
        hit <- grep(paste0("^", field, "[[:space:]]+"), record)
        if (length(hit)) first_field(sub(paste0("^", field, "[[:space:]]+"),
                                        "", record[hit[1L]])) else ""
      }
      id_candidates <- c(field_value("VERSION"), field_value("ACCESSION"), locus_id)
      id <- normalize_id(id_candidates[nzchar(id_candidates)][1L])
      declared_match <- regmatches(record[1L],
        regexpr("[0-9]+[[:space:]]+bp\\b", record[1L], ignore.case = TRUE))
      declared_length <- if (length(declared_match)) {
        as.numeric(sub("[[:space:]]+bp$", "", declared_match, ignore.case = TRUE))
      } else NA_real_
      if (!is.finite(declared_length) || declared_length != nchar(sequence)) {
        stop("GenBank record ", i, " LOCUS length does not match its ORIGIN sequence")
      }
      topology <- topology_from(record[1L], "unknown")
      definition <- grep("^DEFINITION[[:space:]]", record)
      description <- ""
      if (length(definition)) {
        j <- definition[1L]
        description <- sub("^DEFINITION[[:space:]]+", "", record[j])
        while (j + 1L < origin && grepl("^ {12}[^ ]", record[j + 1L])) {
          j <- j + 1L
          description <- paste(description, trimws(record[j]))
        }
      }
      # Preserve all other record text, including feature locations and bases.
      record[1L] <- paste0("LOCUS       ", id,
        sub("^LOCUS[[:space:]]+[^[:space:]]+", "", record[1L]))
    } else {
      header <- trimws(sub("^>", "", record[1L]))
      raw_id <- first_field(header)
      id <- normalize_id(raw_id)
      description <- trimws(sub("^[^[:space:]]+", "", header))
      if (length(record) < 2L) stop("FASTA record ", i, " has no sequence")
      sequence <- toupper(gsub("[[:space:]]", "", paste(record[-1L], collapse = "")))
      topology <- topology_from(description, "circular")
      # Normalizing whitespace does not remove ambiguous nucleotide symbols.
      chunk_starts <- seq.int(1L, max(1L, nchar(sequence)), by = 80L)
      record <- c(paste0(">", id, if (nzchar(description)) paste0(" ", description)),
                  substring(sequence, chunk_starts, chunk_starts + 79L))
    }
    if (is.na(id) || !nzchar(id)) stop("Reference record ", i, " has no sequence ID")
    if (!nzchar(sequence) || grepl("[^ACGTRYSWKMBDHVN]", sequence)) {
      stop("Reference record ", i, " must contain a nonempty nucleotide sequence")
    }
    list(lines = record, info = data.frame(seq_id = id,
      length = nchar(sequence), topology = topology, description = description,
      stringsAsFactors = FALSE))
  })

  info <- do.call(rbind, lapply(records, `[[`, "info"))
  rownames(info) <- NULL
  if (anyDuplicated(info$seq_id)) {
    stop("Reference has duplicate sequence IDs after removing accession versions: ",
         paste(unique(info$seq_id[duplicated(info$seq_id)]), collapse = ", "))
  }
  if (!dir.exists(work_dir) && !dir.create(work_dir, recursive = TRUE)) {
    stop("Could not create reference work directory: ", work_dir)
  }
  work_dir <- normalizePath(work_dir, mustWork = TRUE)
  written <- character()
  complete <- FALSE
  on.exit(if (!complete) unlink(written), add = TRUE)
  replicons <- lapply(seq_along(records), function(i) {
    safe_id <- gsub("[^A-Za-z0-9_.-]", "_", info$seq_id[i])
    destination <- tempfile(paste0(sprintf("%03d", i), "_", safe_id, "_"),
                            tmpdir = work_dir, fileext = if (genbank) ".gbff" else ".fasta")
    written <<- c(written, destination)
    writeLines(records[[i]]$lines, destination, useBytes = TRUE)
    list(reference_file = destination, genome_info = info[i, , drop = FALSE],
         gff_data = NULL)
  })
  names(replicons) <- info$seq_id
  complete <- TRUE
  list(info = info, replicons = replicons)
}
