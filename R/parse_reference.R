#' Parse Reference Genome Information
#'
#' Extracts replicon names, lengths, and topology (linear/circular) from
#' GenBank, FASTA, or SnapGene .dna files. Supports multi-replicon genomes
#' (chromosome + multiple plasmids).
#'
#' @param reference_file Path to reference file (.gbk, .gb, .gbff, .fasta, .fa, .fna, .dna)
#' @return Data frame with columns: seq_id, length, topology, description
#' @export
#' @examples
#' \dontrun{
#' # Get info from GenBank file
#' info <- parse_reference_info("reference.gbk")
#' # Returns:
#' #   seq_id    length topology description
#' # 1 chromosome 4600000 circular E. coli K-12
#' # 2 pBR322       4361 circular Cloning vector
#' # 3 pUC19        2686 circular Cloning vector
#'
#' # Use for normalization
#' genome_lengths <- setNames(info$length, info$seq_id)
#' }
parse_reference_info <- function(reference_file) {

  if (!file.exists(reference_file)) {
    stop("Reference file not found: ", reference_file)
  }

  ext <- tolower(tools::file_ext(reference_file))

  result <- switch(ext,
    "gbk" = , "gb" = , "gbff" = , "genbank" = parse_genbank_info(reference_file),
    "fasta" = , "fa" = , "fna" = , "faa" = parse_fasta_info(reference_file),
    "dna" = parse_dna_info(reference_file),
    stop("Unsupported file format: ", ext,
         "\nSupported: .gbk, .gb, .gbff, .fasta, .fa, .fna, .dna")
  )

  message("Parsed ", nrow(result), " replicon(s) from ", basename(reference_file))
  for (i in seq_len(nrow(result))) {
    message("  ", result$seq_id[i], ": ",
            format(result$length[i], big.mark = ","), " bp (",
            result$topology[i], ")")
  }

  return(result)
}

#' Parse GenBank File for Replicon Information
#'
#' Reads GenBank file(s) and extracts LOCUS information for each replicon.
#'
#' @param gbk_file Path to GenBank file
#' @return Data frame with seq_id, length, topology, description
#' @export
parse_genbank_info <- function(gbk_file) {

  lines <- readLines(gbk_file, warn = FALSE)

  # Find all LOCUS lines
  locus_idx <- grep("^LOCUS", lines)

  if (length(locus_idx) == 0) {
    stop("No LOCUS line found in GenBank file")
  }

  results <- lapply(locus_idx, function(idx) {
    locus_line <- lines[idx]

    # Parse LOCUS line
    # Format: LOCUS       name    length bp    mol_type    topology    division date
    # Example: LOCUS       NC_000913  4641652 bp    DNA     circular BCT 17-SEP-2019

    # Extract sequence ID (first field after LOCUS)
    parts <- strsplit(trimws(locus_line), "\\s+")[[1]]
    seq_id <- parts[2]

    # Extract length (number followed by "bp")
    bp_match <- regmatches(locus_line, regexpr("\\d+\\s*bp", locus_line))
    length <- if (length(bp_match) > 0) {
      as.integer(gsub("\\s*bp", "", bp_match))
    } else NA_integer_

    # Extract topology (circular or linear)
    topology <- if (grepl("circular", locus_line, ignore.case = TRUE)) {
      "circular"
    } else if (grepl("linear", locus_line, ignore.case = TRUE)) {
      "linear"
    } else {
      "unknown"
    }

    # Find DEFINITION line for description
    def_idx <- idx + which(grepl("^DEFINITION", lines[(idx+1):min(idx+10, length(lines))]))[1]
    description <- if (!is.na(def_idx) && def_idx <= length(lines)) {
      sub("^DEFINITION\\s+", "", lines[def_idx])
    } else { "" }

    data.frame(
      seq_id = seq_id,
      length = length,
      topology = topology,
      description = description,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Parse FASTA File for Sequence Information
#'
#' Reads FASTA file and calculates sequence lengths.
#' Topology is inferred from header if possible, otherwise assumed circular for bacteria.
#'
#' @param fasta_file Path to FASTA file
#' @param default_topology Default topology if not specified ("circular" or "linear")
#' @return Data frame with seq_id, length, topology, description
#' @export
parse_fasta_info <- function(fasta_file, default_topology = "circular") {

  lines <- readLines(fasta_file, warn = FALSE)

  # Find header lines
  header_idx <- grep("^>", lines)

  if (length(header_idx) == 0) {
    stop("No FASTA headers found")
  }

  # Add end position for last sequence
  seq_ends <- c(header_idx[-1] - 1, length(lines))

  results <- lapply(seq_along(header_idx), function(i) {
    header <- lines[header_idx[i]]
    seq_lines <- lines[(header_idx[i] + 1):seq_ends[i]]

    # Parse header: >seq_id description
    header_clean <- sub("^>", "", header)
    parts <- strsplit(header_clean, "\\s+", perl = TRUE)[[1]]
    seq_id <- parts[1]
    description <- if (length(parts) > 1) paste(parts[-1], collapse = " ") else ""

    # Calculate sequence length
    seq_length <- sum(nchar(gsub("\\s", "", seq_lines)))

    # Infer topology from description
    topology <- if (grepl("circular|plasmid", description, ignore.case = TRUE)) {
      "circular"
    } else if (grepl("linear|chromosome", description, ignore.case = TRUE)) {
      # Note: bacterial chromosomes are usually circular
      if (grepl("plasmid", description, ignore.case = TRUE)) "circular" else default_topology
    } else {
      default_topology
    }

    data.frame(
      seq_id = seq_id,
      length = seq_length,
      topology = topology,
      description = description,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Parse SnapGene .dna File for Sequence Information
#'
#' Converts .dna to GenBank using SnapGene CLI, then parses the result.
#' Falls back to binary parsing if SnapGene is not available.
#'
#' @param dna_file Path to .dna file
#' @return Data frame with seq_id, length, topology, description
#' @export
parse_dna_info <- function(dna_file) {

  # Try to convert to GenBank first
  if (has_snapgene()) {
    temp_gbk <- tempfile(fileext = ".gbk")
    on.exit(unlink(temp_gbk), add = TRUE)

    tryCatch({
      convert_dna_to_genbank(dna_file, temp_gbk)
      return(parse_genbank_info(temp_gbk))
    }, error = function(e) {
      message("SnapGene conversion failed, trying binary parse: ", e$message)
    })
  }

  # Fallback: parse .dna binary format directly
  # SnapGene .dna format is documented at:
  # https://www.snapgene.com/resources/file-format
  parse_dna_binary(dna_file)
}

#' Parse SnapGene .dna Binary Format
#'
#' Directly parses the SnapGene binary format without needing SnapGene installed.
#'
#' @param dna_file Path to .dna file
#' @return Data frame with seq_id, length, topology, description
#' @keywords internal
parse_dna_binary <- function(dna_file) {

  # Read binary file
  con <- file(dna_file, "rb")
  on.exit(close(con), add = TRUE)

  # Read header
  header <- readBin(con, "raw", n = 1)

  # Check magic byte (0x09 for SnapGene)
  if (as.integer(header) != 9) {
    warning("File may not be a valid SnapGene .dna file")
  }

  # Read file structure
  # SnapGene uses a chunk-based format

  seq_length <- NA_integer_
  topology <- "circular"  # Default assumption
  seq_id <- tools::file_path_sans_ext(basename(dna_file))
  description <- ""

  # Try to read sequence length from file
  tryCatch({
    # Seek to find DNA sequence chunk
    seek(con, 0)
    raw_data <- readBin(con, "raw", n = file.info(dna_file)$size)

    # Look for sequence data (usually after header chunks)
    # The DNA sequence is stored as ASCII in the file
    raw_text <- rawToChar(raw_data[raw_data >= 0x20 & raw_data <= 0x7E])

    # Try to find sequence (ATCG only regions)
    dna_matches <- gregexpr("[ATCGatcg]{100,}", raw_text)[[1]]
    if (dna_matches[1] > 0) {
      # Find the longest match (likely the main sequence)
      match_lengths <- attr(dna_matches, "match.length")
      best_match <- which.max(match_lengths)
      seq_length <- match_lengths[best_match]
    }

    # Check for circular indicator
    if (grepl("circular|\\bcircular\\b", raw_text, ignore.case = TRUE)) {
      topology <- "circular"
    } else if (grepl("linear", raw_text, ignore.case = TRUE)) {
      topology <- "linear"
    }

  }, error = function(e) {
    warning("Could not parse .dna file structure: ", e$message)
  })

  data.frame(
    seq_id = seq_id,
    length = seq_length,
    topology = topology,
    description = description,
    stringsAsFactors = FALSE
  )
}

#' Get Genome Lengths as Named Vector
#'
#' Convenience function to get genome lengths in a format ready for use
#' with normalization functions.
#'
#' @param reference_file Path to reference file
#' @return Named numeric vector (names = seq_id, values = length)
#' @export
#' @examples
#' \dontrun{
#' genome_lengths <- get_genome_lengths("reference.gbk")
#' # chromosome      pBR322       pUC19
#' #    4600000        4361        2686
#' }
get_genome_lengths <- function(reference_file) {
  info <- parse_reference_info(reference_file)
  setNames(info$length, info$seq_id)
}

#' Get Replicon Topology as Named Vector
#'
#' @param reference_file Path to reference file
#' @return Named character vector (names = seq_id, values = topology)
#' @export
get_replicon_topology <- function(reference_file) {
  info <- parse_reference_info(reference_file)
  setNames(info$topology, info$seq_id)
}

#' Check if Replicon is Circular
#'
#' @param seq_id Sequence ID to check
#' @param reference_file Path to reference file
#' @return Logical
#' @export
is_circular <- function(seq_id, reference_file) {
  topology <- get_replicon_topology(reference_file)
  if (seq_id %in% names(topology)) {
    return(topology[seq_id] == "circular")
  }
  # Default: assume circular for bacterial genomes
  TRUE
}

#' Auto-detect seq_id Column in breseq Data
#'
#' breseq output includes a "seq id" column that identifies which replicon
#' (chromosome or plasmid) each mutation is on.
#'
#' @param df Data frame from breseq
#' @return Name of the seq_id column, or NULL if not found
#' @keywords internal
find_seq_id_column <- function(df) {
  candidates <- c("seq id", "seq_id", "seqid", "seq.id",
                  "replicon", "contig", "chromosome")

  for (col in candidates) {
    if (col %in% names(df)) return(col)
  }

  # Try case-insensitive match
  col_lower <- tolower(names(df))
  for (cand in candidates) {
    idx <- which(col_lower == tolower(cand))
    if (length(idx) > 0) return(names(df)[idx[1]])
  }

  NULL
}
