#' Normalize Genome Position for Circular Sequences
#'
#' Handles position normalization for circular genomes (plasmids, bacterial chromosomes).
#' Ensures positions wrap around correctly and handles edge cases near origin.
#'
#' @param position Numeric position(s) to normalize
#' @param genome_length Total length of the circular genome
#' @param origin Position of the origin (default = 1)
#' @return Normalized position(s) within 1:genome_length
#' @export
#' @examples
#' # Position beyond genome length wraps around
#' normalize_position(5100, genome_length = 5000)  # Returns 100
#'
#' # Negative positions wrap from end
#' normalize_position(-50, genome_length = 5000)   # Returns 4950
normalize_position <- function(position, genome_length, origin = 1) {

  if (is.na(position) || is.null(genome_length) || genome_length <= 0) {
    return(NA_integer_)
  }

  # Handle numeric conversion
  pos <- as.numeric(gsub(",", "", as.character(position)))

  if (is.na(pos)) return(NA_integer_)

  # Normalize to 1-based coordinates
  normalized <- ((pos - origin) %% genome_length) + 1

  # Ensure positive
  if (normalized <= 0) {
    normalized <- normalized + genome_length
  }

  return(as.integer(normalized))
}

#' Normalize Position Range for Circular Sequences
#'
#' Normalizes a start-end range for circular genomes. Handles cases where
#' the feature spans the origin.
#'
#' @param start Start position
#' @param end End position
#' @param genome_length Total length of the circular genome
#' @return List with normalized start, end, and spans_origin flag
#' @export
#' @examples
#' # Normal range
#' normalize_range(100, 500, genome_length = 5000)
#'
#' # Range spanning origin
#' normalize_range(4900, 100, genome_length = 5000)
normalize_range <- function(start, end, genome_length) {

  start_norm <- normalize_position(start, genome_length)
  end_norm <- normalize_position(end, genome_length)

  # Check if range spans the origin
  spans_origin <- FALSE
  if (!is.na(start_norm) && !is.na(end_norm)) {
    if (start_norm > end_norm) {
      spans_origin <- TRUE
    }
  }

  list(
    start = start_norm,
    end = end_norm,
    spans_origin = spans_origin,
    # For features spanning origin, provide both segments
    segments = if (spans_origin) {
      list(
        c(start_norm, genome_length),
        c(1, end_norm)
      )
    } else {
      list(c(start_norm, end_norm))
    }
  )
}

#' Parse Position String from breseq Output
#'
#' Parses various position formats from breseq HTML output.
#' Handles formats like: "1234", "1,234", "1234-5678", "1234–5678" (en-dash)
#'
#' @param position_str Position string to parse
#' @return List with start and end positions, or single position if not a range
#' @export
parse_position <- function(position_str) {

  if (is.na(position_str) || is.null(position_str)) {
    return(list(start = NA_integer_, end = NA_integer_))
  }

  # Convert to string and clean
  pos_str <- as.character(position_str)
  pos_str <- gsub(",", "", pos_str)           # Remove commas
  pos_str <- gsub("–", "-", pos_str)          # Replace en-dash with hyphen
  pos_str <- gsub("—", "-", pos_str)          # Replace em-dash with hyphen
  pos_str <- trimws(pos_str)
  pos_str <- gsub("-$", "", pos_str)          # Remove trailing dash

  # Check for range format
  if (grepl("-", pos_str)) {
    parts <- strsplit(pos_str, "-")[[1]]
    if (length(parts) == 2) {
      start <- as.integer(parts[1])
      end <- as.integer(parts[2])
      return(list(start = start, end = end))
    }
  }

  # Single position
  pos <- as.integer(pos_str)
  return(list(start = pos, end = pos))
}

#' Calculate Distance Between Positions on Circular Genome
#'
#' Calculates the shortest distance between two positions on a circular genome.
#'
#' @param pos1 First position
#' @param pos2 Second position
#' @param genome_length Total genome length
#' @return Shortest distance (always positive)
#' @export
circular_distance <- function(pos1, pos2, genome_length) {

  if (any(is.na(c(pos1, pos2, genome_length)))) return(NA_real_)

  direct <- abs(pos2 - pos1)
  wraparound <- genome_length - direct

  min(direct, wraparound)
}

#' Merge Overlapping Mutations on Circular Genome
#'
#' Identifies and optionally merges mutations that are at the same position
#' or overlapping, accounting for circular genome topology.
#'
#' @param df Data frame with position column
#' @param position_col Name of position column
#' @param genome_length Optional genome length for circular handling
#' @param merge_distance Maximum distance to consider as same position (default = 0)
#' @param seq_id_col Optional column for sequence ID (for multi-replicon genomes)
#' @return Data frame with added 'position_group' column for identifying duplicates
#' @export
identify_position_duplicates <- function(df,
                                         position_col = "position",
                                         genome_length = NULL,
                                         merge_distance = 0,
                                         seq_id_col = NULL) {

  if (!position_col %in% names(df)) {
    warning("Position column '", position_col, "' not found.")
    df$position_group <- seq_len(nrow(df))
    return(df)
  }

  # Parse positions
  parsed <- lapply(df[[position_col]], parse_position)
  df$parsed_start <- sapply(parsed, `[[`, "start")

  # If sequence ID column exists, group by it first
  if (!is.null(seq_id_col) && seq_id_col %in% names(df)) {
    df$seq_key <- paste(df[[seq_id_col]], df$parsed_start, sep = "_")
  } else {
    df$seq_key <- as.character(df$parsed_start)
  }

  # Assign position groups
  if (is.null(genome_length) || merge_distance == 0) {
    # Simple exact matching
    df$position_group <- as.integer(factor(df$seq_key))
  } else {
    # Consider circular distance for merging
    positions <- df$parsed_start
    n <- length(positions)
    groups <- rep(0, n)
    current_group <- 0

    for (i in seq_len(n)) {
      if (groups[i] == 0) {
        current_group <- current_group + 1
        groups[i] <- current_group

        for (j in seq_len(n)[-i]) {
          if (groups[j] == 0) {
            dist <- circular_distance(positions[i], positions[j], genome_length)
            if (!is.na(dist) && dist <= merge_distance) {
              groups[j] <- current_group
            }
          }
        }
      }
    }
    df$position_group <- groups
  }

  # Clean up temporary columns
  df$parsed_start <- NULL
  df$seq_key <- NULL

  # Add duplicate flag
  dup_counts <- table(df$position_group)
  df$is_duplicate <- df$position_group %in% names(dup_counts[dup_counts > 1])

  return(df)
}

#' Normalize Positions in Data Frame
#'
#' Applies position normalization to an entire data frame, handling
#' multiple replicons (chromosome + plasmids) if present.
#' Can auto-detect genome lengths from reference file.
#'
#' @param df Data frame with position data
#' @param position_col Name of position column(s)
#' @param genome_lengths Named vector of genome lengths (names = seq_id), OR
#' @param reference_file Path to reference file (.gbk, .fasta, .dna) to auto-parse lengths
#' @param seq_id_col Column containing sequence ID (auto-detected if NULL)
#' @return Data frame with normalized positions
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect from GenBank (recommended)
#' df <- normalize_positions_df(df, reference_file = "reference.gbk")
#'
#' # Or provide lengths manually
#' df <- normalize_positions_df(df, genome_lengths = c(chromosome = 4600000, pBR322 = 4361))
#' }
normalize_positions_df <- function(df,
                                   position_col = "position",
                                   genome_lengths = NULL,
                                   reference_file = NULL,
                                   seq_id_col = NULL) {

  # Auto-parse genome lengths from reference file
  if (is.null(genome_lengths) && !is.null(reference_file)) {
    genome_lengths <- get_genome_lengths(reference_file)
  }

  if (is.null(genome_lengths)) {
    message("No genome lengths provided. Positions not normalized for circularity.")
    return(df)
  }

  # Auto-detect seq_id column
  if (is.null(seq_id_col)) {
    seq_id_col <- find_seq_id_column(df)
  }

  # Normalize each position based on its replicon
  for (col in position_col) {
    if (!col %in% names(df)) next

    new_col <- paste0(col, "_normalized")

    if (!is.null(seq_id_col) && seq_id_col %in% names(df)) {
      df[[new_col]] <- mapply(function(pos, sid) {
        glen <- genome_lengths[sid]
        if (is.na(glen)) return(pos)
        normalize_position(pos, glen)
      }, df[[col]], df[[seq_id_col]])
    } else if (length(genome_lengths) == 1) {
      df[[new_col]] <- sapply(df[[col]], normalize_position,
                               genome_length = genome_lengths[1])
    } else {
      warning("Multiple replicons detected but no seq_id column found. ",
              "Using first genome length for all positions.")
      df[[new_col]] <- sapply(df[[col]], normalize_position,
                               genome_length = genome_lengths[1])
    }
  }

  return(df)
}
