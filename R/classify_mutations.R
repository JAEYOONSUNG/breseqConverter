#' Classify Mutation Types
#'
#' Analyzes mutation strings from breseq output and classifies them into
#' standardized categories: SNP, insertion, deletion, substitution, amplification, inversion, or mobile element.
#'
#' @param mutation_string Character string describing the mutation (e.g., "A→G", "Δ1 bp", "+GC")
#' @param annotation Optional annotation string for additional context
#' @return Character string indicating mutation type
#' @export
#' @examples
#' classify_mutation("A→G")           # "SNP"
#' classify_mutation("Δ500 bp")       # "deletion"
#' classify_mutation("+ATCG")         # "insertion"
#' classify_mutation("(G)8→9")        # "insertion"  (repeat expansion)
#' classify_mutation("(G)9→8")        # "deletion"   (repeat contraction)
classify_mutation <- function(mutation_string, annotation = NULL) {

  if (is.na(mutation_string) || is.null(mutation_string) ||
      tolower(mutation_string) %in% c("na", "nan", "unknown", "")) {
    return("unknown")
  }

  mut <- as.character(mutation_string)

  # Mobile element insertion (IS element)
  if (grepl("IS\\d+|mobile.element|transpos", mut, ignore.case = TRUE)) {
    return("mobile_element")
  }

  # Large deletion (Δ notation)
  if (grepl("^Δ|Δ\\d+", mut)) {
    # Extract size if available
    size_match <- regmatches(mut, regexpr("\\d+", mut))
    if (length(size_match) > 0) {
      size <- as.numeric(gsub(",", "", size_match[1]))
      if (size >= 50) return("large_deletion")
    }
    return("deletion")
  }

  # Simple insertion (+bases)
  if (grepl("^\\+[ATCG]+$", mut, ignore.case = TRUE)) {
    inserted <- gsub("^\\+", "", mut)
    if (nchar(inserted) >= 50) return("large_insertion")
    return("insertion")
  }

  # Simple deletion (-bases)
  if (grepl("^-[ATCG]+$", mut, ignore.case = TRUE)) {
    deleted <- gsub("^-", "", mut)
    if (nchar(deleted) >= 50) return("large_deletion")
    return("deletion")
  }

  # Repeat expansion/contraction: (G)8→9 or (ATCG)5→3
  if (grepl("\\([ATCG]+\\)\\d+.*→.*\\d+", mut, ignore.case = TRUE)) {
    # Extract numbers
    nums <- as.numeric(regmatches(mut, gregexpr("\\d+", mut))[[1]])
    if (length(nums) >= 2) {
      if (nums[2] > nums[1]) {
        diff <- nums[2] - nums[1]
        if (diff >= 50) return("large_insertion")
        return("insertion")  # repeat expansion
      } else {
        diff <- nums[1] - nums[2]
        if (diff >= 50) return("large_deletion")
        return("deletion")   # repeat contraction
      }
    }
    return("indel")
  }

  # SNP: single base change A→G, T→C, etc.
  if (grepl("^[ATCG]→[ATCG]$", mut, ignore.case = TRUE)) {
    return("SNP")
  }

  # Complex substitution: multiple bases → multiple bases
  if (grepl("[ATCG]+→[ATCG]+", mut, ignore.case = TRUE)) {
    parts <- strsplit(mut, "→")[[1]]
    if (length(parts) == 2) {
      old_len <- nchar(gsub("[^ATCGatcg]", "", parts[1]))
      new_len <- nchar(gsub("[^ATCGatcg]", "", parts[2]))
      if (old_len == new_len) {
        if (old_len == 1) return("SNP")
        return("substitution")
      } else if (new_len > old_len) {
        return("insertion")
      } else {
        return("deletion")
      }
    }
  }

  # Amplification
  if (grepl("amplification|duplication|×\\d+", mut, ignore.case = TRUE)) {
    return("amplification")
  }

  # Inversion
  if (grepl("inversion|inv", mut, ignore.case = TRUE)) {
    return("inversion")
  }

  # Check annotation for additional clues
  if (!is.null(annotation) && !is.na(annotation)) {
    ann_lower <- tolower(annotation)
    if (grepl("intergenic|promoter|upstream|downstream", ann_lower)) {
      # Intergenic regions - type determined by mutation itself
    }
    if (grepl("pseudogene", ann_lower)) {
      # Still classify by mutation type

    }
  }

  # Missing coverage (from separate sheet)
  if (grepl("missing.coverage|missing_coverage", mut, ignore.case = TRUE)) {
    return("missing_coverage")
  }

  return("unknown")
}

#' Classify Multiple Mutations (Vectorized)
#'
#' @param mutations Character vector of mutation strings
#' @param annotations Optional character vector of annotations
#' @return Character vector of mutation types
#' @export
classify_mutations <- function(mutations, annotations = NULL) {
  if (is.null(annotations)) {
    annotations <- rep(NA, length(mutations))
  }
  mapply(classify_mutation, mutations, annotations, USE.NAMES = FALSE)
}

#' Get Mutation Type Color
#'
#' Returns a color for each mutation type for visualization.
#'
#' @param mutation_type Character string of mutation type
#' @return Hex color code
#' @export
get_mutation_color <- function(mutation_type) {
  colors <- c(
    "SNP"              = "#3498db",  # Blue
    "insertion"        = "#2ecc71",  # Green
    "deletion"         = "#e74c3c",  # Red
    "large_insertion"  = "#27ae60",  # Dark green
    "large_deletion"   = "#c0392b",  # Dark red
    "substitution"     = "#9b59b6",  # Purple
    "amplification"    = "#f39c12",  # Orange
    "inversion"        = "#1abc9c",  # Teal
    "mobile_element"   = "#e67e22",  # Dark orange
    "missing_coverage" = "#95a5a6",  # Gray
    "indel"            = "#f1c40f",  # Yellow
    "unknown"          = "#bdc3c7"   # Light gray
  )

  ifelse(mutation_type %in% names(colors),
         colors[mutation_type],
         colors["unknown"])
}

#' Add Mutation Type Column to Data Frame
#'
#' Adds a 'mutation_type' column to a data frame containing mutation data.
#'
#' @param df Data frame with a 'mutation' column
#' @param mutation_col Name of the column containing mutation strings
#' @param annotation_col Optional name of annotation column
#' @return Data frame with added 'mutation_type' column
#' @export
add_mutation_types <- function(df, mutation_col = "mutation", annotation_col = "annotation") {

  if (!mutation_col %in% names(df)) {
    # Try to find mutation column
    mut_cols <- grep("mutation|mut\\.", names(df), value = TRUE, ignore.case = TRUE)
    if (length(mut_cols) > 0) {
      mutation_col <- mut_cols[1]
    } else {
      warning("No mutation column found. Returning original data frame.")
      return(df)
    }
  }

  annotations <- if (annotation_col %in% names(df)) df[[annotation_col]] else NULL

  df$mutation_type <- classify_mutations(df[[mutation_col]], annotations)
  df$mutation_color <- sapply(df$mutation_type, get_mutation_color)

  return(df)
}
