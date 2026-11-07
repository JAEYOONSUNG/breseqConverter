#' Compare Named Tables from Two Lists
#'
#' This function compares named tables from two lists (e.g., experimental and control datasets).
#' Comparison keys include the replicon and coordinate. Junction keys include
#' both oriented endpoints; input junction rows must occur in endpoint pairs.
#'
#' @param list_A A named list of data frames, representing the first set of tables (e.g., control data).
#' @param list_B A named list of data frames, representing the second set of tables (e.g., experimental data).
#' @param name_A,name_B Optional sample suffixes used in output column names.
#' @return A named list of data frames containing the comparison results for each table.
#'         Each data frame includes columns from both input tables and a `status` column indicating:
#'         - "New in B": Rows present in `list_B` but missing in `list_A`.
#'         - "Missing in B": Rows present in `list_A` but missing in `list_B`.
#'         - "Unchanged or Modified": Rows that exist in both but may have differences in other columns.
#' @import dplyr
#' @export
#' @examples
#' # Example usage:
#' comparison_results <- compare_named_tables(list_A, list_B)
#' for (name in names(comparison_results)) {
#'   print(comparison_results[[name]])
#' }
compare_named_tables <- function(list_A, list_B, name_A = NULL, name_B = NULL) {

  # Ensure both lists have the same names for comparison
  if (!all(names(list_A) == names(list_B))) {
    stop("Lists A and B must have the same names for comparison.")
  }

  # Use provided names or fall back to deparse
  suffix_A <- if (!is.null(name_A)) name_A else deparse(substitute(list_A))
  suffix_B <- if (!is.null(name_B)) name_B else deparse(substitute(list_B))
  
  # Initialize a list to store comparison results for each named pair
  comparison_results <- list()
  
  # Iterate through each named data frame in list_A and list_B
  for (name in names(list_A)) {
    cat("Comparing:", name, "\n")
    
    # Retrieve the corresponding tables from list_A and list_B
    df_A <- list_A[[name]]
    df_B <- list_B[[name]]
    
    # Clean column names by removing whitespace and making them unique
    colnames(df_A) <- make.unique(.breseq_text(colnames(df_A)))
    colnames(df_B) <- make.unique(.breseq_text(colnames(df_B)))
    
    # Determine the column to use for comparison based on the table name
    comparison_column <- if (name == "Unassigned missing coverage evidence") "start" else "position"
    
    # Ensure the comparison column exists in both data frames
    if (!(comparison_column %in% colnames(df_A)) || !(comparison_column %in% colnames(df_B))) {
      stop(paste("Table", name, "must contain a", comparison_column, "column in both A and B"))
    }
    
    # Keep the legacy coordinate column while joining by replicon-aware keys.
    # Junction identity includes both oriented endpoints, never one coordinate.
    kind <- .breseq_table_kind(name)
    sequence_values <- unique(unlist(lapply(list(df_A, df_B), function(d) {
      at <- intersect(c("seq id", "seq_id", "seqid"), names(d))
      if (length(at)) .breseq_seq(d[[at[1]]]) else character()
    })))
    sequence_values <- sequence_values[!is.na(sequence_values) & nzchar(sequence_values)]
    has_seq <- vapply(list(df_A, df_B), function(d)
      !nrow(d) || any(c("seq id", "seq_id", "seqid") %in% names(d)), logical(1))
    if (length(sequence_values) > 1L && !all(has_seq))
      stop("A table without seq id cannot be compared to multiple replicons.")
    default_seq <- if (length(sequence_values) == 1L) sequence_values else ""
    key_A <- .breseq_source_keys(df_A, kind, default_seq)
    key_B <- .breseq_source_keys(df_B, kind, default_seq)
    if (anyDuplicated(key_A) || anyDuplicated(key_B))
      stop("Duplicate comparison keys within a sample; refusing a many-to-many join.")
    display_keys <- c(key_A, key_B)
    display_values <- c(df_A[[comparison_column]], df_B[[comparison_column]])
    df_A[[comparison_column]] <- key_A
    df_B[[comparison_column]] <- key_B
    comparison <- dplyr::full_join(df_A, df_B, by = comparison_column,
      suffix = c(paste0(".", suffix_A), paste0(".", suffix_B)))
    keys <- comparison[[comparison_column]]
    comparison$status <- ifelse(!keys %in% key_A, paste("New in", suffix_B),
      ifelse(!keys %in% key_B, paste("Missing in", suffix_B), "Unchanged or Modified"))
    comparison[[comparison_column]] <- display_values[match(keys, display_keys)]
    attr(comparison, "breseq_keys") <- keys
    attr(comparison, "breseq_sample_order") <- c(sub("/output/index$", "", suffix_A),
                                                  sub("/output/index$", "", suffix_B))
    
    # Store the comparison result in the output list with the table name
    comparison_results[[name]] <- comparison
  }
  
  return(comparison_results)
}
