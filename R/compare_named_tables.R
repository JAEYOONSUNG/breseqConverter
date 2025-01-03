#' Compare Named Tables from Two Lists
#'
#' This function compares named tables from two lists (e.g., experimental and control datasets).
#' It identifies differences based on a specified column ("position" or "start") depending on the table name.
#'
#' @param list_A A named list of data frames, representing the first set of tables (e.g., control data).
#' @param list_B A named list of data frames, representing the second set of tables (e.g., experimental data).
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
compare_named_tables <- function(list_A, list_B) {
  
  # Ensure both lists have the same names for comparison
  if (!all(names(list_A) == names(list_B))) {
    stop("Lists A and B must have the same names for comparison.")
  }
  
  # Capture the names of list_A and list_B to use as dynamic suffixes
  suffix_A <- deparse(substitute(list_A))
  suffix_B <- deparse(substitute(list_B))
  
  # Initialize a list to store comparison results for each named pair
  comparison_results <- list()
  
  # Iterate through each named data frame in list_A and list_B
  for (name in names(list_A)) {
    cat("Comparing:", name, "\n")
    
    # Retrieve the corresponding tables from list_A and list_B
    df_A <- list_A[[name]]
    df_B <- list_B[[name]]
    
    # Clean column names by removing whitespace and making them unique
    colnames(df_A) <- make.unique(trimws(colnames(df_A)))
    colnames(df_B) <- make.unique(trimws(colnames(df_B)))
    
    # Determine the column to use for comparison based on the table name
    comparison_column <- if (name == "Unassigned missing coverage evidence") "start" else "position"
    
    # Ensure the comparison column exists in both data frames
    if (!(comparison_column %in% colnames(df_A)) || !(comparison_column %in% colnames(df_B))) {
      stop(paste("Table", name, "must contain a", comparison_column, "column in both A and B"))
    }
    
    # Perform the join and compare by the chosen column with dynamic suffixes
    comparison <- full_join(df_A, df_B, by = comparison_column, suffix = c(paste0(".", suffix_A), paste0(".", suffix_B))) %>%
      mutate(
        status = case_when(
          rowSums(is.na(select(., ends_with(paste0(".", suffix_A))))) == length(select(., ends_with(paste0(".", suffix_A)))) ~ paste("New in", suffix_B),
          rowSums(is.na(select(., ends_with(paste0(".", suffix_B))))) == length(select(., ends_with(paste0(".", suffix_B)))) ~ paste("Missing in", suffix_B),
          TRUE ~ "Unchanged or Modified"
        )
      )
    
    # Store the comparison result in the output list with the table name
    comparison_results[[name]] <- comparison
  }
  
  return(comparison_results)
}