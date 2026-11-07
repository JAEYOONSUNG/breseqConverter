#' Save Merged Results to Separate Excel Files
#'
#' This function saves each element of a list of data frames as a separate Excel file.
#' The file names are based on the names of the list elements.
#'
#' @param results A named list of data frames (e.g., the output of `compare_and_merge_multiple_tables`).
#' @param output_dir A character string specifying the directory where the Excel files will be saved.
#'                   Defaults to the current working directory.
#' @param highlight_new Shade every reported call by stage and bold new or changed calls. Default FALSE.
#' @param clean_headers Use readable headers with hidden metadata for canonical restoration.
#' @param sample_labels Named sample-ID to display-label map, used only with clean_headers=TRUE.
#' @param sample_order Optional exact chronological sample order.
#' @return None. The function saves the results as Excel files in the specified directory.
#' @import openxlsx
#' @export
#' @examples
#' # Example usage:
#' save_results_to_excel(merged_results, output_dir = "output")
save_results_to_excel <- function(results, output_dir = ".",
                                   highlight_new = FALSE, clean_headers = FALSE,
                                   sample_labels = NULL, sample_order = NULL) {
  # Load required library
  library(openxlsx)
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each element in the results list
  for (name in names(results)) {
    # Get the data frame
    df <- results[[name]]
    
    # Define the output file path
    file_name <- paste0(name, ".xlsx")
    file_path <- file.path(output_dir, file_name)
    
    # Create a new workbook and write the data frame to it
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
    mutation_names <- names(results)[vapply(names(results), function(n) .breseq_table_kind(n) == "mutation", logical(1))]
    mutation_data <- if (length(mutation_names)) results[[mutation_names[1]]] else NULL
    .breseq_write_presentation(wb, "Sheet1", df, name,
      highlight_new, clean_headers, sample_labels, sample_order, mutation_data)
    
    # Save the workbook
    saveWorkbook(wb, file_path, overwrite = TRUE)
    cat("Saved:", file_path, "\n")
  }
  
  cat("All results have been saved to:", output_dir, "\n")
}
