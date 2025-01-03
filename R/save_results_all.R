#' Save Merged Results to Separate Excel Files
#'
#' This function saves each element of a list of data frames as a separate Excel file.
#' The file names are based on the names of the list elements.
#'
#' @param results A named list of data frames (e.g., the output of `compare_and_merge_multiple_tables`).
#' @param output_dir A character string specifying the directory where the Excel files will be saved.
#'                   Defaults to the current working directory.
#' @return None. The function saves the results as Excel files in the specified directory.
#' @import openxlsx
#' @export
#' @examples
#' # Example usage:
#' save_results_to_excel(merged_results, output_dir = "output")
save_results_to_excel <- function(results, output_dir = ".") {
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
    writeData(wb, "Sheet1", df)
    
    # Save the workbook
    saveWorkbook(wb, file_path, overwrite = TRUE)
    cat("Saved:", file_path, "\n")
  }
  
  cat("All results have been saved to:", output_dir, "\n")
}
