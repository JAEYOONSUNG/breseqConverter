#' Extract and Process Named Tables from HTML Files
#'
#' This function scans a specified directory for HTML files containing "index" and "html" in their names,
#' processes each file, and extracts named tables. Each processed table is saved to the global environment
#' with a name derived from the file name (excluding the .html extension).
#'
#' @param directory Character string specifying the directory to scan for files. Defaults to the working directory.
#' @param file_pattern Character string for matching file names. Defaults to "index.*\\.html".
#' @return None. The processed tables are saved to the global environment.
#' @import rvest
#' @import dplyr
#' @export
#' @examples
#' # Process files in the current directory
#' extract_named_tables()
extract_named_tables <- function(directory = ".", file_pattern = "index.*\\.html") {
  # Load required libraries
  library(rvest)
  library(dplyr)

  # Find files matching the pattern
  html_files <- list.files(directory, pattern = file_pattern, full.names = TRUE)

  # Check if any files match
  if (length(html_files) == 0) {
    stop("No matching files found in the specified directory.")
  }

  # Function to process a single file
  process_file <- function(file_path) {
    # Read the HTML file
    html_data <- rvest::read_html(file_path)
    # Extract all tables
    tables <- html_data %>% rvest::html_table(fill = TRUE)
    # Remove the first table if unnecessary
    tables <- tables[-1]
    # Assign new names
    names(tables) <- c("Mutation predictions", "Unassigned missing coverage evidence", "Unassigned new junction evidence")

    # Process each table
    tables <- lapply(names(tables), function(name) {
      df <- tables[[name]]
      colnames(df) <- as.character(df[1, ])  # Set first row as headers
      df <- df[-1, ]  # Remove the header row

      # Special handling for specific table types
      if (name == "Unassigned missing coverage evidence") {
        colnames(df)[1:3] <- c("Forward aligned reads", "Reverse aligned reads", "Read Coverage Depth")
      }
      if (name == "Unassigned new junction evidence") {
        colnames(df)[1:2] <- c("aligned reads1", "aligned reads2")
      }

      return(df)
    })
    names(tables) <- c("Mutation predictions", "Unassigned missing coverage evidence", "Unassigned new junction evidence")
    return(tables)
  }

  # Process all matched files and save to the environment
  for (file_path in html_files) {
    # Extract the base file name without the .html extension
    base_name <- tools::file_path_sans_ext(basename(file_path))

    # Process the file and assign the result to the environment
    processed_tables <- process_file(file_path)
    assign(base_name, processed_tables, envir = .GlobalEnv)  # Save to global environment

    cat("Processed tables saved to environment with name:", base_name, "\n")
  }

  cat("All files processed and saved to the environment.\n")
}
