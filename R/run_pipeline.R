#' Run Complete breseqConverter Pipeline
#'
#' Executes the full analysis pipeline: extract mutations from breseq HTML,
#' compare across evolution stages, classify mutation types, optionally
#' create annotated GenBank files for SnapGene, and generate visualizations.
#'
#' @param input_dir Directory containing breseq output folders (each with index.html)
#' @param output_dir Directory for output files
#' @param reference_file Path to reference genome (.dna, .gbk, .fasta).
#'   Genome lengths are automatically extracted from this file.
#'   Also used for SnapGene export if export_snapgene = TRUE.
#' @param genome_lengths Optional: manually specify genome lengths as named vector.
#'   If NULL and reference_file is provided, lengths are auto-detected.
#' @param create_plots Create visualization plots (default: TRUE)
#' @param export_snapgene Export annotated GenBank for SnapGene (default: TRUE if reference_file provided)
#' @param color_by How to color in plots/GenBank: "type" or "stage"
#' @return List containing all results: comparisons, plots, and output file paths
#' @export
#' @examples
#' \dontrun{
#' # Simplest usage - just provide reference file
#' results <- run_breseq_pipeline(
#'   input_dir = "breseq_results/",
#'   output_dir = "analysis_output/",
#'   reference_file = "reference.gbk"  # Genome lengths auto-detected!
#' )
#'
#' # Multi-plasmid genome - all lengths auto-detected from GenBank
#' results <- run_breseq_pipeline(
#'   input_dir = "breseq_results/",
#'   reference_file = "multi_replicon.gbk"  # Contains chromosome + plasmids
#' )
#' }
run_breseq_pipeline <- function(input_dir = ".",
                                output_dir = "breseq_analysis",
                                reference_file = NULL,
                                genome_lengths = NULL,
                                create_plots = TRUE,
                                export_snapgene = NULL,
                                color_by = c("type", "stage")) {

  color_by <- match.arg(color_by)

  # Use current directory if "."
  if (input_dir == ".") {
    input_dir <- getwd()
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("\n", strrep("=", 60))
  message("breseqConverter Pipeline")
  message(strrep("=", 60))
  message("Working directory: ", getwd())

  # Auto-detect GenBank reference file if not provided
  if (is.null(reference_file)) {
    gb_files <- list.files(input_dir, pattern = "\\.(gb|gbk|gbff)$",
                           full.names = TRUE, ignore.case = TRUE)
    # Also check current directory if different
    if (length(gb_files) == 0 && input_dir != getwd()) {
      gb_files <- list.files(".", pattern = "\\.(gb|gbk|gbff)$",
                             full.names = TRUE, ignore.case = TRUE)
    }

    if (length(gb_files) > 0) {
      reference_file <- gb_files[1]
      message("Auto-detected GenBank: ", basename(reference_file))
      if (length(gb_files) > 1) {
        message("  (Found ", length(gb_files), " files, using first one)")
      }
    }
  }

  # Set export_snapgene default
  if (is.null(export_snapgene)) {
    export_snapgene <- !is.null(reference_file)
  }

  # Auto-detect genome lengths from reference file
  replicon_info <- NULL
  if (!is.null(reference_file) && file.exists(reference_file)) {
    message("\n[0/7] Parsing reference genome...")
    tryCatch({
      replicon_info <- parse_reference_info(reference_file)
      if (is.null(genome_lengths)) {
        genome_lengths <- setNames(replicon_info$length, replicon_info$seq_id)
      }
    }, error = function(e) {
      warning("  Could not parse reference file: ", e$message)
    })
  }

  results <- list()

  # Step 1: Extract tables from HTML files
  message("\n[1/7] Extracting mutation tables from breseq HTML files...")

  html_files <- list.files(input_dir, pattern = "index.*\\.html",
                           full.names = TRUE, recursive = TRUE)

  if (length(html_files) == 0) {
    stop("No breseq HTML files found in: ", input_dir)
  }

  message("  Found ", length(html_files), " HTML files")

  extracted_tables <- list()
  for (html_file in html_files) {
    # Get sample name - try filename first, then directory
    filename <- tools::file_path_sans_ext(basename(html_file))

    # Try to extract sample name from filename (e.g., "index_45SJY3" -> "45SJY3")
    sample_match <- regmatches(filename, regexpr("(45)?SJY\\d+", filename, ignore.case = TRUE))
    if (length(sample_match) > 0 && nchar(sample_match) > 0) {
      sample_name <- sample_match
    } else {
      # Fall back to directory name or filename
      dir_name <- basename(dirname(html_file))
      if (dir_name %in% c("output", ".", "extdata", "inst")) {
        sample_name <- filename
      } else {
        sample_name <- dir_name
      }
    }

    message("  Processing: ", sample_name)

    tryCatch({
      tables <- extract_single_html(html_file)
      extracted_tables[[sample_name]] <- tables
    }, error = function(e) {
      warning("  Failed to process ", sample_name, ": ", e$message)
    })
  }

  results$extracted_tables <- extracted_tables

  # Sort samples by natural order (SJY3 -> SJY7 -> SJY8 -> SJY10 -> SJY16)
  sample_names <- names(extracted_tables)
  if (requireNamespace("gtools", quietly = TRUE)) {
    sample_names <- sample_names[gtools::mixedorder(sample_names)]
  } else {
    sample_names <- sort(sample_names)
  }

  # Reorder extracted_tables
  extracted_tables <- extracted_tables[sample_names]
  results$extracted_tables <- extracted_tables
  results$sample_order <- sample_names

  # Show detected order and ask for confirmation
  message("\n[2/7] Detected evolution order:")
  message("  ", paste(sample_names, collapse = " -> "))

  # Interactive confirmation (only in interactive mode)
  if (interactive()) {
    message("\n  Is this the correct order of evolution stages?")
    message("  (Mutations will be tracked cumulatively in this order)")
    response <- readline("  Press ENTER to confirm, or type the correct order (comma-separated): ")

    if (nchar(trimws(response)) > 0) {
      # User provided custom order
      new_order <- trimws(strsplit(response, ",")[[1]])
      new_order <- trimws(new_order)

      # Validate
      if (!all(new_order %in% sample_names)) {
        missing <- setdiff(new_order, sample_names)
        stop("Unknown samples: ", paste(missing, collapse = ", "),
             "\nAvailable: ", paste(sample_names, collapse = ", "))
      }

      sample_names <- new_order
      extracted_tables <- extracted_tables[sample_names]
      results$extracted_tables <- extracted_tables
      results$sample_order <- sample_names
      message("  Using custom order: ", paste(sample_names, collapse = " -> "))
    } else {
      message("  Confirmed!")
    }
  }

  # Step 2b: Compare tables sequentially
  message("\n[3/7] Comparing mutation tables across samples...")

  if (length(extracted_tables) >= 2) {
    comparisons <- compare_multiple_tables(extracted_tables)
    results$comparisons <- comparisons
    message("  Generated ", length(comparisons), " comparisons")
  } else {
    message("  Skipping comparison (need at least 2 samples)")
    comparisons <- NULL
  }

  # Step 3: Merge and classify mutations
  message("\n[4/7] Merging results and classifying mutation types...")

  merged_results <- merge_comparison_results(comparisons)

  # Add mutation type classification
  if ("Mutation predictions" %in% names(merged_results)) {
    merged_results$`Mutation predictions` <- add_mutation_types(
      merged_results$`Mutation predictions`
    )
    message("  Classified ", nrow(merged_results$`Mutation predictions`), " mutations")
  }

  # Handle plasmid/circular genome positions (multi-replicon aware)
  if (!is.null(genome_lengths) && length(genome_lengths) > 0) {
    message("  Normalizing positions for ", length(genome_lengths), " replicon(s)...")
    for (name in names(merged_results)) {
      df <- merged_results[[name]]
      pos_col <- if ("position" %in% names(df)) "position" else "start"
      seq_id_col <- find_seq_id_column(df)

      # Use seq_id for multi-replicon normalization
      if (!is.null(seq_id_col) && length(genome_lengths) > 1) {
        merged_results[[name]] <- identify_position_duplicates(
          df,
          position_col = pos_col,
          genome_length = NULL,  # Will use seq_id grouping
          seq_id_col = seq_id_col
        )
        message("    ", name, ": using seq_id column '", seq_id_col, "' for multi-replicon")
      } else {
        # Single replicon or no seq_id
        merged_results[[name]] <- identify_position_duplicates(
          df,
          position_col = pos_col,
          genome_length = genome_lengths[1]
        )
      }
    }
  }

  results$merged_results <- merged_results

  # Step 4: Save to Excel
  message("\n[5/7] Saving results to Excel...")

  excel_file <- file.path(output_dir, "merged_mutation_results.xlsx")
  save_merged_results_to_excel(merged_results, excel_file)
  results$excel_file <- excel_file
  message("  Saved: ", excel_file)

  # Step 5: Create plots
  if (create_plots) {
    message("\n[6/7] Creating visualizations...")

    plot_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir)

    if ("Mutation predictions" %in% names(merged_results)) {
      mut_df <- merged_results$`Mutation predictions`

      # Get total genome length for plotting (sum of all replicons)
      total_genome_length <- if (!is.null(genome_lengths)) sum(genome_lengths, na.rm = TRUE) else NULL

      # Try circlize-based Circos plot first (preferred)
      circos_available <- requireNamespace("circlize", quietly = TRUE) &&
                          requireNamespace("ComplexHeatmap", quietly = TRUE)

      if (circos_available && !is.null(genome_lengths) && length(genome_lengths) > 0) {
        message("  Creating Circos-style circular map (circlize)...")

        # For multi-replicon, create separate circos plots
        if (length(genome_lengths) > 1) {
          seq_id_col <- find_seq_id_column(mut_df)
          for (rep_name in names(genome_lengths)) {
            tryCatch({
              if (!is.null(seq_id_col) && seq_id_col %in% names(mut_df)) {
                rep_df <- mut_df[mut_df[[seq_id_col]] == rep_name, ]
              } else {
                rep_df <- mut_df
              }
              if (nrow(rep_df) > 0) {
                fname <- paste0("circos_", gsub("[^a-zA-Z0-9]", "_", rep_name), ".pdf")
                plot_circos_genome(
                  mutation_data = rep_df,
                  reference_file = reference_file,
                  output_file = file.path(plot_dir, fname)
                )
                results$plots[[paste0("circos_", rep_name)]] <- file.path(plot_dir, fname)
              }
            }, error = function(e) warning("  Circos plot for ", rep_name, " failed: ", e$message))
          }
        } else {
          tryCatch({
            plot_circos_genome(
              mutation_data = mut_df,
              reference_file = reference_file,
              output_file = file.path(plot_dir, "circos_genome_map.pdf")
            )
            results$plots$circos <- file.path(plot_dir, "circos_genome_map.pdf")
          }, error = function(e) warning("  Circos plot failed: ", e$message))
        }
      }

      # Also create simple ggplot plots as backup/alternative
      if (requireNamespace("ggplot2", quietly = TRUE)) {
        # Linear plot
        tryCatch({
          p1 <- plot_mutations(mut_df, genome_length = total_genome_length,
                              plot_type = "linear", color_by = color_by,
                              title = "Mutation Landscape (Linear)")
          ggplot2::ggsave(file.path(plot_dir, "mutations_linear.png"), p1,
                          width = 12, height = 4, dpi = 150)
          results$plots$linear <- p1
          message("  Saved: mutations_linear.png")
        }, error = function(e) warning("  Linear plot failed: ", e$message))

        # Timeline plot
        tryCatch({
          p3 <- plot_mutations(mut_df, plot_type = "timeline", color_by = color_by,
                              title = "Mutation Accumulation Timeline")
          ggplot2::ggsave(file.path(plot_dir, "mutations_timeline.png"), p3,
                          width = 10, height = 6, dpi = 150)
          results$plots$timeline <- p3
          message("  Saved: mutations_timeline.png")
        }, error = function(e) warning("  Timeline plot failed: ", e$message))

        # Summary bar chart
        tryCatch({
          p4 <- plot_mutation_summary(mut_df)
          ggplot2::ggsave(file.path(plot_dir, "mutation_types.png"), p4,
                          width = 8, height = 5, dpi = 150)
          results$plots$summary <- p4
          message("  Saved: mutation_types.png")
        }, error = function(e) warning("  Summary plot failed: ", e$message))
      }
    }
  } else {
    message("\n[6/7] Skipping visualizations (create_plots = FALSE)")
  }

  # Step 6: Export to SnapGene
  if (export_snapgene && !is.null(reference_file)) {
    message("\n[7/7] Creating annotated GenBank for SnapGene...")

    if (!file.exists(reference_file)) {
      warning("  Reference file not found: ", reference_file)
    } else {
      tryCatch({
        genbank_file <- file.path(output_dir, "annotated_mutations.gbk")
        create_annotated_genbank(
          reference_file = reference_file,
          mutation_data = excel_file,
          output_file = genbank_file,
          color_by = color_by
        )
        results$genbank_file <- genbank_file
        message("  Saved: ", genbank_file)
      }, error = function(e) {
        warning("  GenBank export failed: ", e$message)
      })
    }
  } else {
    message("\n[7/7] Skipping SnapGene export")
  }

  message("\n", strrep("=", 60))
  message("Pipeline complete!")
  message("Output directory: ", normalizePath(output_dir))
  message(strrep("=", 60), "\n")

  return(invisible(results))
}

#' Extract Tables from Single HTML File
#'
#' Internal function to extract tables from a single breseq HTML file.
#'
#' @param html_file Path to HTML file
#' @return Named list of data frames
#' @keywords internal
extract_single_html <- function(html_file) {

  if (!requireNamespace("rvest", quietly = TRUE)) {
    stop("Package 'rvest' required. Install with: install.packages('rvest')")
  }

  html_data <- rvest::read_html(html_file)
  all_tables <- html_data |> rvest::html_table(fill = TRUE)

  if (length(all_tables) == 0) {
    warning("No tables found in ", basename(html_file))
    return(list())
  }

  # Remove first table (usually summary/header info)
  if (length(all_tables) > 1) {
    all_tables <- all_tables[-1]
  }

  # Standard breseq table names
  expected_names <- c("Mutation predictions",
                     "Unassigned missing coverage evidence",
                     "Unassigned new junction evidence")

  # Assign names based on available tables
  n_tables <- min(length(all_tables), length(expected_names))
  table_names <- expected_names[1:n_tables]

  # Process each table
  processed_tables <- list()

  for (i in seq_len(n_tables)) {
    name <- table_names[i]
    df <- all_tables[[i]]

    if (is.null(df) || nrow(df) == 0) {
      processed_tables[[name]] <- data.frame()
      next
    }

    # First row is usually the header
    colnames(df) <- as.character(df[1, ])
    df <- df[-1, , drop = FALSE]

    # Clean column names
    colnames(df) <- trimws(colnames(df))

    # Special column handling for specific table types
    if (name == "Unassigned missing coverage evidence" && ncol(df) >= 3) {
      colnames(df)[1:3] <- c("Forward aligned reads", "Reverse aligned reads", "Read Coverage Depth")
    }
    if (name == "Unassigned new junction evidence" && ncol(df) >= 2) {
      colnames(df)[1:2] <- c("aligned reads1", "aligned reads2")
    }

    processed_tables[[name]] <- df
  }

  return(processed_tables)
}

#' Merge Comparison Results
#'
#' Combines multiple pairwise comparisons into unified tables.
#'
#' @param comparisons List of comparison results from compare_multiple_tables
#' @return Named list of merged data frames
#' @keywords internal
merge_comparison_results <- function(comparisons) {

  if (is.null(comparisons) || length(comparisons) == 0) {
    return(list())
  }

  # Get all table names
  all_tables <- unique(unlist(lapply(comparisons, names)))

  merged <- list()

  for (tbl_name in all_tables) {
    # Collect this table from all comparisons
    tbl_list <- lapply(comparisons, function(comp) {
      if (tbl_name %in% names(comp)) comp[[tbl_name]] else NULL
    })
    tbl_list <- Filter(Negate(is.null), tbl_list)

    if (length(tbl_list) > 0) {
      # Get position column
      pos_col <- if ("position" %in% names(tbl_list[[1]])) "position" else "start"

      # Merge by position - use empty suffix to avoid duplication
      merged_tbl <- tbl_list[[1]]
      for (i in seq_along(tbl_list)[-1]) {
        # Get columns that already exist in merged_tbl (excluding position)
        existing_cols <- setdiff(names(merged_tbl), pos_col)
        new_cols <- setdiff(names(tbl_list[[i]]), pos_col)

        # Only add columns that don't already exist
        cols_to_add <- setdiff(new_cols, existing_cols)
        if (length(cols_to_add) > 0) {
          add_df <- tbl_list[[i]][, c(pos_col, cols_to_add), drop = FALSE]
          merged_tbl <- merge(merged_tbl, add_df, by = pos_col, all = TRUE)
        }
      }

      # Remove status columns and duplicates
      merged_tbl <- merged_tbl[, !grepl("^status", names(merged_tbl))]

      merged[[tbl_name]] <- merged_tbl
    }
  }

  return(merged)
}

#' Save Merged Results to Excel
#'
#' Saves merged comparison results to a multi-sheet Excel file.
#'
#' @param merged_results Named list of data frames
#' @param output_file Path to output Excel file
#' @keywords internal
save_merged_results_to_excel <- function(merged_results, output_file) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' required")
  }

  wb <- openxlsx::createWorkbook()

  sheet_names <- c("Mutation predictions" = "Mutations",
                  "Unassigned missing coverage evidence" = "Missing",
                  "Unassigned new junction evidence" = "New junction")

  for (name in names(merged_results)) {
    sheet_name <- if (name %in% names(sheet_names)) sheet_names[name] else name
    sheet_name <- substr(sheet_name, 1, 31)  # Excel limit

    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeData(wb, sheet_name, merged_results[[name]])

    # Add conditional formatting for mutation types if present
    if ("mutation_type" %in% names(merged_results[[name]])) {
      type_col <- which(names(merged_results[[name]]) == "mutation_type")
      n_rows <- nrow(merged_results[[name]])

      openxlsx::conditionalFormatting(wb, sheet_name,
        cols = type_col, rows = 2:(n_rows + 1),
        rule = "SNP",
        style = openxlsx::createStyle(bgFill = "#3498db")
      )
      openxlsx::conditionalFormatting(wb, sheet_name,
        cols = type_col, rows = 2:(n_rows + 1),
        rule = "deletion",
        style = openxlsx::createStyle(bgFill = "#e74c3c")
      )
      openxlsx::conditionalFormatting(wb, sheet_name,
        cols = type_col, rows = 2:(n_rows + 1),
        rule = "insertion",
        style = openxlsx::createStyle(bgFill = "#2ecc71")
      )
    }
  }

  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
}

#' Quick Analysis Wrapper
#'
#' Simplified wrapper for common analysis scenarios.
#'
#' @param ... Arguments passed to run_breseq_pipeline
#' @return Analysis results
#' @export
analyze_breseq <- function(...) {
  run_breseq_pipeline(...)
}
