#' Run Full breseqConverter Analysis Pipeline
#'
#' Complete end-to-end workflow from breseq HTML files to Circos visualization.
#' This function combines all steps: extraction, comparison, processing, and visualization.
#'
#' @param input_dir Directory containing breseq HTML files (index.html)
#' @param output_dir Directory for all output files
#' @param reference_file Path to reference genome (.gbk, .gbff, .fasta, or .dna)
#' @param excel_file Optional: pre-existing Excel file with mutation data (skips extraction)
#' @param gff_file Optional: GFF3 file for gene annotations (required for genome map)
#' @param eggnog_file Optional: eggNOG-mapper output file for COG annotations (.xlsx or .tsv)
#' @param add_genome_map Add genome map with COG annotations (NULL=ask, TRUE/FALSE=skip prompt)
#' @param create_circos Create Circos-style circular genome map (default: TRUE)
#' @param create_genbank Create annotated GenBank file for SnapGene (default: TRUE)
#' @return List containing all results and file paths
#' @export
#' @examples
#' \dontrun{
#' # Full analysis from breseq HTML files
#' results <- run_full_analysis(
#'   input_dir = "breseq_outputs/",
#'   output_dir = "analysis/",
#'   reference_file = "reference.gbk"
#' )
#'
#' # Or from existing Excel file
#' results <- run_full_analysis(
#'   excel_file = "mutation_results.xlsx",
#'   output_dir = "analysis/",
#'   reference_file = "reference.gbk"
#' )
#' }
run_full_analysis <- function(input_dir = ".",
                              output_dir = "breseq_analysis",
                              reference_file = NULL,
                              excel_file = NULL,
                              gff_file = NULL,
                              eggnog_file = NULL,
                              add_genome_map = NULL,
                              create_circos = TRUE,
                              create_genbank = TRUE) {

  # Use current directory as default
  if (input_dir == ".") {
    input_dir <- getwd()
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  message("\n", strrep("=", 70))
  message("breseqConverter Full Analysis Pipeline")
  message(strrep("=", 70))
  message("Working directory: ", getwd())

  # Auto-detect GenBank reference file if not provided
  if (is.null(reference_file)) {
    gb_files <- list.files(input_dir, pattern = "\\.(gb|gbk|gbff)$",
                           full.names = TRUE, ignore.case = TRUE)
    if (length(gb_files) == 0 && input_dir != getwd()) {
      gb_files <- list.files(".", pattern = "\\.(gb|gbk|gbff)$",
                             full.names = TRUE, ignore.case = TRUE)
    }
    if (length(gb_files) > 0) {
      reference_file <- gb_files[1]
      message("Auto-detected GenBank: ", basename(reference_file))
    }
  }

  # =============================================================
  # GENOME MAP (GFF + eggNOG) - Interactive check
  # =============================================================

  # Always auto-detect GFF and eggNOG files (regardless of add_genome_map setting)
  message("\nChecking for genome annotation files...")

  # Auto-detect GFF file
  if (is.null(gff_file)) {
    gff_files <- list.files(input_dir, pattern = "\\.(gff|gff3)$",
                            full.names = TRUE, ignore.case = TRUE)
    if (length(gff_files) == 0 && input_dir != getwd()) {
      gff_files <- list.files(".", pattern = "\\.(gff|gff3)$",
                              full.names = TRUE, ignore.case = TRUE)
    }
    if (length(gff_files) > 0) {
      gff_file <- gff_files[1]
      message("  [OK] Auto-detected GFF: ", basename(gff_file))
    }
  } else if (file.exists(gff_file)) {
    message("  [OK] GFF file: ", basename(gff_file))
  }

  # Auto-detect eggNOG file (emapper.annotations.xlsx)
  if (is.null(eggnog_file)) {
    eggnog_files <- list.files(input_dir, pattern = "emapper\\.annotations\\.xlsx$",
                               full.names = TRUE, ignore.case = TRUE)
    if (length(eggnog_files) == 0 && input_dir != getwd()) {
      eggnog_files <- list.files(".", pattern = "emapper\\.annotations\\.xlsx$",
                                 full.names = TRUE, ignore.case = TRUE)
    }
    # Fallback: also check for .tsv format
    if (length(eggnog_files) == 0) {
      eggnog_files <- list.files(input_dir, pattern = "emapper\\.annotations(\\.tsv)?$",
                                 full.names = TRUE, ignore.case = TRUE)
      if (length(eggnog_files) == 0 && input_dir != getwd()) {
        eggnog_files <- list.files(".", pattern = "emapper\\.annotations(\\.tsv)?$",
                                   full.names = TRUE, ignore.case = TRUE)
      }
    }
    if (length(eggnog_files) > 0) {
      eggnog_file <- eggnog_files[1]
      message("  [OK] Auto-detected eggNOG: ", basename(eggnog_file))
    }
  } else if (file.exists(eggnog_file)) {
    message("  [OK] eggNOG file: ", basename(eggnog_file))
  }

  # Check status of annotation files
  gff_ok <- !is.null(gff_file) && file.exists(gff_file)
  eggnog_ok <- !is.null(eggnog_file) && file.exists(eggnog_file)

  if (!gff_ok && !eggnog_ok) {
    message("  [INFO] No annotation files found - COG/RNA tracks will be skipped")
  } else if (!gff_ok) {
    message("  [WARNING] GFF file not found - COG/RNA tracks require GFF")
  } else if (!eggnog_ok) {
    message("  [WARNING] eggNOG file not found - COG colors will use defaults")
  } else {
    message("  [OK] All annotation files found!")
  }

  # =============================================================
  # INTERACTIVE PROMPT: Add genome map? (COG, RNA, GC tracks)
  # Only ask if add_genome_map is NULL and annotation files exist
  # =============================================================
  if (is.null(add_genome_map) && gff_ok) {
    message("\n")
    message(strrep("-", 50))
    message("Genome annotation files detected!")
    message("  - GFF: ", basename(gff_file))
    if (eggnog_ok) message("  - eggNOG: ", basename(eggnog_file))
    message(strrep("-", 50))
    message("\nAdd genome map tracks? (COG 5'->3', COG 3'->5', RNA, GC skew, GC ratio)")
    message("  [Y] Yes - Include genome annotation tracks (recommended)")
    message("  [N] No  - Only show breseq mutation tracks")
    message("")

    # Read user input
    user_input <- readline(prompt = "Enter choice [Y/n]: ")
    user_input <- trimws(toupper(user_input))

    if (user_input == "" || user_input == "Y" || user_input == "YES") {
      add_genome_map <- TRUE
      message("  -> Including genome map tracks\n")
    } else {
      add_genome_map <- FALSE
      message("  -> Skipping genome map tracks (breseq mutations only)\n")
    }
  } else if (is.null(add_genome_map)) {
    # No GFF file found, default to FALSE
    add_genome_map <- FALSE
  }

  # Auto-detect Excel file if not provided and no HTML extraction needed
  if (is.null(excel_file)) {
    xlsx_files <- list.files(input_dir, pattern = "mutant_results\\.xlsx$",
                             full.names = TRUE, ignore.case = TRUE)
    if (length(xlsx_files) == 0) {
      xlsx_files <- list.files(".", pattern = "mutant_results\\.xlsx$",
                               full.names = TRUE, ignore.case = TRUE)
    }
    if (length(xlsx_files) > 0) {
      excel_file <- xlsx_files[1]
      message("Auto-detected Excel: ", basename(excel_file))
    }
  }

  results <- list()

  # Get genome info from reference
  genome_info <- NULL
  if (!is.null(reference_file) && file.exists(reference_file)) {
    message("\n[1/5] Parsing reference genome...")
    tryCatch({
      genome_info <- parse_reference_info(reference_file)
      results$genome_info <- genome_info
    }, error = function(e) {
      warning("Could not parse reference: ", e$message)
    })
  } else {
    message("\n[1/5] No reference file provided, skipping genome parsing...")
  }

  # Step 2: Get mutation data (from HTML or Excel)
  mutation_data <- NULL

  if (!is.null(excel_file)) {
    message("\n[2/5] Loading Excel file...")
    message("    Path: ", excel_file)
    message("    Exists: ", file.exists(excel_file))

    if (file.exists(excel_file)) {
      mutation_data <- load_mutation_excel(excel_file)
      results$excel_file <- excel_file
    } else {
      # Try to find in current directory
      base_name <- basename(excel_file)
      if (file.exists(base_name)) {
        message("    Found in current directory: ", base_name)
        mutation_data <- load_mutation_excel(base_name)
        results$excel_file <- base_name
      } else {
        stop("Excel file not found: ", excel_file,
             "\n  Current working directory: ", getwd())
      }
    }
  } else {
    # No Excel file - extract from HTML files
    message("\n[2/5] Extracting mutations from breseq HTML files...")
    message("    Input directory: ", input_dir)

    # Check for HTML files
    html_files <- list.files(input_dir, pattern = "index.*\\.html$",
                             recursive = TRUE, full.names = TRUE)
    if (length(html_files) == 0) {
      stop("No breseq HTML files (index*.html) found in: ", input_dir,
           "\n  And no existing mutant_results.xlsx found.",
           "\n  Please check your working directory or provide excel_file path.")
    }

    message("    Found ", length(html_files), " HTML files")

    # Run extraction pipeline
    pipeline_results <- run_breseq_pipeline(
      input_dir = input_dir,
      output_dir = output_dir,
      reference_file = reference_file,
      create_plots = FALSE,
      export_snapgene = FALSE
    )

    excel_file <- pipeline_results$excel_file
    mutation_data <- load_mutation_excel(excel_file)
    results$excel_file <- excel_file
    results$pipeline_results <- pipeline_results
  }

  if (is.null(mutation_data)) {
    stop("Failed to load mutation data")
  }

  # Step 3: Process mutation data for visualization
  message("\n[3/5] Processing mutation data...")
  processed_data <- process_mutation_data_for_viz(mutation_data)
  results$mutation_list <- processed_data$mutation_list
  results$missing_data <- processed_data$missing_data
  results$samples <- processed_data$samples

  # Step 4: Create Circos visualization
  if (create_circos) {
    message("\n[4/5] Creating Circos visualization...")

    plot_dir <- file.path(output_dir, "plots")
    if (!dir.exists(plot_dir)) dir.create(plot_dir)

    # Initialize annotation data (only load if add_genome_map == TRUE)
    gff_data <- NULL
    eggnog_data <- NULL
    gc_data <- NULL

    if (isTRUE(add_genome_map)) {
      # Parse GFF file if provided
      if (!is.null(gff_file) && file.exists(gff_file)) {
        message("  Parsing GFF file...")
        gff_data <- parse_gff_file(gff_file)
        results$gff_data <- gff_data
      }

      # Parse eggNOG file if provided
      if (!is.null(eggnog_file) && file.exists(eggnog_file)) {
        message("  Parsing eggNOG file...")
        eggnog_data <- parse_eggnog_file(eggnog_file)
        results$eggnog_data <- eggnog_data
      }

      # Calculate GC content if reference file is available
      if (!is.null(reference_file) && file.exists(reference_file)) {
        message("  Calculating GC content...")
        gc_data <- calculate_gc_content(reference_file, window_size = 10000)
        results$gc_data <- gc_data
      }
    } else {
      message("  Skipping genome annotation tracks (add_genome_map = FALSE)")
    }

    tryCatch({
      circos_file <- file.path(plot_dir, "circos_genome_map.pdf")
      create_full_circos_plot(
        mutation_list = processed_data$mutation_list,
        missing_data = processed_data$missing_data,
        samples = processed_data$samples,
        genome_info = genome_info,
        gff_data = gff_data,         # NULL if add_genome_map = FALSE
        eggnog_data = eggnog_data,   # NULL if add_genome_map = FALSE
        gc_data = gc_data,           # NULL if add_genome_map = FALSE
        output_file = circos_file
      )
      results$circos_file <- circos_file
      message("  Saved: ", circos_file)
    }, error = function(e) {
      warning("Circos plot failed: ", e$message)
    })
  } else {
    message("\n[4/5] Skipping Circos visualization")
  }

  # Step 5: Create annotated GenBank
  if (create_genbank && !is.null(reference_file)) {
    message("\n[5/5] Creating annotated GenBank file...")

    tryCatch({
      genbank_file <- file.path(output_dir, "annotated_mutations.gbk")
      create_annotated_genbank(
        reference_file = reference_file,
        mutation_data = excel_file,
        output_file = genbank_file
      )
      results$genbank_file <- genbank_file
    }, error = function(e) {
      warning("GenBank export failed: ", e$message)
    })
  } else {
    message("\n[5/5] Skipping GenBank export")
  }

  message("\n", strrep("=", 70))
  message("Analysis complete!")
  message("Output directory: ", normalizePath(output_dir))
  message(strrep("=", 70), "\n")

  return(invisible(results))
}


#' Load Mutation Data from Excel File
#'
#' Reads mutation Excel file with Mutations, Missing, New junction sheets.
#'
#' @param excel_file Path to Excel file
#' @return List with data frames for each sheet
#' @keywords internal
load_mutation_excel <- function(excel_file) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' required")
  }

  sheet_names <- c("Mutations", "Missing", "New junction")

  data_list <- list()
  for (sheet in sheet_names) {
    tryCatch({
      data_list[[sheet]] <- openxlsx::read.xlsx(excel_file, sheet = sheet)
      message("  Loaded sheet '", sheet, "': ", nrow(data_list[[sheet]]), " rows")
    }, error = function(e) {
      message("  Sheet '", sheet, "' not found or empty")
      data_list[[sheet]] <- NULL
    })
  }

  return(data_list)
}


#' Parse GFF3 File for Gene Annotations
#'
#' Reads GFF3 file and extracts gene/CDS features with strand information.
#'
#' @param gff_file Path to GFF3 file
#' @return Data frame with seqid, start, end, strand, type, gene, product
#' @keywords internal
parse_gff_file <- function(gff_file) {

  if (!file.exists(gff_file)) {
    warning("GFF file not found: ", gff_file)
    return(NULL)
  }

  # Read GFF file (skip comment lines)
  lines <- readLines(gff_file, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]

  if (length(lines) == 0) {
    warning("GFF file is empty or contains only comments")
    return(NULL)
  }

  # Parse GFF columns
  gff_data <- do.call(rbind, lapply(lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 9) {
      # Parse attributes
      attrs <- fields[9]
      get_attr <- function(attr_name) {
        pattern <- paste0(attr_name, "=([^;]+)")
        m <- regmatches(attrs, regexpr(pattern, attrs))
        if (length(m) > 0) {
          return(gsub(paste0(attr_name, "="), "", m))
        }
        return(NA)
      }

      data.frame(
        seqid = fields[1],
        source = fields[2],
        type = fields[3],
        start = as.numeric(fields[4]),
        end = as.numeric(fields[5]),
        score = fields[6],
        strand = fields[7],
        phase = fields[8],
        gene = get_attr("gene"),
        locus_tag = get_attr("locus_tag"),
        product = get_attr("product"),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }))

  if (is.null(gff_data) || nrow(gff_data) == 0) {
    warning("No valid features found in GFF file")
    return(NULL)
  }

  # Filter for relevant features (CDS, RNA types, gene)
  relevant_types <- c("CDS", "gene", "mRNA", "tRNA", "rRNA", "tmRNA", "ncRNA", "region")
  gff_filtered <- gff_data[gff_data$type %in% relevant_types, ]

  message("    Parsed ", nrow(gff_filtered), " features from GFF")
  message("    Feature types: ", paste(names(table(gff_filtered$type)), collapse = ", "))
  if (any(gff_data$type == "CDS")) {
    message("    CDS with locus_tag: ", sum(!is.na(gff_filtered$locus_tag) & gff_filtered$type == "CDS"))
  }

  return(gff_filtered)
}


#' Parse eggNOG-mapper Output for COG Annotations
#'
#' Reads eggNOG-mapper annotations file and extracts COG categories.
#' Supports both .xlsx and .tsv/.annotations formats.
#' Query column is cleaned to match GFF locus_tag format (GenomeDrawer method).
#'
#' @param eggnog_file Path to eggNOG-mapper output file (.xlsx or .tsv)
#' @return Data frame with query, COG_category, description, preferred_name, COG_color
#' @keywords internal
parse_eggnog_file <- function(eggnog_file) {

  if (!file.exists(eggnog_file)) {
    warning("eggNOG file not found: ", eggnog_file)
    return(NULL)
  }

  eggnog_data <- NULL

  # Check file extension and parse accordingly
  if (grepl("\\.xlsx$", eggnog_file, ignore.case = TRUE)) {
    # Parse Excel format
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' required for reading Excel files")
    }
    message("    Reading eggNOG Excel file...")

    tryCatch({
      # GenomeDrawer uses startRow = 3 for eggNOG Excel files
      eggnog_raw <- openxlsx::read.xlsx(eggnog_file, sheet = 1, startRow = 3, colNames = TRUE)

      # Remove last 3 rows if necessary (GenomeDrawer standard)
      if (nrow(eggnog_raw) > 3) {
        eggnog_raw <- eggnog_raw[1:(nrow(eggnog_raw) - 3), ]
      }

      # Standardize column names (handle variations)
      col_mapping <- c(
        "query" = "query",
        "#query" = "query",
        "COG_category" = "COG_category",
        "COG.category" = "COG_category",
        "cog_category" = "COG_category",
        "Description" = "description",
        "description" = "description",
        "Preferred_name" = "preferred_name",
        "preferred_name" = "preferred_name"
      )

      # Rename columns if needed
      for (old_name in names(col_mapping)) {
        if (old_name %in% names(eggnog_raw)) {
          names(eggnog_raw)[names(eggnog_raw) == old_name] <- col_mapping[old_name]
        }
      }

      # Select relevant columns
      required_cols <- c("query", "COG_category")
      if (!all(required_cols %in% names(eggnog_raw))) {
        # Try startRow = 1 as fallback
        message("    Trying alternative Excel parsing (startRow = 1)...")
        eggnog_raw <- openxlsx::read.xlsx(eggnog_file, sheet = 1, startRow = 1, colNames = TRUE)
        for (old_name in names(col_mapping)) {
          if (old_name %in% names(eggnog_raw)) {
            names(eggnog_raw)[names(eggnog_raw) == old_name] <- col_mapping[old_name]
          }
        }
        if (!all(required_cols %in% names(eggnog_raw))) {
          warning("Required columns not found in eggNOG Excel file")
          return(NULL)
        }
      }

      eggnog_data <- eggnog_raw[, intersect(names(eggnog_raw),
                                             c("query", "COG_category", "description", "preferred_name"))]

    }, error = function(e) {
      warning("Failed to read eggNOG Excel file: ", e$message)
      return(NULL)
    })

  } else {
    # Parse TSV/text format
    lines <- readLines(eggnog_file, warn = FALSE)

    # Find header line (starts with #query or just query)
    header_idx <- which(grepl("^#?query", lines, ignore.case = TRUE))[1]
    if (is.na(header_idx)) {
      data_lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]
      if (length(data_lines) == 0) {
        warning("Cannot parse eggNOG file format")
        return(NULL)
      }
      header_idx <- 0
    }

    # Parse header
    if (header_idx > 0) {
      header <- gsub("^#", "", lines[header_idx])
      col_names <- strsplit(header, "\t")[[1]]
      data_lines <- lines[(header_idx + 1):length(lines)]
      data_lines <- data_lines[!grepl("^#", data_lines) & nchar(data_lines) > 0]
    } else {
      col_names <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                     "max_annot_lvl", "COG_category", "Description",
                     "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway",
                     "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
                     "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")
      data_lines <- lines[!grepl("^#", lines) & nchar(lines) > 0]
    }

    if (length(data_lines) == 0) {
      warning("No data found in eggNOG file")
      return(NULL)
    }

    # Parse data
    eggnog_data <- do.call(rbind, lapply(data_lines, function(line) {
      fields <- strsplit(line, "\t")[[1]]
      if (length(fields) >= 9) {
        data.frame(
          query = fields[1],
          COG_category = if (length(fields) >= 7) fields[7] else NA,
          description = if (length(fields) >= 8) fields[8] else NA,
          preferred_name = if (length(fields) >= 9) fields[9] else NA,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
  }

  if (is.null(eggnog_data) || nrow(eggnog_data) == 0) {
    warning("No valid annotations found in eggNOG file")
    return(NULL)
  }

  # ============================================================
  # GenomeDrawer query cleanup logic (EggNOG_annotations 함수와 동일)
  # Clean up the 'query' column to match GFF locus_tag format
  # ============================================================
  message("    Cleaning query column (GenomeDrawer method)...")

  clean_query <- function(q) {
    if (is.na(q) || q == "") return(q)
    # Remove gnl| prefix
    if (grepl("gnl\\|", q)) {
      q <- gsub("gnl\\|", "", q)
    }
    # Replace | with :
    if (grepl("\\|", q)) {
      q <- gsub("\\|", ":", q)
    }
    # Remove lcl:NC/AC/BG/NT/NW/NZ prefix pattern
    if (grepl("lcl:(AC|NC|BG|NT|NW|NZ)_([a-zA-Z]+)?[0-9]+\\.[0-9]+_prot_", q)) {
      q <- gsub("lcl:(AC|NC|BG|NT|NW|NZ)_([a-zA-Z]+)?[0-9]+\\.[0-9]+_prot_", "", q)
    }
    # Remove _nnnn suffix (protein sequence number)
    if (grepl("_[0-9]{1,4}$", q)) {
      q <- gsub("_[0-9]{1,4}$", "", q)
    }
    # Remove extdb: prefix
    if (grepl("^extdb:", q)) {
      q <- gsub("^extdb:", "", q)
    }
    return(q)
  }

  eggnog_data$query_original <- eggnog_data$query
  eggnog_data$query <- sapply(eggnog_data$query, clean_query)

  # Filter out numeric-only queries
  eggnog_data <- eggnog_data[!grepl("^[0-9]+$", eggnog_data$query), ]

  # Remove rows with empty COG category
  eggnog_data <- eggnog_data[!is.na(eggnog_data$COG_category) &
                               eggnog_data$COG_category != "" &
                               eggnog_data$COG_category != "-", ]

  message("    Parsed ", nrow(eggnog_data), " annotations with COG categories")

  # ============================================================
  # Add COG colors (GenomeDrawer standard colors)
  # ============================================================
  cog_colors <- get_cog_colors()

  # Create COG_category_for_plot (first letter only)
  eggnog_data$COG_category_for_plot <- substr(eggnog_data$COG_category, 1, 1)

  # Assign COG colors using GenomeDrawer color scheme
  eggnog_data$COG_color <- sapply(eggnog_data$COG_category_for_plot, function(cog) {
    if (is.na(cog) || cog == "" || cog == "-") return("#CCCCCC")
    if (cog %in% names(cog_colors)) {
      return(cog_colors[cog])
    }
    return("#CCCCCC")  # Default gray for unknown COG
  })

  # Debug: Show sample query transformations
  if (nrow(eggnog_data) > 0) {
    sample_idx <- utils::head(which(!is.na(eggnog_data$query_original) &
                               eggnog_data$query_original != eggnog_data$query), 3)
    if (length(sample_idx) > 0) {
      message("    Query cleanup examples:")
      for (i in sample_idx) {
        message("      '", eggnog_data$query_original[i], "' -> '", eggnog_data$query[i], "'")
      }
    }
    message("    Sample cleaned queries: ", paste(utils::head(eggnog_data$query, 5), collapse = ", "))
  }

  return(eggnog_data)
}


#' Get COG Category Colors (GenomeDrawer Standard)
#'
#' Returns GenomeDrawer's standard COG category color scheme.
#' These colors match exactly with GenomeDrawer/EggNOG_annotations() function.
#'
#' @return Named vector of colors for each COG category
#' @export
#' @examples
#' colors <- get_cog_colors()
#' colors["J"]  # "#ff0000" - Translation
get_cog_colors <- function() {
  # GenomeDrawer 표준 COG 색상 (EggNOG_annotations 함수와 동일)
  c(
    "J" = "#ff0000",  # Translation, ribosomal structure and biogenesis
    "A" = "#c2af58",  # RNA processing and modification
    "K" = "#ff9900",  # Transcription
    "L" = "#ffff00",  # Replication, recombination and repair
    "B" = "#ffc600",  # Chromatin structure and dynamics
    "D" = "#99ff00",  # Cell cycle control, cell division, chromosome partitioning
    "Y" = "#493126",  # Nuclear structure
    "V" = "#ff008a",  # Defense mechanisms
    "T" = "#0000ff",  # Signal transduction mechanisms
    "M" = "#9ec928",  # Cell wall/membrane/envelope biogenesis
    "N" = "#006633",  # Cell motility
    "Z" = "#660099",  # Cytoskeleton
    "W" = "#336699",  # Extracellular structures
    "U" = "#33cc99",  # Intracellular trafficking, secretion, and vesicular transport
    "O" = "#00ffff",  # Posttranslational modification, protein turnover, chaperones
    "C" = "#9900ff",  # Energy production and conversion
    "G" = "#805642",  # Carbohydrate transport and metabolism
    "E" = "#ff00ff",  # Amino acid transport and metabolism
    "F" = "#99334d",  # Nucleotide transport and metabolism
    "H" = "#727dcc",  # Coenzyme transport and metabolism
    "I" = "#5c5a1b",  # Lipid transport and metabolism
    "P" = "#0099ff",  # Inorganic ion transport and metabolism
    "Q" = "#ffcc99",  # Secondary metabolites biosynthesis, transport and catabolism
    "R" = "#ff9999",  # General function prediction only
    "S" = "#d6aadf"   # Function unknown
  )
}


#' Get COG Category Legends (GenomeDrawer Standard)
#'
#' Returns GenomeDrawer's standard COG category descriptions.
#' These match exactly with GenomeDrawer/EggNOG_annotations() function.
#'
#' @return Named vector of legend descriptions for each COG category
#' @export
#' @examples
#' legends <- get_cog_legends()
#' legends["J"]  # "[J] Translation, ribosomal structure and biogenesis"
get_cog_legends <- function() {
  c(
    "J" = "[J] Translation, ribosomal structure and biogenesis",
    "A" = "[A] RNA processing and modification",
    "K" = "[K] Transcription",
    "L" = "[L] Replication, recombination and repair",
    "B" = "[B] Chromatin structure and dynamics",
    "D" = "[D] Cell cycle control, cell division, chromosome partitioning",
    "Y" = "[Y] Nuclear structure",
    "V" = "[V] Defense mechanisms",
    "T" = "[T] Signal transduction mechanisms",
    "M" = "[M] Cell wall/membrane/envelope biogenesis",
    "N" = "[N] Cell motility",
    "Z" = "[Z] Cytoskeleton",
    "W" = "[W] Extracellular structures",
    "U" = "[U] Intracellular trafficking, secretion, and vesicular transport",
    "O" = "[O] Posttranslational modification, protein turnover, chaperones",
    "C" = "[C] Energy production and conversion",
    "G" = "[G] Carbohydrate transport and metabolism",
    "E" = "[E] Amino acid transport and metabolism",
    "F" = "[F] Nucleotide transport and metabolism",
    "H" = "[H] Coenzyme transport and metabolism",
    "I" = "[I] Lipid transport and metabolism",
    "P" = "[P] Inorganic ion transport and metabolism",
    "Q" = "[Q] Secondary metabolites biosynthesis, transport and catabolism",
    "R" = "[R] General function prediction only",
    "S" = "[S] Function unknown"
  )
}


#' Calculate GC Content Statistics from Reference Genome
#'
#' Calculates GC skew and GC ratio in sliding windows across the genome.
#' Uses GenomeDrawer formula: GC skew = (C - G) / (C + G), GC ratio = (G + C) / total
#' Note: GenomeDrawer uses (C - G) not (G - C) for consistency with origin of replication analysis.
#'
#' @param reference_file Path to reference genome file (GenBank, FASTA, or .dna)
#' @param window_size Size of sliding window (default: 10000 bp, GenomeDrawer standard)
#' @return Data frame with start, end, gc_skew, gc_ratio, gc_skew_minus_average, gc_ratio_minus_average, and colors
#' @keywords internal
calculate_gc_content <- function(reference_file, window_size = 10000) {

  if (!file.exists(reference_file)) {
    warning("Reference file not found: ", reference_file)
    return(NULL)
  }

  # Read sequence based on file type
  sequence <- NULL
  seqid <- "Genome"

  if (grepl("\\.(gb|gbk|gbff)$", reference_file, ignore.case = TRUE)) {
    # Parse GenBank format
    lines <- readLines(reference_file, warn = FALSE)

    # Get LOCUS name
    locus_line <- lines[grepl("^LOCUS", lines)][1]
    if (!is.na(locus_line)) {
      seqid <- strsplit(trimws(sub("^LOCUS\\s+", "", locus_line)), "\\s+")[[1]][1]
    }

    # Extract ORIGIN sequence
    origin_start <- which(grepl("^ORIGIN", lines))
    end_line <- which(grepl("^//", lines))

    if (length(origin_start) > 0 && length(end_line) > 0) {
      seq_lines <- lines[(origin_start[1] + 1):(end_line[1] - 1)]
      sequence <- toupper(paste(gsub("[^ATCGatcg]", "", seq_lines), collapse = ""))
    }

  } else if (grepl("\\.(fa|fasta|fna)$", reference_file, ignore.case = TRUE)) {
    # Parse FASTA format
    lines <- readLines(reference_file, warn = FALSE)
    header_idx <- which(grepl("^>", lines))

    if (length(header_idx) > 0) {
      seqid <- sub("^>\\s*", "", lines[header_idx[1]])
      seqid <- strsplit(seqid, "\\s+")[[1]][1]

      if (length(header_idx) == 1) {
        seq_lines <- lines[(header_idx[1] + 1):length(lines)]
      } else {
        seq_lines <- lines[(header_idx[1] + 1):(header_idx[2] - 1)]
      }
      sequence <- toupper(paste(seq_lines, collapse = ""))
    }
  }

  if (is.null(sequence) || nchar(sequence) < window_size) {
    warning("Could not extract sequence or sequence too short")
    return(NULL)
  }

  genome_length <- nchar(sequence)
  message("    Calculating GC content for ", seqid, " (", format(genome_length, big.mark = ","), " bp)")
  message("    Using GenomeDrawer formula: GC skew = (C - G) / (C + G)")

  # Calculate GC in sliding windows
  n_windows <- ceiling(genome_length / window_size)
  gc_data <- data.frame(
    seqid = character(n_windows),
    start = numeric(n_windows),
    end = numeric(n_windows),
    gc_skew = numeric(n_windows),
    gc_ratio = numeric(n_windows),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_windows)) {
    win_start <- (i - 1) * window_size + 1
    win_end <- min(i * window_size, genome_length)
    win_seq <- substr(sequence, win_start, win_end)

    n_G <- nchar(gsub("[^G]", "", win_seq))
    n_C <- nchar(gsub("[^C]", "", win_seq))
    n_total <- nchar(win_seq)

    gc_data$seqid[i] <- seqid
    gc_data$start[i] <- win_start
    gc_data$end[i] <- win_end

    # GC skew = (C - G) / (C + G) - GenomeDrawer formula (nC - nG)/(nC + nG)
    if ((n_G + n_C) > 0) {
      gc_data$gc_skew[i] <- (n_C - n_G) / (n_C + n_G)
    } else {
      gc_data$gc_skew[i] <- 0
    }

    # GC ratio = (G + C) / total
    gc_data$gc_ratio[i] <- (n_G + n_C) / n_total
  }

  # Calculate deviation from average (GenomeDrawer naming: _minus_average)
  avg_skew <- mean(gc_data$gc_skew, na.rm = TRUE)
  avg_ratio <- mean(gc_data$gc_ratio, na.rm = TRUE)

  # GenomeDrawer naming convention: gc_skew_minus_average, gc_ratio_minus_average
  gc_data$gc_skew_minus_average <- gc_data$gc_skew - avg_skew
  gc_data$gc_ratio_minus_average <- gc_data$gc_ratio - avg_ratio

  # Also keep deviation naming for compatibility
  gc_data$gc_skew_deviation <- gc_data$gc_skew_minus_average
  gc_data$gc_ratio_deviation <- gc_data$gc_ratio_minus_average

  # Assign colors (트랙 출력 순서에 맞게 스왑됨)
  # 결과: GC skew 트랙 = + 파란색(#0066CC), - 노란색(#FFD700)
  # 결과: GC ratio 트랙 = + 초록색(#228B22), - 빨간색(#DC143C)
  gc_data$gc_skew_color <- ifelse(gc_data$gc_skew_minus_average >= 0, "#228B22", "#DC143C")
  gc_data$gc_ratio_color <- ifelse(gc_data$gc_ratio_minus_average >= 0, "#0066CC", "#FFD700")

  message("    GC skew range: ", round(min(gc_data$gc_skew), 3), " to ", round(max(gc_data$gc_skew), 3))
  message("    GC ratio range: ", round(min(gc_data$gc_ratio), 3), " to ", round(max(gc_data$gc_ratio), 3))
  message("    Average GC skew: ", round(avg_skew, 4), ", Average GC ratio: ", round(avg_ratio, 4))

  return(gc_data)
}


#' Process Mutation Data for Visualization
#'
#' Reshapes mutation data from wide format to long format,
#' identifies new mutations at each ALE stage, and processes missing regions.
#'
#' @param mutation_data List from load_mutation_excel()
#' @return List with mutation_list, missing_data, and samples
#' @keywords internal
process_mutation_data_for_viz <- function(mutation_data) {

  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' required. Install with: install.packages('reshape2')")
  }
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' required. Install with: install.packages('gtools')")
  }

  mutations_df <- mutation_data[["Mutations"]]
  missing_df <- mutation_data[["Missing"]]

  mutation_list <- list()
  missing_processed <- NULL
  samples <- character()

  # Process Mutations
  if (!is.null(mutations_df) && nrow(mutations_df) > 0) {
    message("  Processing mutations...")

    # Clean position column
    if ("position" %in% names(mutations_df)) {
      mutations_df$position <- as.numeric(gsub(",", "", as.character(mutations_df$position)))
    }

    # Identify columns
    base_cols <- c("position")
    value_cols <- c("evidence", "mutation", "annotation", "gene", "description")
    all_cols <- colnames(mutations_df)
    tag_cols <- setdiff(all_cols, c(base_cols, value_cols))

    # Extract sample names from column names
    get_sample_name <- function(col_name) {
      if (col_name %in% value_cols) return("None")

      # Pattern 1: 45SJY followed by numbers (e.g., 45SJY3, 45SJY16)
      m <- regmatches(col_name, regexpr("45SJY\\d+", col_name, ignore.case = TRUE))
      if (length(m) > 0 && nchar(m) > 0) return(m)

      # Pattern 2: SJY followed by numbers (e.g., SJY3, SJY16)
      m <- regmatches(col_name, regexpr("SJY\\d+", col_name, ignore.case = TRUE))
      if (length(m) > 0 && nchar(m) > 0) return(m)

      # Pattern 3: Extract from dot-separated format (mutation.SAMPLE/output/index)
      parts <- strsplit(col_name, "\\.")[[1]]
      if (length(parts) > 1) {
        tag <- sub("/output.*", "", parts[2])
        tag <- sub("/.*", "", tag)  # Remove anything after first /
        if (nchar(tag) > 0 && tag != col_name) return(tag)
      }

      return("Unknown")
    }

    # Debug: Print column names and detected samples
    message("    Column names: ", paste(utils::head(all_cols, 10), collapse = ", "),
            if(length(all_cols) > 10) "..." else "")
    test_samples <- unique(sapply(tag_cols, get_sample_name))
    message("    Detected samples from columns: ", paste(test_samples, collapse = ", "))

    # Melt to long format
    melted <- reshape2::melt(
      mutations_df,
      id.vars = "position",
      measure.vars = setdiff(all_cols, "position"),
      variable.name = "variable",
      value.name = "value"
    )

    # Extract sample and variable type
    melted$ALE_version <- sapply(as.character(melted$variable), get_sample_name)
    melted$variable_type <- sapply(as.character(melted$variable), function(v) {
      if (grepl("^evidence", v)) return("evidence")
      if (grepl("^mutation", v)) return("mutation")
      if (grepl("^annotation", v)) return("annotation")
      if (grepl("^gene", v)) return("gene")
      if (grepl("^description", v)) return("description")
      return("other")
    })

    # Filter and cast
    melted <- melted[!is.na(melted$value) & melted$variable_type %in% value_cols, ]

    if (nrow(melted) > 0) {
      mutation_casted <- reshape2::dcast(
        melted,
        position + ALE_version ~ variable_type,
        value.var = "value",
        fun.aggregate = function(x) paste(stats::na.omit(x), collapse = ";")
      )

      # Ensure required columns exist
      for (col in c("mutation", "annotation", "gene", "description")) {
        if (!col %in% names(mutation_casted)) mutation_casted[[col]] <- NA
      }

      # Create label
      mutation_casted$label <- paste(
        ifelse(is.na(mutation_casted$mutation), "", mutation_casted$mutation),
        ifelse(is.na(mutation_casted$annotation), "", mutation_casted$annotation),
        sep = " | "
      )
      mutation_casted$label <- trimws(gsub("^\\s*\\|\\s*|\\s*\\|\\s*$", "", mutation_casted$label))

      # Get sorted sample names
      samples <- unique(mutation_casted$ALE_version)
      samples <- samples[samples != "None" & samples != "Unknown"]
      samples <- samples[gtools::mixedorder(samples)]

      # Create sample colors
      n_samples <- length(samples)
      if (n_samples > 0) {
        sample_colors <- setNames(
          RColorBrewer::brewer.pal(max(n_samples, 3), "Set2")[1:n_samples],
          samples
        )

        # Track new mutations per sample
        previous_mutations <- data.frame(position = numeric(), mutation = character())

        for (sample in samples) {
          sample_data <- mutation_casted[mutation_casted$ALE_version == sample, ]

          if (nrow(sample_data) > 0) {
            # Mark new mutations
            sample_data$new <- !paste(sample_data$position, sample_data$mutation, sep = ":") %in%
              paste(previous_mutations$position, previous_mutations$mutation, sep = ":")

            sample_data$version_color <- sample_colors[sample]

            mutation_list[[sample]] <- sample_data
            previous_mutations <- rbind(
              previous_mutations,
              sample_data[, c("position", "mutation")]
            )
          }
        }
      }

      message("    Found ", length(samples), " samples: ", paste(samples, collapse = " -> "))

      # Debug: Show mutation counts per sample
      for (s in samples) {
        if (s %in% names(mutation_list)) {
          total <- nrow(mutation_list[[s]])
          new_count <- sum(mutation_list[[s]]$new, na.rm = TRUE)
          message("      ", s, ": ", total, " total, ", new_count, " new")
        }
      }
    }
  }

  # Process Missing Coverage
  if (!is.null(missing_df) && nrow(missing_df) > 0) {
    message("  Processing missing coverage regions...")
    missing_processed <- process_missing_regions(missing_df, samples)
  }

  return(list(
    mutation_list = mutation_list,
    missing_data = missing_processed,
    samples = samples
  ))
}


#' Process Missing Coverage Regions
#'
#' @param missing_df Missing coverage data frame
#' @param samples Sample names
#' @return Processed missing data with certain/uncertain regions
#' @keywords internal
process_missing_regions <- function(missing_df, samples) {

  if (!requireNamespace("reshape2", quietly = TRUE)) return(NULL)

  # Similar processing to mutations
  base_cols <- c("start")
  value_cols <- c("seq id", "end", "size", "gene", "description",
                  "Forward aligned reads", "Reverse aligned reads", "Read Coverage Depth")
  all_cols <- colnames(missing_df)

  # Get sample names
  get_sample_name <- function(col_name) {
    if (col_name %in% value_cols) return("None")

    # Pattern 1: 45SJY followed by numbers
    m <- regmatches(col_name, regexpr("45SJY\\d+", col_name, ignore.case = TRUE))
    if (length(m) > 0 && nchar(m) > 0) return(m)

    # Pattern 2: SJY followed by numbers
    m <- regmatches(col_name, regexpr("SJY\\d+", col_name, ignore.case = TRUE))
    if (length(m) > 0 && nchar(m) > 0) return(m)

    # Fallback
    parts <- strsplit(col_name, "\\.")[[1]]
    if (length(parts) > 1) {
      tag <- sub("/output.*", "", parts[2])
      tag <- sub("/.*", "", tag)
      if (nchar(tag) > 0) return(tag)
    }
    return("Unknown")
  }

  # Melt
  melted <- reshape2::melt(
    missing_df,
    id.vars = "start",
    measure.vars = setdiff(all_cols, "start"),
    variable.name = "variable",
    value.name = "value"
  )

  melted$ALE_version <- sapply(as.character(melted$variable), get_sample_name)
  melted$variable_type <- sapply(as.character(melted$variable), function(v) {
    for (val_col in value_cols) {
      if (grepl(paste0("^", gsub("([.|()\\^{}+$*?])", "\\\\\\1", val_col)), v)) {
        return(val_col)
      }
    }
    return("other")
  })

  melted <- melted[!is.na(melted$value), ]

  if (nrow(melted) == 0) return(NULL)

  # Cast
  missing_casted <- reshape2::dcast(
    melted,
    start + ALE_version ~ variable_type,
    value.var = "value",
    fun.aggregate = function(x) paste(na.omit(x), collapse = ";")
  )

  # Parse coordinate ranges for certain/uncertain regions
  parse_coord_range <- function(coord_str) {
    if (is.na(coord_str) || coord_str == "") return(list(min = NA, max = NA))
    coord_str <- as.character(coord_str)
    if (!grepl("[–-]", coord_str)) {
      val <- as.numeric(gsub(",", "", coord_str))
      return(list(min = val, max = val))
    }
    parts <- strsplit(coord_str, "[–-]")[[1]]
    if (length(parts) == 2) {
      return(list(
        min = as.numeric(gsub(",", "", parts[1])),
        max = as.numeric(gsub(",", "", parts[2]))
      ))
    }
    return(list(min = NA, max = NA))
  }

  # Add parsed coordinates
  start_parsed <- lapply(missing_casted$start, parse_coord_range)
  missing_casted$start_min <- sapply(start_parsed, `[[`, "min")
  missing_casted$start_max <- sapply(start_parsed, `[[`, "max")

  if ("end" %in% names(missing_casted)) {
    end_parsed <- lapply(missing_casted$end, parse_coord_range)
    missing_casted$end_min <- sapply(end_parsed, `[[`, "min")
    missing_casted$end_max <- sapply(end_parsed, `[[`, "max")
  } else {
    missing_casted$end_min <- missing_casted$start_min
    missing_casted$end_max <- missing_casted$start_max
  }

  # Determine region type and calculate certain/uncertain regions
  missing_casted$region_type <- ifelse(
    !is.na(missing_casted$end_min) & !is.na(missing_casted$end_max) &
      missing_casted$end_min > missing_casted$end_max,
    "MC_inversion", "MC_repeat"
  )

  # Certain region
  missing_casted$certain_start <- missing_casted$start_max
  missing_casted$certain_end <- ifelse(
    missing_casted$region_type == "MC_inversion",
    missing_casted$end_max,
    missing_casted$end_min
  )

  # Uncertain regions
  missing_casted$uncertain_left_start <- ifelse(
    missing_casted$start_min != missing_casted$start_max,
    missing_casted$start_min, NA
  )
  missing_casted$uncertain_left_end <- ifelse(
    missing_casted$start_min != missing_casted$start_max,
    missing_casted$start_max, NA
  )

  missing_casted$uncertain_right_start <- ifelse(
    missing_casted$end_min != missing_casted$end_max,
    ifelse(missing_casted$region_type == "MC_inversion",
           missing_casted$end_max, missing_casted$end_min), NA
  )
  missing_casted$uncertain_right_end <- ifelse(
    missing_casted$end_min != missing_casted$end_max,
    ifelse(missing_casted$region_type == "MC_inversion",
           missing_casted$end_min, missing_casted$end_max), NA
  )

  # Add seqid if available
  if ("seq id" %in% names(missing_casted)) {
    missing_casted$seqid <- missing_casted$`seq id`
  } else {
    missing_casted$seqid <- NA
  }

  message("    Found ", nrow(missing_casted), " missing regions")

  return(missing_casted)
}


#' Create Full Circos Plot
#'
#' Creates comprehensive Circos-style genome map with mutation and missing tracks.
#' Optional annotation tracks from GFF and eggNOG files can be added to the outer edge.
#'
#' @param mutation_list List of mutation data frames by sample
#' @param missing_data Processed missing coverage data
#' @param samples Sample names in order
#' @param genome_info Genome info from parse_reference_info()
#' @param gff_data Parsed GFF data for gene annotation tracks (optional)
#' @param eggnog_data Parsed eggNOG data for COG annotation tracks (optional)
#' @param output_file Output PDF file path
#' @return Path to created file
#' @keywords internal
create_full_circos_plot <- function(mutation_list,
                                    missing_data = NULL,
                                    samples,
                                    genome_info = NULL,
                                    gff_data = NULL,
                                    eggnog_data = NULL,
                                    gc_data = NULL,
                                    output_file = "circos_map.pdf") {

  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' required")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' required")
  }

  library(circlize)

  # Get genome length (note: circlize functions used without namespace inside library() context)
  if (!is.null(genome_info) && nrow(genome_info) > 0) {
    genome_length <- genome_info$length[1]
    seqid <- genome_info$seq_id[1]
  } else {
    # Estimate from data
    all_positions <- unlist(lapply(mutation_list, function(x) x$position))
    genome_length <- max(all_positions, na.rm = TRUE) * 1.1
    seqid <- "Genome"
  }

  message("    Genome: ", seqid, " (", format(genome_length, big.mark = ","), " bp)")
  message("    Samples: ", paste(samples, collapse = " → "))

  # Create sample colors
  n_samples <- length(samples)
  sample_colors <- setNames(
    RColorBrewer::brewer.pal(max(n_samples, 3), "Set2")[1:n_samples],
    samples
  )

  # Mutation type colors (사용자 지정 팔레트)
  type_colors <- c(
    "INS" = "#FF5C8D",   # Bright Pink (삽입)
    "DEL" = "#03045E",   # Dark Navy (결실)
    "SNP" = "#139487",   # Teal (SNP)
    "REP" = "#E5C87B",   # Golden Yellow (치환)
    "UNKNOWN" = "gray"   # Gray
  )

  # Missing region colors
  missing_colors <- c("MC_inversion" = "#FF6B6B", "MC_repeat" = "#4ECDC4")

  # Open PDF
  grDevices::pdf(file = output_file, width = 14, height = 14)

  # Initialize circos with proper canvas settings
  circlize::circos.clear()

  # Calculate total tracks needed for canvas sizing
  n_annotation_tracks <- 0
  if (!is.null(gff_data) && nrow(gff_data) > 0) n_annotation_tracks <- n_annotation_tracks + 2  # 5'->3' and 3'->5'
  if (!is.null(eggnog_data) && nrow(eggnog_data) > 0) n_annotation_tracks <- n_annotation_tracks + 2  # COG tracks

  # circos.par settings (GenomeDrawer/breesq_drawer.R 설정값 그대로)
  circlize::circos.par(
    "start.degree" = 90,
    "track.height" = 0.8,
    "gap.degree" = 0,
    "cell.padding" = c(0, 0, 0, 0),
    "canvas.xlim" = c(-1.2, 1.2),
    "canvas.ylim" = c(-1.2, 1.2)
  )

  circlize::circos.initialize(factors = seqid, xlim = c(1, genome_length))

  # Axis settings
  step_size <- genome_length / 12
  step_size <- ceiling(step_size / (10^(nchar(as.character(as.integer(step_size))) - 1))) *
    10^(nchar(as.character(as.integer(step_size))) - 1)
  brk <- seq(0, genome_length, by = step_size)

  scale_factor <- if (genome_length >= 1e6) 1e6 else 1e3
  label_unit <- if (genome_length >= 1e6) "Mb" else "kb"

  # =============================================================
  # GenomeDrawer Track Order (Outside → Inside):
  # 1. Mutation labels (outside)
  # 2. Mutation tracks (outer)
  # 3. MC track
  # 4. Axis track
  # 5. Genome base circle
  # 6. COG 3'→5' track (reverse strand, inner)
  # 7. COG 5'→3' track (forward strand, inner)
  # 8. Center text
  # =============================================================

  # Prepare COG data for later use (inner tracks)
  gff_forward <- NULL
  gff_reverse <- NULL
  rna_data <- NULL

  if (!is.null(gff_data) && nrow(gff_data) > 0) {
    message("    Processing GFF data: ", nrow(gff_data), " features")
    message("    Feature types: ", paste(unique(gff_data$type), collapse = ", "))

    # Extract CDS features for COG tracks
    gff_cds <- gff_data[gff_data$type == "CDS", ]
    message("    CDS features: ", nrow(gff_cds))

    # Extract RNA features (tRNA, rRNA) for RNA track
    rna_types <- c("tRNA", "rRNA", "tmRNA", "ncRNA")
    rna_data <- gff_data[gff_data$type %in% rna_types, ]
    if (nrow(rna_data) > 0) {
      message("    RNA features: ", nrow(rna_data), " (", paste(table(rna_data$type), collapse = ", "), ")")
      # Assign colors by RNA type (GenomeDrawer style)
      rna_data$color <- sapply(rna_data$type, function(t) {
        switch(t,
          "tRNA" = "#E41A1C",   # Red
          "rRNA" = "#377EB8",   # Blue
          "tmRNA" = "#4DAF4A",  # Green
          "ncRNA" = "#984EA3",  # Purple
          "#999999"             # Default gray
        )
      })
    } else {
      message("    No RNA features found in GFF")
    }

    if (nrow(gff_cds) > 0) {
      # If eggNOG data is available, merge COG info
      if (!is.null(eggnog_data) && nrow(eggnog_data) > 0) {
        message("    Merging eggNOG COG annotations (GenomeDrawer method)...")

        # Debug: Show sample data formats
        sample_gff_lt <- utils::head(stats::na.omit(gff_cds$locus_tag), 5)
        sample_eggnog_q <- utils::head(eggnog_data$query, 5)
        message("    [DEBUG] Sample GFF locus_tags: ", paste(sample_gff_lt, collapse = ", "))
        message("    [DEBUG] Sample eggNOG queries (cleaned): ", paste(sample_eggnog_q, collapse = ", "))

        # Get COG colors from get_cog_colors() (GenomeDrawer standard)
        cog_color_map <- get_cog_colors()

        # =================================================================
        # GenomeDrawer matching strategy:
        # 1. Direct match: locus_tag == query (after query cleanup)
        # 2. Clean both sides and try matching
        # =================================================================

        # Create lookup tables from cleaned eggNOG queries
        eggnog_lookup <- setNames(eggnog_data$COG_category, eggnog_data$query)
        eggnog_color_lookup <- setNames(eggnog_data$COG_color, eggnog_data$query)

        # Also create lookups without version suffix for fallback
        eggnog_data$query_no_version <- gsub("\\.[0-9]+$", "", eggnog_data$query)
        eggnog_lookup_no_version <- setNames(eggnog_data$COG_category, eggnog_data$query_no_version)
        eggnog_color_lookup_no_version <- setNames(eggnog_data$COG_color, eggnog_data$query_no_version)

        # Match function with multiple strategies
        match_cog <- function(lt) {
          if (is.na(lt) || lt == "") return(list(category = NA, color = "#CCCCCC"))

          # Strategy 1: Direct match
          if (lt %in% names(eggnog_lookup)) {
            cat <- eggnog_lookup[lt]
            col <- eggnog_color_lookup[lt]
            return(list(category = cat, color = col))
          }

          # Strategy 2: Remove version suffix from locus_tag
          lt_no_version <- gsub("\\.[0-9]+$", "", lt)
          if (lt_no_version %in% names(eggnog_lookup_no_version)) {
            cat <- eggnog_lookup_no_version[lt_no_version]
            col <- eggnog_color_lookup_no_version[lt_no_version]
            return(list(category = cat, color = col))
          }

          # Strategy 3: Partial match (locus_tag contained in query or vice versa)
          partial_match_idx <- which(grepl(lt, eggnog_data$query, fixed = TRUE) |
                                       grepl(lt_no_version, eggnog_data$query, fixed = TRUE))
          if (length(partial_match_idx) > 0) {
            cat <- eggnog_data$COG_category[partial_match_idx[1]]
            col <- eggnog_data$COG_color[partial_match_idx[1]]
            return(list(category = cat, color = col))
          }

          return(list(category = NA, color = "#CCCCCC"))
        }

        # Apply matching to all CDS features
        match_results <- lapply(gff_cds$locus_tag, match_cog)
        gff_cds$COG_category <- sapply(match_results, `[[`, "category")
        gff_cds$COG_color <- sapply(match_results, `[[`, "color")

        # Validate COG_color (use color map for consistency)
        gff_cds$COG_color <- sapply(seq_len(nrow(gff_cds)), function(i) {
          cog <- gff_cds$COG_category[i]
          if (is.na(cog) || cog == "" || cog == "-") return("#CCCCCC")
          first_cog <- substr(cog, 1, 1)
          if (first_cog %in% names(cog_color_map)) {
            return(cog_color_map[first_cog])
          }
          return("#CCCCCC")
        })

        matched_cog <- sum(!is.na(gff_cds$COG_category) & gff_cds$COG_category != "" & gff_cds$COG_category != "-")
        message("    COG matched: ", matched_cog, " / ", nrow(gff_cds), " CDS features (",
                round(matched_cog / nrow(gff_cds) * 100, 1), "%)")

        # Show COG category distribution and color assignment
        if (matched_cog > 0) {
          cog_table <- table(substr(gff_cds$COG_category[!is.na(gff_cds$COG_category) & gff_cds$COG_category != ""], 1, 1))
          message("    COG categories found: ", paste(names(cog_table), collapse = ", "))

          # Debug: Show a few color assignments
          colored_idx <- which(!is.na(gff_cds$COG_category) & gff_cds$COG_category != "")[1:min(3, matched_cog)]
          message("    Sample COG color assignments:")
          for (idx in colored_idx) {
            message("      ", gff_cds$locus_tag[idx], " -> COG=", gff_cds$COG_category[idx],
                    " -> Color=", gff_cds$COG_color[idx])
          }
        }
      } else {
        # No eggNOG - use default colors by strand
        message("    No eggNOG data - using default strand colors")
        gff_cds$COG_color <- ifelse(gff_cds$strand == "+", "#6495ED", "#FF7F50")
      }

      gff_forward <- gff_cds[gff_cds$strand == "+", ]
      gff_reverse <- gff_cds[gff_cds$strand == "-", ]
      message("    Forward strand CDS: ", nrow(gff_forward), ", Reverse strand CDS: ", nrow(gff_reverse))
    }
  }

  # =============================================================
  # MUTATION PROCESSING
  # =============================================================

  # Parse mutation type
  parse_mut_type <- function(mut) {
    if (is.na(mut) || nchar(as.character(mut)) == 0) return("UNKNOWN")
    mut <- as.character(mut)
    # Insertion: +ATCG
    if (grepl("^\\+[ATCG]+$", mut, ignore.case = TRUE)) return("INS")
    # Homopolymer repeat: (T)6→7
    if (grepl("^\\([ATCG]+\\)\\d+→\\d+$", mut)) {
      nums <- as.numeric(regmatches(mut, gregexpr("\\d+", mut))[[1]])
      if (length(nums) >= 2) {
        return(if (nums[2] > nums[1]) "INS" else "DEL")
      }
    }
    # Replacement: 10 bp→ATCG
    if (grepl("\\d+\\s*bp.*→", mut)) return("REP")
    # SNP: A→G
    if (grepl("^[ATCG]→[ATCG]$", mut, ignore.case = TRUE)) return("SNP")
    # Deletion: Δ100 bp
    if (grepl("^[∆Δ]", mut)) return("DEL")
    return("UNKNOWN")
  }

  # Parse mutation end position
  parse_mut_end <- function(mut, pos) {
    if (is.na(mut) || nchar(as.character(mut)) == 0) return(pos)
    mut <- as.character(mut)
    # Insertion
    if (grepl("^\\+[ATCG]+$", mut, ignore.case = TRUE)) {
      return(pos + nchar(gsub("\\+", "", mut)))
    }
    # Homopolymer repeat
    if (grepl("^\\([ATCG]+\\)\\d+→\\d+$", mut)) {
      nums <- as.numeric(regmatches(mut, gregexpr("\\d+", mut))[[1]])
      if (length(nums) >= 2) return(pos + max(nums))
    }
    # Deletion: Δ2,340 bp
    if (grepl("^[∆Δ][0-9,]+", mut)) {
      num_str <- regmatches(mut, regexpr("[0-9,]+", mut))
      num <- as.numeric(gsub(",", "", num_str))
      if (!is.na(num)) return(pos + num - 1)
    }
    return(pos)
  }

  # Combine all mutation data
  all_mutations <- do.call(rbind, lapply(names(mutation_list), function(sample) {
    df <- mutation_list[[sample]]
    df$sample <- sample
    df$type <- sapply(df$mutation, parse_mut_type)
    df$type_color <- ifelse(df$type %in% names(type_colors),
                            type_colors[df$type], type_colors["UNKNOWN"])
    df$start <- as.numeric(gsub(",", "", df$position))
    # Calculate proper end position based on mutation type
    df$end <- mapply(parse_mut_end, df$mutation, df$start)
    # Ensure end >= start and visible on plot (minimum 500bp for visibility)
    df$end <- pmax(df$end, df$start + 500)
    df
  }))

  # Filter new mutations for plotting
  new_mutations <- all_mutations[all_mutations$new == TRUE, ]

  # Add labels for new mutations
  if (nrow(new_mutations) > 0 && any(!is.na(new_mutations$label) & new_mutations$label != "")) {
    label_data <- new_mutations[!is.na(new_mutations$label) & new_mutations$label != "", ]
    if (nrow(label_data) > 0) {
      label_bed <- data.frame(
        seqid = seqid,
        start = label_data$start,
        end = label_data$end,
        label = label_data$label,
        color = label_data$version_color,
        stringsAsFactors = FALSE
      )

      tryCatch({
        circlize::circos.genomicLabels(
          bed = label_bed,
          labels.column = 4,
          side = "outside",
          col = label_bed$color,
          line_col = label_bed$color,
          cex = 0.5,
          connection_height = circlize::mm_h(5),
          line_lwd = 0.5
        )
      }, error = function(e) message("    Warning: Labels failed"))
    }
  }

  # Mutation tracks (cumulative by sample)
  # REVERSE ORDER: outermost track = most accumulated (last sample)
  track_idx <- circlize::get.current.track.index() + 1
  reversed_samples <- rev(samples)  # 45SJY16 -> 45SJY10 -> ... -> 45SJY3
  mutation_track_start <- track_idx

  for (i in seq_along(reversed_samples)) {
    sample <- reversed_samples[i]
    # Get index in original order for cumulative calculation
    orig_idx <- which(samples == sample)
    # Cumulative: all mutations up to this sample in original order
    cumulative_samples <- samples[1:orig_idx]
    cum_data <- all_mutations[all_mutations$sample %in% cumulative_samples, ]

    bg_color <- adjustcolor(sample_colors[sample], alpha.f = 0.2)

    message("    Track ", i, " (", sample, "): ", nrow(cum_data), " cumulative mutations")

    if (nrow(cum_data) > 0) {
      # Pre-compute data for panel.fun closure
      plot_data <- cum_data

      circlize::circos.track(
        factors = seqid,
        ylim = c(0, 1),
        track.index = track_idx,
        bg.col = bg_color,
        bg.border = FALSE,
        track.height = 0.03,
        track.margin = c(0.001, 0.001),
        panel.fun = function(x, y) {
          for (j in seq_len(nrow(plot_data))) {
            # COLOR BY MUTATION TYPE (not version!)
            circlize::circos.rect(
              xleft = plot_data$start[j],
              xright = plot_data$end[j],
              ybottom = 0.1, ytop = 0.9,
              col = plot_data$type_color[j],  # Type color!
              border = NA
            )
          }
        }
      )
      track_idx <- track_idx + 1
    } else {
      # Empty track
      circlize::circos.track(
        factors = seqid,
        ylim = c(0, 1),
        track.index = track_idx,
        bg.col = bg_color,
        bg.border = FALSE,
        track.height = 0.03,
        track.margin = c(0.001, 0.001)
      )
      track_idx <- track_idx + 1
    }
  }

  # =============================================================
  # MISSING COVERAGE TRACK - Single consolidated track
  # =============================================================
  if (!is.null(missing_data) && nrow(missing_data) > 0) {
    # Filter to certain regions only
    certain_data <- missing_data[!is.na(missing_data$certain_start) &
                                   !is.na(missing_data$certain_end), ]

    if (nrow(certain_data) > 0) {
      certain_data$start_num <- as.numeric(gsub(",", "", certain_data$certain_start))
      certain_data$end_num <- as.numeric(gsub(",", "", certain_data$certain_end))

      # Track which missing regions are "new" at each version and assign color
      seen_mc <- list()  # Store mc_id -> first_sample mapping
      certain_data$first_sample <- NA
      certain_data$sample_color <- NA

      if ("ALE_version" %in% names(certain_data)) {
        for (s in samples) {
          s_rows <- which(certain_data$ALE_version == s)
          for (r in s_rows) {
            mc_id <- paste(certain_data$start_num[r], certain_data$end_num[r], sep = ":")
            if (!(mc_id %in% names(seen_mc))) {
              seen_mc[[mc_id]] <- s
              certain_data$first_sample[r] <- s
              certain_data$sample_color[r] <- sample_colors[s]
            }
          }
        }
      }

      # Get unique MC regions (first appearance only)
      unique_mc <- certain_data[!is.na(certain_data$first_sample), ]
      unique_mc <- unique_mc[!duplicated(paste(unique_mc$start_num, unique_mc$end_num)), ]

      message("    Missing coverage: ", nrow(unique_mc), " unique regions")

      # Draw SINGLE MC track with all regions
      if (nrow(unique_mc) > 0) {
        mc_plot <- unique_mc
        circlize::circos.track(
          factors = seqid,
          ylim = c(0, 1),
          track.index = track_idx,
          bg.col = adjustcolor("#F5F5F5", alpha.f = 0.5),
          bg.border = FALSE,
          track.height = 0.025,
          track.margin = c(0.002, 0.002),
          panel.fun = function(x, y) {
            for (j in seq_len(nrow(mc_plot))) {
              # Use sample color for when it first appeared
              col <- mc_plot$sample_color[j]
              if (is.null(col) || is.na(col)) col <- "#808080"
              # Add pattern or border to distinguish MC type
              border_col <- missing_colors[mc_plot$region_type[j]]
              if (is.null(border_col) || is.na(border_col)) border_col <- "black"
              circlize::circos.rect(
                xleft = mc_plot$start_num[j],
                xright = mc_plot$end_num[j],
                ybottom = 0.1, ytop = 0.9,
                col = adjustcolor(col, alpha.f = 0.6),
                border = border_col, lwd = 1
              )
            }
          }
        )
        track_idx <- track_idx + 1
      }
    }
  }

  # Axis track
  circlize::circos.track(
    track.index = circlize::get.current.track.index(),
    panel.fun = function(x, y) {
      circlize::circos.axis(
        h = "top",
        major.at = brk,
        labels = paste(round(brk / scale_factor, 1), label_unit),
        labels.cex = 0.7,
        col = "grey40",
        labels.col = "grey40",
        lwd = 0.7,
        labels.facing = "clockwise"
      )
    },
    bg.border = FALSE
  )

  # Genome base track
  circlize::circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      circlize::circos.text(
        circlize::CELL_META$xlim[1] + diff(circlize::CELL_META$xlim) / 2,
        0.5, seqid,
        cex = 0.4, col = "grey40",
        facing = "bending.inside", niceFacing = TRUE
      )
    },
    bg.col = "grey70",
    bg.border = FALSE,
    track.height = 0.01
  )

  # =============================================================
  # COG ANNOTATION TRACKS (INNER - after genome base, GenomeDrawer style)
  # =============================================================
  track_idx <- circlize::get.current.track.index() + 1

  if (!is.null(gff_reverse) && nrow(gff_reverse) > 0) {
    # COG 3' to 5' direction (reverse strand)
    message("    Adding COG track (3'→5'): ", nrow(gff_reverse), " genes")
    rev_data <- gff_reverse
    circlize::circos.track(
      factors = seqid,
      ylim = c(0, 1),
      track.index = track_idx,
      bg.col = "grey90",
      bg.border = FALSE,
      track.height = 0.03,
      track.margin = c(0.001, 0.001),
      panel.fun = function(x, y) {
        for (j in seq_len(nrow(rev_data))) {
          circlize::circos.rect(
            xleft = rev_data$start[j],
            xright = rev_data$end[j],
            ybottom = 0.1, ytop = 0.9,
            col = rev_data$COG_color[j],
            border = NA
          )
        }
      }
    )
    track_idx <- track_idx + 1
  }

  if (!is.null(gff_forward) && nrow(gff_forward) > 0) {
    # COG 5' to 3' direction (forward strand)
    message("    Adding COG track (5'→3'): ", nrow(gff_forward), " genes")
    fwd_data <- gff_forward
    circlize::circos.track(
      factors = seqid,
      ylim = c(0, 1),
      track.index = track_idx,
      bg.col = "grey90",
      bg.border = FALSE,
      track.height = 0.03,
      track.margin = c(0.001, 0.001),
      panel.fun = function(x, y) {
        for (j in seq_len(nrow(fwd_data))) {
          circlize::circos.rect(
            xleft = fwd_data$start[j],
            xright = fwd_data$end[j],
            ybottom = 0.1, ytop = 0.9,
            col = fwd_data$COG_color[j],
            border = NA
          )
        }
      }
    )
    track_idx <- track_idx + 1
  }

  # =============================================================
  # RNA TRACK (tRNA, rRNA - GenomeDrawer style)
  # =============================================================
  if (!is.null(rna_data) && nrow(rna_data) > 0) {
    message("    Adding RNA track: ", nrow(rna_data), " features")
    rna_plot <- rna_data
    circlize::circos.track(
      factors = seqid,
      ylim = c(0, 1),
      track.index = track_idx,
      bg.col = "grey90",
      bg.border = FALSE,
      track.height = 0.03,
      track.margin = c(0.001, 0.001),
      panel.fun = function(x, y) {
        for (j in seq_len(nrow(rna_plot))) {
          circlize::circos.rect(
            xleft = rna_plot$start[j],
            xright = rna_plot$end[j],
            ybottom = 0.1, ytop = 0.9,
            col = rna_plot$color[j],
            border = NA
          )
        }
      }
    )
    track_idx <- track_idx + 1
  }

  # =============================================================
  # GC SKEW AND GC RATIO TRACKS (GenomeDrawer/breesq_drawer.R 방식: circos.barplot 사용)
  # GenomeDrawer 네이밍: gc_skew_minus_average, gc_ratio_minus_average
  # =============================================================
  if (!is.null(gc_data) && nrow(gc_data) > 0) {
    message("    Adding GC skew track (circos.barplot, GenomeDrawer style)...")

    gc_plot <- gc_data

    # Check for GenomeDrawer column naming (gc_skew_minus_average)
    skew_col <- if ("gc_skew_minus_average" %in% names(gc_plot)) "gc_skew_minus_average" else "gc_skew_deviation"
    ratio_col <- if ("gc_ratio_minus_average" %in% names(gc_plot)) "gc_ratio_minus_average" else "gc_ratio_deviation"

    max_skew <- max(abs(gc_plot[[skew_col]]), na.rm = TRUE)

    # Force evaluation of colors before closure (closure 문제 방지)
    skew_values <- gc_plot[[skew_col]]
    skew_positions <- (gc_plot$start + gc_plot$end) / 2
    skew_colors <- gc_plot$gc_skew_color  # + 파란색 (#0066CC), - 노란색 (#FFD700)

    # GC Skew track using circos.barplot (breesq_drawer.R line 386-396)
    circlize::circos.track(
      factors = seqid,
      ylim = c(-max_skew, max_skew),
      track.index = track_idx,
      bg.col = NA,
      bg.border = FALSE,
      track.height = 0.06,
      panel.fun = function(x, y) {
        circlize::circos.barplot(
          value = skew_values,
          pos = skew_positions,
          col = skew_colors,
          border = NA,
          bar_width = 10000
        )
      }
    )
    track_idx <- track_idx + 1

    message("    Adding GC ratio track (circos.barplot, GenomeDrawer style)...")

    max_ratio <- max(abs(gc_plot[[ratio_col]]), na.rm = TRUE)

    # Force evaluation of colors before closure (closure 문제 방지)
    ratio_values <- gc_plot[[ratio_col]]
    ratio_positions <- (gc_plot$start + gc_plot$end) / 2
    ratio_colors <- gc_plot$gc_ratio_color  # + 초록색 (#228B22), - 빨간색 (#DC143C)

    # GC Ratio track using circos.barplot (breesq_drawer.R style)
    circlize::circos.track(
      factors = seqid,
      ylim = c(-max_ratio, max_ratio),
      track.index = track_idx,
      bg.col = NA,
      bg.border = FALSE,
      track.height = 0.06,
      panel.fun = function(x, y) {
        circlize::circos.barplot(
          value = ratio_values,
          pos = ratio_positions,
          col = ratio_colors,
          border = NA,
          bar_width = 10000
        )
      }
    )
    track_idx <- track_idx + 1
  }

  # Center text
  graphics::text(0, 0, paste(
    seqid,
    paste(format(genome_length, big.mark = ","), "bp"),
    sep = "\n"
  ), cex = 1.2)

  # =============================================================
  # LEGENDS - 분리하여 배치
  # Left-bottom: Mutation, RNA, Sample, MC legends
  # Right-bottom: COG legend (별도 배치)
  # =============================================================

  lgd_list_left <- list()   # 왼쪽 아래 legends
  lgd_cog <- NULL           # COG legend (오른쪽 아래)

  # COG legend (if eggNOG data present) - GenomeDrawer 스타일 -> 오른쪽 아래
  if (!is.null(eggnog_data) && nrow(eggnog_data) > 0) {
    cog_colors <- get_cog_colors()
    used_cogs <- unique(substr(eggnog_data$COG_category, 1, 1))
    used_cogs <- used_cogs[used_cogs %in% names(cog_colors)]
    used_cogs <- used_cogs[!is.na(used_cogs)]
    if (length(used_cogs) > 0) {
      # GenomeDrawer COG category descriptions (EggNOG_annotations 함수와 동일)
      cog_descriptions <- c(
        "J" = "Translation, ribosomal structure and biogenesis",
        "A" = "RNA processing and modification",
        "K" = "Transcription",
        "L" = "Replication, recombination and repair",
        "B" = "Chromatin structure and dynamics",
        "D" = "Cell cycle control, cell division",
        "Y" = "Nuclear structure",
        "V" = "Defense mechanisms",
        "T" = "Signal transduction mechanisms",
        "M" = "Cell wall/membrane/envelope biogenesis",
        "N" = "Cell motility",
        "Z" = "Cytoskeleton",
        "W" = "Extracellular structures",
        "U" = "Intracellular trafficking, secretion",
        "O" = "PTM, protein turnover, chaperones",
        "C" = "Energy production and conversion",
        "G" = "Carbohydrate transport and metabolism",
        "E" = "Amino acid transport and metabolism",
        "F" = "Nucleotide transport and metabolism",
        "H" = "Coenzyme transport and metabolism",
        "I" = "Lipid transport and metabolism",
        "P" = "Inorganic ion transport and metabolism",
        "Q" = "Secondary metabolites biosynthesis",
        "R" = "General function prediction only",
        "S" = "Function unknown"
      )
      cog_labels <- sapply(used_cogs, function(c) {
        if (c %in% names(cog_descriptions)) {
          return(paste0("[", c, "] ", cog_descriptions[c]))
        }
        return(paste0("[", c, "]"))
      })
      lgd_cog <- ComplexHeatmap::Legend(
        labels = cog_labels,
        title = "COG Category",
        legend_gp = grid::gpar(fill = cog_colors[used_cogs]),
        ncol = 1
      )
    }
  }

  # RNA type legend -> 왼쪽 아래
  if (!is.null(rna_data) && nrow(rna_data) > 0) {
    rna_types_used <- unique(rna_data$type)
    rna_colors <- c(
      "tRNA" = "#E41A1C",
      "rRNA" = "#377EB8",
      "tmRNA" = "#4DAF4A",
      "ncRNA" = "#984EA3"
    )
    rna_colors_used <- rna_colors[rna_types_used[rna_types_used %in% names(rna_colors)]]
    if (length(rna_colors_used) > 0) {
      lgd_rna <- ComplexHeatmap::Legend(
        labels = names(rna_colors_used),
        title = "RNA Type",
        legend_gp = grid::gpar(fill = rna_colors_used),
        ncol = 1
      )
      lgd_list_left <- c(lgd_list_left, list(lgd_rna))
    }
  }

  # Mutation type legend -> 왼쪽 아래
  if (nrow(all_mutations) > 0) {
    used_types <- unique(all_mutations$type)
    used_types <- used_types[used_types %in% names(type_colors)]
    if (length(used_types) > 0) {
      lgd_type <- ComplexHeatmap::Legend(
        labels = used_types,
        title = "Mutation Type",
        legend_gp = grid::gpar(fill = type_colors[used_types]),
        ncol = 1
      )
      lgd_list_left <- c(lgd_list_left, list(lgd_type))
    }
  }

  # Sample legend (track background colors with transparency) -> 왼쪽 아래
  lgd_sample <- ComplexHeatmap::Legend(
    labels = rev(samples),  # Match track order
    title = "ALE Version (Track)",
    legend_gp = grid::gpar(fill = adjustcolor(rev(sample_colors), alpha.f = 0.3)),
    ncol = 1
  )
  lgd_list_left <- c(lgd_list_left, list(lgd_sample))

  # Missing coverage legend (MC type shown as border color) -> 왼쪽 아래
  if (!is.null(missing_data) && nrow(missing_data) > 0) {
    lgd_missing <- ComplexHeatmap::Legend(
      labels = c("MC_inversion (border)", "MC_repeat (border)"),
      title = "Missing Coverage Type",
      legend_gp = grid::gpar(fill = missing_colors),
      ncol = 1
    )
    lgd_list_left <- c(lgd_list_left, list(lgd_missing))
  }

  # Draw LEFT legends (Mutation, RNA, Sample, MC)
  if (length(lgd_list_left) > 0) {
    pd_left <- ComplexHeatmap::packLegend(list = lgd_list_left)
    ComplexHeatmap::draw(pd_left, x = grid::unit(4, "mm"), y = grid::unit(4, "mm"),
                         just = c("left", "bottom"))
  }

  # Draw RIGHT legend (COG only) - 오른쪽 아래
  if (!is.null(lgd_cog)) {
    ComplexHeatmap::draw(lgd_cog, x = grid::unit(1, "npc") - grid::unit(4, "mm"),
                         y = grid::unit(4, "mm"),
                         just = c("right", "bottom"))
  }

  grDevices::dev.off()
  circlize::circos.clear()

  return(output_file)
}
