#' Convert SnapGene .dna File to GenBank Format
#'
#' Uses SnapGene CLI to convert .dna files to .gbk format.
#' Requires SnapGene to be installed.
#'
#' @param dna_file Path to .dna file
#' @param output_file Optional output path (default: same name with .gbk extension)
#' @param timeout Timeout in seconds (default: 300)
#' @param snapgene_path Optional custom path to SnapGene executable
#' @return Path to the generated GenBank file
#' @export
convert_dna_to_genbank <- function(dna_file,
                                   output_file = NULL,
                                   timeout = 300,
                                   snapgene_path = NULL) {

  # Verify input file exists
  if (!file.exists(dna_file)) {
    stop("Input file not found: ", dna_file)
  }

  # If already a GenBank file, return as-is
  if (grepl("\\.(gb|gbk|genbank)$", dna_file, ignore.case = TRUE)) {
    message("File is already in GenBank format: ", dna_file)
    return(normalizePath(dna_file))
  }

  # Find SnapGene
  sg_path <- if (!is.null(snapgene_path)) snapgene_path else find_snapgene()

  # Set output file path
  if (is.null(output_file)) {
    output_file <- sub("\\.dna$", ".gbk", dna_file, ignore.case = TRUE)
  }

  # Skip if output already exists
  if (file.exists(output_file)) {
    message("Output file already exists, using existing: ", output_file)
    return(normalizePath(output_file))
  }

  # Ensure absolute paths
  dna_file <- normalizePath(dna_file)
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_file <- file.path(normalizePath(output_dir), basename(output_file))

  message("Converting ", basename(dna_file), " to GenBank format...")

  # Build command
  args <- c("--convert", "GenBank - SnapGene",
            "--input", dna_file,
            "--output", output_file)

  # Run SnapGene
  result <- tryCatch({
    if (requireNamespace("processx", quietly = TRUE)) {
      processx::run(sg_path, args = args, timeout = timeout,
                    error_on_status = FALSE)
    } else {
      # Fallback to system2
      status <- system2(sg_path, args = shQuote(args),
                       stdout = TRUE, stderr = TRUE, timeout = timeout)
      list(status = attr(status, "status") %||% 0, stdout = status)
    }
  }, error = function(e) {
    stop("SnapGene conversion failed: ", e$message)
  })

  # Verify output
  if (!file.exists(output_file)) {
    stop("Conversion failed: output file not created. ",
         "SnapGene may require a GUI session on some systems.")
  }

  message("Successfully converted to: ", output_file)
  return(normalizePath(output_file))
}

#' Generate Color for Evolution Stage
#'
#' Creates a color that fades from saturated magenta (early) to pale (late stages).
#'
#' @param stage_index Index of the current stage (0-based)
#' @param total_stages Total number of stages
#' @param base_hue Base hue value (default: 300 = magenta)
#' @return Hex color code
#' @export
get_stage_color <- function(stage_index, total_stages, base_hue = 300) {

  # Saturation decreases with stage (100% -> 20%)
  saturation <- max(20, 100 - stage_index * (80 / max(1, total_stages - 1)))

  # Convert HSV to RGB to Hex
  rgb_vals <- grDevices::hsv(base_hue / 360, saturation / 100, 1)
  return(rgb_vals)
}


#' Get Mutation Type Color for SnapGene
#'
#' Returns hex color code for mutation type (same as circos visualization).
#'
#' @param mutation_type Mutation type (INS, DEL, SNP, REP, UNKNOWN)
#' @return Hex color code
#' @export
get_mutation_color <- function(mutation_type) {
  colors <- c(
    "INS" = "#FF5C8D",   # Bright Pink
    "DEL" = "#03045E",   # Dark Navy
    "SNP" = "#139487",   # Teal
    "REP" = "#E5C87B",   # Golden Yellow
    "UNKNOWN" = "#808080", # Gray
    "insertion" = "#FF5C8D",
    "deletion" = "#03045E",
    "snp" = "#139487",
    "replacement" = "#E5C87B",
    "unknown" = "#808080"
  )

  type_upper <- toupper(mutation_type)
  type_lower <- tolower(mutation_type)

  if (type_upper %in% names(colors)) {
    return(colors[type_upper])
  } else if (type_lower %in% names(colors)) {
    return(colors[type_lower])
  }

  return("#808080")  # Default gray
}


#' Classify Mutation Type from Mutation String
#'
#' Determines mutation type (INS, DEL, SNP, REP) from breseq mutation string.
#'
#' @param mutation Mutation string from breseq
#' @return Mutation type as string
#' @export
classify_mutation <- function(mutation) {
  if (is.na(mutation) || is.null(mutation) || nchar(as.character(mutation)) == 0) {
    return("UNKNOWN")
  }

  mut <- as.character(mutation)

  # Insertion: +ATCG or (G)8→9 where new > old

  if (grepl("^\\+[ATCG]+$", mut, ignore.case = TRUE)) {
    return("INS")
  }

  # Homopolymer repeat change: (G)8→9
  if (grepl("\\([ATCG]+\\)\\d+.*→.*\\d+", mut, ignore.case = TRUE)) {
    nums <- as.numeric(regmatches(mut, gregexpr("\\d+", mut))[[1]])
    if (length(nums) >= 2) {
      if (nums[2] > nums[1]) return("INS")
      if (nums[2] < nums[1]) return("DEL")
    }
  }

  # Deletion: Δxxx bp or -ATCG
  if (grepl("^[ΔΔ]", mut) || grepl("^-[ATCG]+$", mut, ignore.case = TRUE)) {
    return("DEL")
  }

  # SNP: A→G (single base change)
  if (grepl("^[ATCG]→[ATCG]$", mut, ignore.case = TRUE)) {
    return("SNP")
  }

  # Replacement: xx bp→ATCG
  if (grepl("\\d+\\s*bp.*→", mut)) {
    return("REP")
  }

  return("UNKNOWN")
}

#' Format GenBank Qualifier Line
#'
#' Formats a qualifier value to fit within GenBank line width limits.
#' Uses standard GenBank format compatible with SnapGene.
#' - /label=value (NO quotes)
#' - /note="value" (WITH quotes)
#'
#' @param key Qualifier key (e.g., "/label", "/note")
#' @param value Qualifier value
#' @param max_width Maximum line width (default: 80)
#' @param indent Indentation for qualifier lines (default: 21)
#' @return Character vector of formatted lines
format_qualifier <- function(key, value, max_width = 80, indent = 21) {

  # Standard GenBank format: 21 spaces + /key=value or /key="value"
  prefix <- strrep(" ", indent)
  key_str <- paste0(key, "=")

  # Clean the value
  value <- trimws(as.character(value))
  if (nchar(value) == 0) return(character())

  # /label does NOT use quotes, /note DOES use quotes
  use_quotes <- !grepl("^/label", key)

  if (use_quotes) {
    # With quotes: /note="value"
    full_line <- paste0(prefix, key_str, '"', value, '"')
  } else {
    # Without quotes: /label=value
    full_line <- paste0(prefix, key_str, value)
  }

  if (nchar(full_line) <= max_width) {
    return(full_line)
  }

  # For long values, split across lines
  lines <- character()
  remaining <- value
  first_line <- TRUE

  while (nchar(remaining) > 0) {
    if (first_line) {
      if (use_quotes) {
        available <- max_width - nchar(prefix) - nchar(key_str) - 1
        chunk <- substr(remaining, 1, available)
        lines <- c(lines, paste0(prefix, key_str, '"', chunk))
      } else {
        available <- max_width - nchar(prefix) - nchar(key_str)
        chunk <- substr(remaining, 1, available)
        lines <- c(lines, paste0(prefix, key_str, chunk))
      }
      remaining <- substr(remaining, available + 1, nchar(remaining))
      first_line <- FALSE
    } else {
      available <- max_width - nchar(prefix)
      chunk <- substr(remaining, 1, available)

      if (nchar(remaining) <= available) {
        if (use_quotes) {
          lines <- c(lines, paste0(prefix, chunk, '"'))
        } else {
          lines <- c(lines, paste0(prefix, chunk))
        }
      } else {
        lines <- c(lines, paste0(prefix, chunk))
      }
      remaining <- substr(remaining, available + 1, nchar(remaining))
    }
  }

  return(lines)
}


#' Format SnapGene-compatible Color Qualifiers
#'
#' Creates color qualifiers in ApE format that SnapGene can import.
#' Uses both /ApEinfo_fwdcolor and /ApEinfo_revcolor for compatibility.
#'
#' @param color Hex color code (e.g., "#FF5C8D")
#' @param indent Indentation (default: 21)
#' @return Character vector of color qualifier lines
#' @keywords internal
format_snapgene_color <- function(color, indent = 21) {
  prefix <- strrep(" ", indent)

  # Convert hex to decimal for SnapGene native format
  # SnapGene uses decimal RGB values
  hex_clean <- gsub("^#", "", color)
  r <- strtoi(substr(hex_clean, 1, 2), 16)
  g <- strtoi(substr(hex_clean, 3, 4), 16)
  b <- strtoi(substr(hex_clean, 5, 6), 16)
  decimal_color <- r + g * 256 + b * 65536

  c(
    # ApE format (for ApE/SnapGene import)
    paste0(prefix, '/ApEinfo_fwdcolor="', color, '"'),
    paste0(prefix, '/ApEinfo_revcolor="', color, '"'),
    # Note format (backup)
    paste0(prefix, '/note="color: ', color, '"')
  )
}

#' Format Number with Comma Separators
#'
#' @param x Numeric value
#' @return Character string with commas (no extra spaces)
#' @keywords internal
format_number_with_commas <- function(x) {
  # Use prettyNum to avoid extra padding spaces from format()
  prettyNum(as.integer(x), big.mark = ",", scientific = FALSE)
}

#' Parse Mutation Range from breseq Data
#'
#' Extracts start and end positions from mutation data.
#'
#' @param position Position value from breseq
#' @param mutation Mutation string (for calculating deletion size, etc.)
#' @param end_position Optional explicit end position
#' @param mutation_type Type of mutation (from classify_mutation)
#' @return List with start and end positions
#' @export
parse_mutation_range <- function(position, mutation = NULL,
                                 end_position = NULL, mutation_type = NULL) {

  parsed <- parse_position(position)
  start <- parsed$start
  end <- parsed$end

  if (is.na(start)) return(list(ranges = list()))

  # If end position is explicitly provided

  if (!is.null(end_position) && !is.na(end_position)) {
    end_parsed <- parse_position(end_position)
    if (!is.na(end_parsed$start)) {
      end <- end_parsed$end
    }
  }

  # Infer end from mutation string if needed
  if (!is.null(mutation) && !is.na(mutation) && start == end) {
    mut <- as.character(mutation)

    # Deletion: Δ500 bp or Δ12,345 bp (with comma)
    if (grepl("Δ", mut)) {
      # Match digits with optional commas: 12,345 or 12345
      size_match <- regmatches(mut, regexpr("[0-9,]+", mut))
      if (length(size_match) > 0) {
        delta <- as.integer(gsub(",", "", size_match[1]))
        end <- start + delta - 1
      }
    }

    # Repeat change: (G)8→9
    if (grepl("\\([ATCG]+\\)\\d+.*→.*\\d+", mut, ignore.case = TRUE)) {
      nums <- as.numeric(regmatches(mut, gregexpr("\\d+", mut))[[1]])
      if (length(nums) >= 2) {
        repeat_unit <- regmatches(mut, regexpr("\\([ATCG]+\\)", mut, ignore.case = TRUE))
        unit_len <- nchar(gsub("[()]", "", repeat_unit))
        diff <- abs(nums[2] - nums[1])
        end <- start + (diff * unit_len) - 1
      }
    }

    # Insertion: +ATCG
    if (grepl("^\\+[ATCG]+$", mut, ignore.case = TRUE)) {
      inserted <- gsub("^\\+", "", mut)
      end <- start + nchar(inserted) - 1
    }
  }

  # Ensure start <= end
  if (end < start) {
    tmp <- start
    start <- end
    end <- tmp
  }

  return(list(ranges = list(c(start, end))))
}

#' Create Annotated GenBank File with Mutations
#'
#' Takes a reference sequence (GenBank or SnapGene .dna) and adds mutation
#' annotations from breseq comparison results. For multi-replicon genomes,
#' only mutations on the reference sequence (typically chromosome) are included;
#' plasmid mutations are ignored.
#'
#' @param reference_file Path to reference file (.dna or .gbk)
#' @param mutation_data Data frame with mutations (from compare_named_tables) or path to Excel file
#' @param output_file Path for output GenBank file
#' @param sheets Sheet names to process (default: c("Mutations", "Missing"))
#' @param color_by How to color features: "stage" (evolution round) or "type" (mutation type)
#' @param target_seq_id Optional: only include mutations from this seq_id.
#'   If NULL, auto-detected from reference file (first/largest replicon).
#' @return Path to the created GenBank file
#' @export
#' @examples
#' \dontrun{
#' # From data frame
#' create_annotated_genbank(
#'   reference_file = "reference.dna",
#'   mutation_data = comparison_results,
#'   output_file = "annotated.gbk"
#' )
#'
#' # From Excel file - only chromosome mutations
#' create_annotated_genbank(
#'   reference_file = "chromosome.dna",
#'   mutation_data = "mutations.xlsx",
#'   output_file = "annotated.gbk"
#' )
#' }
create_annotated_genbank <- function(reference_file,
                                     mutation_data,
                                     output_file,
                                     sheets = c("Mutations", "Missing"),
                                     color_by = c("stage", "type"),
                                     target_seq_id = NULL) {

  color_by <- match.arg(color_by)

  # Convert .dna to .gbk if needed
  if (grepl("\\.dna$", reference_file, ignore.case = TRUE)) {
    temp_gbk <- convert_dna_to_genbank(reference_file)
  } else {
    temp_gbk <- reference_file
  }

  # Get reference seq_id for filtering mutations
  ref_info <- tryCatch({
    parse_genbank_info(temp_gbk)
  }, error = function(e) NULL)

  if (is.null(target_seq_id) && !is.null(ref_info)) {
    target_seq_id <- ref_info$seq_id[1]
    message("Target replicon for GenBank export: ", target_seq_id)
  }

  # Read GenBank file
  gb_lines <- readLines(temp_gbk)

  # Remove CONTIG lines (causes SnapGene issues)
  contig_lines <- grep("^CONTIG", gb_lines)
  if (length(contig_lines) > 0) {
    message("Removing CONTIG lines for SnapGene compatibility")
    gb_lines <- gb_lines[-contig_lines]
  }

  # Find FEATURES and ORIGIN sections
  feature_start <- grep("^FEATURES", gb_lines)[1]
  origin_start <- grep("^ORIGIN", gb_lines)[1]

  if (is.na(feature_start) || is.na(origin_start)) {
    stop("Invalid GenBank format: missing FEATURES or ORIGIN section")
  }

  header_lines <- gb_lines[1:feature_start]
  feature_lines <- gb_lines[(feature_start + 1):(origin_start - 1)]
  origin_lines <- gb_lines[origin_start:length(gb_lines)]

  # Load mutation data
  if (is.character(mutation_data) && file.exists(mutation_data)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' required. Install with: install.packages('openxlsx')")
    }
    mutation_list <- lapply(sheets, function(sheet) {
      tryCatch(
        openxlsx::read.xlsx(mutation_data, sheet = sheet),
        error = function(e) NULL
      )
    })
    names(mutation_list) <- sheets
    mutation_list <- Filter(Negate(is.null), mutation_list)
  } else if (is.data.frame(mutation_data)) {
    mutation_list <- list(Mutations = mutation_data)
  } else if (is.list(mutation_data)) {
    mutation_list <- mutation_data
  } else {
    stop("mutation_data must be a data frame, list, or path to Excel file")
  }

  # Generate new feature lines
  new_features <- character()
  skipped_plasmid <- 0

  for (sheet_name in names(mutation_list)) {
    df <- mutation_list[[sheet_name]]
    if (is.null(df) || nrow(df) == 0) next

    message("Processing ", sheet_name, " (", nrow(df), " rows, ", ncol(df), " columns)")

    # Filter to target seq_id only (chromosome)
    seq_id_col <- find_seq_id_column(df)
    if (!is.null(seq_id_col) && !is.null(target_seq_id) && seq_id_col %in% names(df)) {
      original_rows <- nrow(df)
      df <- df[df[[seq_id_col]] == target_seq_id | is.na(df[[seq_id_col]]), ]
      skipped_plasmid <- skipped_plasmid + (original_rows - nrow(df))
    }

    # Find position column
    pos_col <- if ("position" %in% names(df)) "position" else
               if ("start" %in% names(df)) "start" else NULL

    if (is.null(pos_col)) {
      warning("No position column found in ", sheet_name, ". Skipping.")
      next
    }

    # Detect sample/stage columns (format: "mutation.SJY3/output/index" or "evidence.45SJY3")
    sample_pattern <- "\\.(45)?SJY\\d+|\\.[^/]+/output"
    sample_cols <- grep(sample_pattern, names(df), value = TRUE, ignore.case = TRUE)

    # Extract unique sample names
    extract_sample_name <- function(col) {
      # Try pattern: "mutation.SJY3/output/index" -> "SJY3"
      m <- regmatches(col, regexpr("(45)?SJY\\d+", col, ignore.case = TRUE))
      if (length(m) > 0 && nchar(m) > 0) return(m)
      # Try pattern: "mutation.sample_name" -> "sample_name"
      parts <- strsplit(col, "\\.")[[1]]
      if (length(parts) >= 2) return(sub("/.*", "", parts[2]))
      return(NA)
    }

    stage_names <- unique(na.omit(sapply(sample_cols, extract_sample_name)))
    stage_names <- stage_names[order(as.numeric(gsub("[^0-9]", "", stage_names)))]

    message("  Detected stages: ", paste(stage_names, collapse = ", "))

    # Find mutation columns (any column starting with "mutation")
    mut_cols <- grep("^mutation", names(df), value = TRUE, ignore.case = TRUE)
    # Find annotation columns
    ann_cols <- grep("^annotation", names(df), value = TRUE, ignore.case = TRUE)
    # Find gene columns
    gene_cols <- grep("^gene", names(df), value = TRUE, ignore.case = TRUE)

    message("  Found ", length(mut_cols), " mutation cols, ", length(ann_cols), " annotation cols")

    # Process each row
    for (i in seq_len(nrow(df))) {
      row <- df[i, ]
      position <- row[[pos_col]]

      if (is.na(position)) next

      # Remove commas from position
      position <- as.numeric(gsub(",", "", as.character(position)))
      if (is.na(position)) next

      # Get mutation info from first non-NA mutation column
      mutation <- NA
      first_stage <- NA

      for (mut_col in mut_cols) {
        val <- row[[mut_col]]
        if (!is.null(val) && !is.na(val) && nchar(as.character(val)) > 0) {
          mutation <- as.character(val)
          first_stage <- extract_sample_name(mut_col)
          break
        }
      }

      # Get annotation from first non-NA annotation column
      annotation <- "Unknown"
      for (ann_col in ann_cols) {
        val <- row[[ann_col]]
        if (!is.null(val) && !is.na(val) && nchar(as.character(val)) > 0) {
          annotation <- as.character(val)
          break
        }
      }

      # Also try mutation_type column if it exists
      if ("mutation_type" %in% names(row)) {
        mut_type <- as.character(row[["mutation_type"]])
        if (is.na(mut_type) || nchar(mut_type) == 0) {
          mut_type <- if (!is.na(mutation)) classify_mutation(mutation) else "unknown"
        }
      } else {
        mut_type <- if (!is.na(mutation)) classify_mutation(mutation) else "unknown"
      }

      # Handle Missing sheet (end position)
      end_position <- position
      if (sheet_name == "Missing" || sheet_name == "Unassigned missing coverage evidence") {
        end_cols <- grep("^end", names(df), value = TRUE, ignore.case = TRUE)
        for (end_col in end_cols) {
          val <- row[[end_col]]
          if (!is.null(val) && !is.na(val)) {
            end_val <- as.numeric(gsub(",", "", as.character(val)))
            if (!is.na(end_val) && end_val > position) {
              end_position <- end_val
              break
            }
          }
        }
      }

      # Parse mutation range
      range_info <- parse_mutation_range(position, mutation, end_position, mut_type)

      if (length(range_info$ranges) == 0) next

      # Determine color
      if (color_by == "stage" && !is.na(first_stage)) {
        stage_idx <- match(first_stage, stage_names)
        if (is.na(stage_idx)) stage_idx <- 1
        color <- get_stage_color(stage_idx - 1, length(stage_names))
      } else {
        color <- get_mutation_color(mut_type)
      }

      # Add feature for each range
      stage_str <- if(!is.na(first_stage)) first_stage else "Unknown"

      for (rng in range_info$ranges) {
        start_pos <- as.integer(rng[1])
        end_pos <- as.integer(rng[2])

        if (is.na(start_pos) || is.na(end_pos)) next
        if (start_pos > end_pos) {
          tmp <- start_pos; start_pos <- end_pos; end_pos <- tmp
        }

        if (grepl("Missing", sheet_name, ignore.case = TRUE)) {
          delta <- end_pos - start_pos + 1
          # Format: [MC:SJY3] (13,543 bp) annotation - no Δ so SnapGene won't parse it
          label <- sprintf("[MC:%s] (%s bp) %s", stage_str, format_number_with_commas(delta), annotation)
        } else {
          mut_str <- if (!is.na(mutation)) as.character(mutation) else ""
          # Keep commas in label for readability (coordinates are calculated separately)
          # Format: [SNP:SJY3] A→G annotation  or  [DEL:SJY3] Δ13,543 bp annotation
          label <- sprintf("[%s:%s] %s %s", mut_type, stage_str, mut_str, annotation)
        }

        # Clean label
        label <- gsub("\\s+", " ", label)
        label <- trimws(label)

        # Build feature entry (Python format for SnapGene)
        feature_entry <- c(
          sprintf("     misc_feature    %d..%d", start_pos, end_pos),
          format_qualifier("/label", label),
          format_snapgene_color(color)
        )

        new_features <- c(new_features, feature_entry)
      }
    }
  }

  # Combine all sections
  all_lines <- c(header_lines, feature_lines, new_features, origin_lines)

  # Write output
  output_dir <- dirname(output_file)
  if (nchar(output_dir) > 0 && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  writeLines(all_lines, output_file)

  # Report results
  n_features <- length(grep("^     misc_feature", new_features))
  message("Annotated GenBank file saved: ", output_file)
  message("  Added ", n_features, " mutation features")
  if (skipped_plasmid > 0) {
    message("  Skipped ", skipped_plasmid, " plasmid mutations (chromosome only)")
  }

  return(normalizePath(output_file))
}


#' Export Mutations to SnapGene-compatible GenBank
#'
#' Standalone function to add breseq mutation features to a GenBank file for
#' visualization in SnapGene. This is separate from the main pipeline and can
#' be run independently.
#'
#' Features:
#' - Converts SnapGene .dna files automatically via CLI
#' - Adds mutation features with type-based colors (INS, DEL, SNP, REP)
#' - Adds Missing Coverage (MC) regions
#' - Preserves original GenBank annotations
#' - Auto-detects breseq_analysis folder for input/output
#'
#' Color scheme (same as Circos):
#' - INS (Insertion): #FF5C8D (Pink)
#' - DEL (Deletion): #03045E (Dark Navy)
#' - SNP: #139487 (Teal)
#' - REP (Replacement): #E5C87B (Golden Yellow)
#'
#' @param reference_file Path to reference file (.dna or .gbk)
#' @param mutation_xlsx Path to Excel file (default: auto-detect breseq_analysis/merged_mutation_results.xlsx)
#' @param output_file Path for output GenBank file (default: breseq_analysis/<ref>_with_mutations.gbk)
#' @param sheets Sheet names to process (default: c("Mutations", "Missing"))
#' @param color_by How to color features: "type" (mutation type) or "stage" (ALE version)
#' @return Path to the created GenBank file
#' @export
#' @examples
#' \dontrun{
#' # Simplest usage - auto-detect mutation file from breseq_analysis folder
#' export_mutations_to_genbank(
#'   reference_file = "reference.gbk"
#' )
#'
#' # From SnapGene .dna file (auto-detects mutation xlsx)
#' export_mutations_to_genbank(
#'   reference_file = "reference.dna"
#' )
#'
#' # Color by ALE stage instead of mutation type
#' export_mutations_to_genbank(
#'   reference_file = "reference.gbk",
#'   color_by = "stage"
#' )
#'
#' # Manual paths
#' export_mutations_to_genbank(
#'   reference_file = "reference.gbk",
#'   mutation_xlsx = "custom_path/mutations.xlsx",
#'   output_file = "output/annotated.gbk"
#' )
#' }
export_mutations_to_genbank <- function(reference_file,
                                        mutation_xlsx = NULL,
                                        output_file = NULL,
                                        sheets = c("Mutations", "Missing"),
                                        color_by = c("type", "stage")) {

  color_by <- match.arg(color_by)

  # Check if SnapGene is running (macOS)
  snapgene_running <- tryCatch({
    result <- system("pgrep -x SnapGene", intern = TRUE, ignore.stderr = TRUE)
    length(result) > 0
  }, error = function(e) FALSE)

  if (snapgene_running) {
    message("\n⚠️  SnapGene이 실행 중입니다. CLI 변환을 위해 종료해야 할 수 있습니다.")
    answer <- readline("SnapGene을 종료할까요? [Y/n]: ")
    if (tolower(answer) != "n") {
      message("SnapGene 종료 중...")
      system("pkill -x SnapGene", ignore.stderr = TRUE)
      Sys.sleep(2)  # Wait for clean shutdown
      message("SnapGene이 종료되었습니다.")
    } else {
      message("SnapGene을 종료하지 않고 계속 진행합니다.")
    }
  }

  # Validate reference file
  if (!file.exists(reference_file)) {
    stop("Reference file not found: ", reference_file)
  }

  # Get reference directory for searching

  ref_dir <- dirname(normalizePath(reference_file))

  # Auto-detect mutation xlsx if not provided
  if (is.null(mutation_xlsx)) {
    # Search for breseq_analysis folder containing merged_mutation_results.xlsx
    search_paths <- c(
      file.path(ref_dir, "breseq_analysis", "merged_mutation_results.xlsx"),
      file.path(ref_dir, "breseq_analysis", "mutant_results.xlsx"),
      file.path(dirname(ref_dir), "breseq_analysis", "merged_mutation_results.xlsx"),
      file.path(dirname(ref_dir), "breseq_analysis", "mutant_results.xlsx"),
      file.path(getwd(), "breseq_analysis", "merged_mutation_results.xlsx"),
      file.path(getwd(), "breseq_analysis", "mutant_results.xlsx")
    )

    for (path in search_paths) {
      if (file.exists(path)) {
        mutation_xlsx <- normalizePath(path)
        message("Auto-detected mutation file: ", mutation_xlsx)
        break
      }
    }

    if (is.null(mutation_xlsx)) {
      stop("Could not find mutation Excel file. Searched in:\n",
           "  - breseq_analysis/merged_mutation_results.xlsx\n",
           "  - breseq_analysis/mutant_results.xlsx\n",
           "Please provide mutation_xlsx parameter explicitly.")
    }
  }

  if (!file.exists(mutation_xlsx)) {
    stop("Mutation Excel file not found: ", mutation_xlsx)
  }

  # Determine breseq_analysis folder for output
  breseq_analysis_dir <- dirname(mutation_xlsx)

  message("\n", strrep("=", 60))
  message("Export Mutations to GenBank for SnapGene")
  message(strrep("=", 60))
  message("Reference: ", basename(reference_file))
  message("Mutations: ", basename(mutation_xlsx))
  message("Color by: ", color_by)

  # Convert .dna to .gbk if needed
  if (grepl("\\.dna$", reference_file, ignore.case = TRUE)) {
    message("\nConverting SnapGene .dna to GenBank...")
    temp_gbk <- convert_dna_to_genbank(reference_file)
    cleanup_temp <- TRUE
  } else {
    temp_gbk <- reference_file
    cleanup_temp <- FALSE
  }

  # Set output file name - default to breseq_analysis folder
  if (is.null(output_file)) {
    base_name <- sub("\\.(dna|gb|gbk|genbank)$", "", basename(reference_file), ignore.case = TRUE)
    output_file <- file.path(breseq_analysis_dir, paste0(base_name, "_with_mutations.gbk"))
  }

  # Read GenBank file
  gb_lines <- readLines(temp_gbk)

  # Remove CONTIG lines (causes SnapGene issues)
  contig_lines <- grep("^CONTIG", gb_lines)
  if (length(contig_lines) > 0) {
    message("Removing CONTIG lines for SnapGene compatibility")
    gb_lines <- gb_lines[-contig_lines]
  }

  # Find FEATURES and ORIGIN sections
  feature_start <- grep("^FEATURES", gb_lines)[1]
  origin_start <- grep("^ORIGIN", gb_lines)[1]

  if (is.na(feature_start) || is.na(origin_start)) {
    stop("Invalid GenBank format: missing FEATURES or ORIGIN section")
  }

  header_lines <- gb_lines[1:feature_start]
  existing_features <- gb_lines[(feature_start + 1):(origin_start - 1)]
  origin_lines <- gb_lines[origin_start:length(gb_lines)]

  # Load mutation data from Excel
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' required. Install with: install.packages('openxlsx')")
  }

  mutation_list <- lapply(sheets, function(sheet) {
    tryCatch(
      openxlsx::read.xlsx(mutation_xlsx, sheet = sheet),
      error = function(e) {
        message("  Sheet '", sheet, "' not found, skipping...")
        NULL
      }
    )
  })
  names(mutation_list) <- sheets
  mutation_list <- Filter(Negate(is.null), mutation_list)

  if (length(mutation_list) == 0) {
    stop("No valid sheets found in Excel file")
  }

  # Process mutations
  new_features <- character()
  mutation_count <- 0
  mc_count <- 0

  for (sheet_name in names(mutation_list)) {
    df <- mutation_list[[sheet_name]]
    if (is.null(df) || nrow(df) == 0) next

    message("\nProcessing ", sheet_name, " (", nrow(df), " rows)...")

    # Find position column
    pos_col <- if ("position" %in% names(df)) "position" else
               if ("start" %in% names(df)) "start" else NULL
    if (is.null(pos_col)) {
      message("  No position column found, skipping...")
      next
    }

    # Count for this sheet
    sheet_processed <- 0
    sheet_skipped <- 0

    # Detect stage columns for stage-based coloring
    sample_cols <- grep("\\.(45)?SJY\\d+|\\.[^/]+/output", names(df), value = TRUE, ignore.case = TRUE)
    extract_sample <- function(col) {
      m <- regmatches(col, regexpr("(45)?SJY\\d+", col, ignore.case = TRUE))
      if (length(m) > 0 && nchar(m) > 0) return(m)
      parts <- strsplit(col, "\\.")[[1]]
      if (length(parts) >= 2) return(sub("/.*", "", parts[2]))
      return(NA)
    }
    stage_names <- unique(na.omit(sapply(sample_cols, extract_sample)))
    stage_names <- stage_names[order(as.numeric(gsub("[^0-9]", "", stage_names)))]

    if (length(stage_names) > 0) {
      message("  Detected stages: ", paste(stage_names, collapse = " -> "))
    }

    # Find mutation/annotation columns
    mut_cols <- grep("^mutation", names(df), value = TRUE, ignore.case = TRUE)
    ann_cols <- grep("^annotation", names(df), value = TRUE, ignore.case = TRUE)

    for (i in seq_len(nrow(df))) {
      row <- df[i, ]
      position <- as.numeric(gsub(",", "", as.character(row[[pos_col]])))
      if (is.na(position)) {
        sheet_skipped <- sheet_skipped + 1
        next
      }

      # Get mutation info (first non-NA)
      mutation <- NA
      first_stage <- NA
      for (mut_col in mut_cols) {
        val <- row[[mut_col]]
        if (!is.null(val) && !is.na(val) && nchar(as.character(val)) > 0) {
          mutation <- as.character(val)
          first_stage <- extract_sample(mut_col)
          break
        }
      }

      # Get annotation (first non-NA)
      annotation <- ""
      for (ann_col in ann_cols) {
        val <- row[[ann_col]]
        if (!is.null(val) && !is.na(val) && nchar(as.character(val)) > 0) {
          annotation <- as.character(val)
          break
        }
      }

      # Classify mutation type
      mut_type <- classify_mutation(mutation)

      # Calculate end position
      end_position <- position
      if (grepl("Missing", sheet_name, ignore.case = TRUE)) {
        # For Missing sheet, look for end column
        end_cols <- grep("^end", names(df), value = TRUE, ignore.case = TRUE)
        for (end_col in end_cols) {
          val <- row[[end_col]]
          if (!is.null(val) && !is.na(val)) {
            # Handle range format (XXX-YYY)
            val_str <- gsub(",", "", as.character(val))
            val_str <- gsub("–", "-", val_str)  # em-dash to hyphen
            if (grepl("-", val_str)) {
              parts <- strsplit(val_str, "-")[[1]]
              end_position <- as.numeric(parts[length(parts)])
            } else {
              end_position <- as.numeric(val_str)
            }
            if (!is.na(end_position) && end_position > position) break
          }
        }
      } else {
        # Calculate from mutation string
        range_info <- parse_mutation_range(position, mutation)
        if (length(range_info$ranges) > 0) {
          end_position <- range_info$ranges[[1]][2]
        }
      }

      # Ensure valid range
      if (is.na(end_position) || end_position < position) {
        end_position <- position
      }

      # Determine color
      if (color_by == "type") {
        color <- get_mutation_color(mut_type)
      } else if (!is.na(first_stage) && length(stage_names) > 0) {
        stage_idx <- match(first_stage, stage_names)
        if (is.na(stage_idx)) stage_idx <- 1
        color <- get_stage_color(stage_idx - 1, length(stage_names))
      } else {
        color <- get_mutation_color(mut_type)
      }

      # Create label with stage visible
      stage_str <- if(!is.na(first_stage)) first_stage else "Unknown"
      is_mc <- grepl("Missing", sheet_name, ignore.case = TRUE)
      if (is_mc) {
        delta <- end_position - position + 1
        # Format: [MC:SJY3] (13,543 bp) annotation - no Δ so SnapGene won't parse it
        label <- sprintf("[MC:%s] (%s bp) %s", stage_str, format_number_with_commas(delta), annotation)
        mc_count <- mc_count + 1
      } else {
        mut_str <- if (!is.na(mutation)) as.character(mutation) else ""
        # Keep commas in label for readability (coordinates are calculated separately)
        # Format: [SNP:SJY3] A→G annotation  or  [DEL:SJY3] Δ13,543 bp annotation
        label <- sprintf("[%s:%s] %s %s", mut_type, stage_str, mut_str, annotation)
        mutation_count <- mutation_count + 1
      }
      label <- trimws(gsub("\\s+", " ", label))

      # Build GenBank feature entry (SnapGene compatible format)
      # ApE color format for SnapGene import
      feature_entry <- c(
        sprintf("     misc_feature    %d..%d", as.integer(position), as.integer(end_position)),
        format_qualifier("/label", label),
        format_snapgene_color(color)
      )

      new_features <- c(new_features, feature_entry)
      sheet_processed <- sheet_processed + 1
    }

    message("  -> Processed: ", sheet_processed, ", Skipped (no position): ", sheet_skipped)
  }

  # Combine all sections
  all_lines <- c(header_lines, existing_features, new_features, origin_lines)

  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (nchar(output_dir) > 0 && output_dir != "." && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Write output file
  writeLines(all_lines, output_file)

  # Summary
  message("\n", strrep("=", 60))
  message("GenBank file created: ", output_file)
  message("  Mutations added: ", mutation_count)
  message("  Missing Coverage added: ", mc_count)
  message("  Total features: ", mutation_count + mc_count)
  message(strrep("=", 60))
  message("\nOpen in SnapGene to visualize mutations on the genome map!")

  return(normalizePath(output_file))
}
