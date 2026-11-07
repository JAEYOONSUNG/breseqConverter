# breseqConverter standalone legacy workflow (2026-11-07)
# 사용법: 아래 설정을 바꾼 뒤 기존처럼 실행하세요. breseqConverter 패키지 설치는 필요 없습니다.
# Use the existing workflow; this file adds no automatic rendering call.
# HIDE_INTERGENIC_LABELS=TRUE hides intergenic labels only; all mutation marks remain.
# GC colors follow the original reference figure: skew green/red, ratio blue/yellow.
# Positive/negative mean above/below the existing *_minus_average baseline.
# GC formula, averaging, coordinates, mutations, and cumulative tracks are unchanged.
# The original workbook convention (*mutant_results.xlsx in getwd()) and libraries
# below still apply. Prepare the same GenomeDrawer globals before rendering:
# gff, contig_length, <seqid>_length, nucleotide, window_<seqid>, Genome_summary,
# RNA, mutation_list, samples, and version_colors (genbank_table is optional).
# Call breseqVisualization(...) as before; both options can also be overridden there.
HIDE_INTERGENIC_LABELS <- FALSE
GC_COLORS <- c(
  skew_positive = "#00FF00", skew_negative = "#FF0000",
  ratio_positive = "#0000FF", ratio_negative = "#FFFF00"
)

.breseq_one_shot_filter_intergenic_labels <- function(data, hide_intergenic_labels = FALSE) {
  if (!isTRUE(hide_intergenic_labels) || !"annotation" %in% names(data)) return(data)
  annotation <- trimws(gsub("\u00a0", " ", as.character(data$annotation), fixed = TRUE))
  intergenic <- !is.na(annotation) &
    grepl("^intergenic([[:space:]]|$)", annotation, ignore.case = TRUE)
  data[!intergenic, , drop = FALSE]
}

.breseq_one_shot_gc_colors <- function(gc_colors = NULL) {
  defaults <- c(skew_positive = "#00FF00", skew_negative = "#FF0000",
                ratio_positive = "#0000FF", ratio_negative = "#FFFF00")
  if (is.null(gc_colors)) return(defaults)
  if (!is.character(gc_colors) || !length(gc_colors) || is.null(names(gc_colors)) ||
      anyNA(gc_colors) || anyNA(names(gc_colors)) ||
      any(!nzchar(gc_colors)) || any(!nzchar(names(gc_colors))) || anyDuplicated(names(gc_colors))) {
    stop("gc_colors must be a named character vector with unique, non-missing entries.")
  }
  unknown <- setdiff(names(gc_colors), names(defaults))
  if (length(unknown)) stop("Unknown gc_colors entries: ", paste(unknown, collapse = ", "))
  tryCatch(grDevices::col2rgb(gc_colors), error = function(e) {
    stop("gc_colors contains an invalid R color: ", conditionMessage(e))
  })
  defaults[names(gc_colors)] <- gc_colors
  defaults
}

.breseq_one_shot_gc_track_colors <- function(type, deviation, gc_colors = NULL) {
  palette <- .breseq_one_shot_gc_colors(gc_colors)
  prefix <- switch(type, gc_skew = "skew", gc_ratio = "ratio",
                   stop("Unknown GC track type: ", type))
  ifelse(deviation >= 0, palette[[paste0(prefix, "_positive")]],
         palette[[paste0(prefix, "_negative")]])
}

# 필요한 라이브러리 로드
library(circlize)
library(readr)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(rlang)
library(gtools)
library(reshape2)
library(ComplexHeatmap)
library(stringr)

# Set the working directory
breseq_result_dir <- getwd()

# List breseq result files, excluding temporary files
breseq_result_files <- list.files(
  path = breseq_result_dir,
  pattern = "^[^~].*mutant_results\\.xlsx$",
  full.names = TRUE
)

# Handle cases with multiple or zero matching files
if (length(breseq_result_files) > 1) {
  warning("Multiple mutant_results.xlsx files found. Only the first file will be used.")
} else if (length(breseq_result_files) == 0) {
  stop("No mutant_results.xlsx file found in the directory.")
}

# Use the first matching file
mutation_excel_file_name <- breseq_result_files[1]
cat("Using breseq result file:", mutation_excel_file_name, "\n")

# Define sheet names to extract
sheet_names <- c("Mutations", "Missing", "New junction")

# Read data from each sheet into a list
mutation_data_list <- list()
for (sheet in sheet_names) {
  mutation_data_list[[sheet]] <- readxl::read_excel(mutation_excel_file_name, sheet = sheet)
}

# Assign each sheet to a separate variable
mutations_df <- mutation_data_list[["Mutations"]]
missing_df <- mutation_data_list[["Missing"]]
new_junction_df <- mutation_data_list[["New junction"]]

###############################################################################################################################################
# Process mutation data from mutations_df
if (!is.null(mutations_df)) {
  message("Processing mutation data from mutations_df...")

  # Clean position column
  mutations_df$position <- as.numeric(gsub(",", "", mutations_df$position))

  # Identify columns for reshaping
  base_cols <- c("position")
  value_cols <- c("evidence", "mutation", "annotation", "gene", "description")
  all_cols <- colnames(mutations_df)
  tag_cols <- setdiff(all_cols, c(base_cols, value_cols))

  # Function to extract group tag
  get_group_tag <- function(col_name) {
    if (col_name %in% value_cols) {
      return("None")
    }
    parts <- unlist(strsplit(col_name, "\\."))
    if (length(parts) > 1) {
      tag <- sub("/output/index", "", parts[2])
      return(tag)
    }
    return("Unknown")
  }

  # Melt the data to long format
  melted <- reshape2::melt(
    mutations_df,
    id.vars = "position",
    measure.vars = setdiff(all_cols, "position"),
    variable.name = "variable",
    value.name = "value"
  )

  # Extract group tag and variable type
  melted <- dplyr::mutate(
    melted,
    ALE_version = sapply(as.character(variable), get_group_tag),
    variable_type = dplyr::case_when(
      grepl("^evidence", variable) ~ "evidence",
      grepl("^mutation", variable) ~ "mutation",
      grepl("^annotation", variable) ~ "annotation",
      grepl("^gene", variable) ~ "gene",
      grepl("^description", variable) ~ "description",
      TRUE ~ "other"
    )
  )

  # Filter out rows with NA values and non-relevant variable types
  melted <- dplyr::filter(melted, !is.na(value), variable_type %in% value_cols)

  # Cast to wide format
  mutation_casted <- reshape2::dcast(
    melted,
    position + ALE_version ~ variable_type,
    value.var = "value",
    fun.aggregate = function(x) paste(x, collapse = ";")
  )

  # Standardize columns
  required_cols <- c("position", "ALE_version", "mutation", "annotation", "gene", "description")
  for (col in required_cols) {
    if (!col %in% colnames(mutation_casted)) {
      mutation_casted[[col]] <- NA
    }
  }

  # Select and add label column
  mutation_data <- dplyr::select(mutation_casted, dplyr::all_of(required_cols)) %>%
    dplyr::filter(!(is.na(mutation) & is.na(annotation) & is.na(gene) & is.na(description))) %>%
    dplyr::mutate(
      mutation = as.character(mutation),
      annotation = as.character(annotation),
      label = paste(mutation, annotation, sep = " | ")
    ) %>%
    as.data.frame()

  # Identify new mutations
  samples <- tag_cols %>% sub("^[^.]+\\.([^/]+?)/.*$", "\\1", .) %>% unique()
  samples <- samples[samples %in% unique(mutation_data$ALE_version)]
  mutation_list <- list()
  previous_mutations <- data.frame(position = numeric(), mutation = character(), stringsAsFactors = FALSE)

  for (sample in samples) {
    data <- dplyr::filter(mutation_data, ALE_version == sample) %>%
      dplyr::select(position, ALE_version, mutation, annotation, gene, description, label) %>%
      as.data.frame()

    if (nrow(data) == 0) {
      cat("Warning: No mutations for sample", sample, "\n")
      next
    }

    data$new <- !paste(data$position, data$mutation, sep = ":") %in%
      paste(previous_mutations$position, previous_mutations$mutation, sep = ":")

    # Define colors for ALE_version
    ALE_versions <- unique(mutation_data$ALE_version)
    ALE_versions <- ALE_versions[gtools::mixedorder(ALE_versions)]
    n_versions <- length(ALE_versions)
    version_colors <- setNames(
      RColorBrewer::brewer.pal(min(n_versions, 8), "Set2")[1:n_versions],
      ALE_versions
    )
    assign("version_colors", version_colors, envir = .GlobalEnv)

    # Add version_color column
    data <- data %>%
      dplyr::mutate(
        version_color = version_colors[ALE_version]
      )

    mutation_list[[sample]] <- data
    previous_mutations <- rbind(previous_mutations, data[, c("position", "mutation")])
  }

  # Ensure all mutation_list have new column as logical
  if (length(mutation_list) > 0) {
    mutation_list <- lapply(mutation_list, function(df) {
      if (!"new" %in% colnames(df)) {
        df$new <- FALSE
      }
      df$new <- as.logical(df$new)
      return(df)
    })
  }
} # mutations_df 처리 종료


###############################################################################################################################################
# Process missing data from missing_df
if (!is.null(missing_df)) {
  message("Processing missing data from missing_df...")

  # Identify columns for reshaping
  base_cols <- c("start")
  value_cols <- c("Forward aligned reads", "Reverse aligned reads", "Read Coverage Depth",
                  "seq id", "end", "size", "←reads", "reads→", "gene", "description")
  all_cols <- colnames(missing_df)
  tag_cols <- setdiff(all_cols, c(base_cols, value_cols))

  # Function to extract group tag
  get_group_tag <- function(col_name) {
    if (col_name %in% value_cols) {
      return("None")
    }
    parts <- unlist(strsplit(col_name, "\\."))
    if (length(parts) > 1) {
      tag <- sub("/output/index", "", parts[2])
      return(tag)
    }
    return("Unknown")
  }

  # Melt the data to long format
  melted <- reshape2::melt(
    missing_df,
    id.vars = "start",
    measure.vars = setdiff(all_cols, "start"),
    variable.name = "variable",
    value.name = "value"
  )

  # Extract group tag and variable type
  melted <- dplyr::mutate(
    melted,
    ALE_version = sapply(as.character(variable), get_group_tag),
    variable_type = dplyr::case_when(
      grepl("^Forward aligned reads", variable) ~ "Forward aligned reads",
      grepl("^Reverse aligned reads", variable) ~ "Reverse aligned reads",
      grepl("^Read Coverage Depth", variable) ~ "Read Coverage Depth",
      grepl("^seq id", variable) ~ "seq id",
      grepl("^end", variable) ~ "end",
      grepl("^size", variable) ~ "size",
      grepl("^←reads", variable) ~ "←reads",
      grepl("^reads→", variable) ~ "reads→",
      grepl("^gene", variable) ~ "gene",
      grepl("^description", variable) ~ "description",
      TRUE ~ "other"
    )
  )

  # Filter out rows with NA values and non-relevant variable types
  melted <- dplyr::filter(melted, !is.na(value), variable_type %in% value_cols)

  # Cast to wide format
  missing_casted <- reshape2::dcast(
    melted,
    start + ALE_version ~ variable_type,
    value.var = "value",
    fun.aggregate = function(x) paste(x, collapse = ";")
  )

  # Standardize columns
  required_cols <- c("start", "Forward aligned reads", "Reverse aligned reads",
                     "Read Coverage Depth", "seq id", "end", "size", "←reads",
                     "reads→", "gene", "description")
  for (col in required_cols) {
    if (!col %in% colnames(missing_casted)) {
      missing_casted[[col]] <- NA
    }
  }

  # Function to parse coordinate ranges (e.g., "824661–826118" or "34385")
  parse_coordinate_range <- function(coord_str) {
    if (is.na(coord_str) || coord_str == "" || is.null(coord_str)) {
      return(list(min = NA, max = NA))
    }

    coord_str <- as.character(coord_str)

    # Handle single coordinate (no range)
    if (!grepl("–", coord_str) && !grepl("-", coord_str)) {
      coord <- as.numeric(coord_str)
      if (is.na(coord)) {
        return(list(min = NA, max = NA))
      }
      return(list(min = coord, max = coord))
    }

    # Handle range coordinates
    if (grepl("–", coord_str)) {
      coords <- unlist(strsplit(coord_str, "–"))
    } else if (grepl("-", coord_str)) {
      coords <- unlist(strsplit(coord_str, "-"))
    } else {
      return(list(min = NA, max = NA))
    }

    if (length(coords) == 2) {
      min_val <- as.numeric(coords[1])
      max_val <- as.numeric(coords[2])

      if (is.na(min_val) || is.na(max_val)) {
        return(list(min = NA, max = NA))
      }

      return(list(min = min_val, max = max_val))
    }

    return(list(min = NA, max = NA))
  }

  # Function to determine region type (inversion vs repeat)
  determine_region_type <- function(start_min, start_max, end_min, end_max) {
    if (!is.na(end_min) && !is.na(end_max) && end_min > end_max) {
      return("MC_inversion")
    } else {
      return("MC_repeat")
    }
  }

  # Function to calculate certain and uncertain regions
  calculate_regions <- function(start_min, start_max, end_min, end_max, region_type) {
    start_is_single <- (!is.na(start_min) && !is.na(start_max) && start_min == start_max)
    end_is_single <- (!is.na(end_min) && !is.na(end_max) && end_min == end_max)

    if (region_type == "MC_inversion") {
      certain_start <- start_max
      certain_end <- end_max

      if (start_is_single) {
        uncertain_left_start <- NA
        uncertain_left_end <- NA
      } else {
        uncertain_left_start <- start_min
        uncertain_left_end <- start_max
      }

      if (end_is_single) {
        uncertain_right_start <- NA
        uncertain_right_end <- NA
      } else {
        uncertain_right_start <- end_max
        uncertain_right_end <- end_min
      }

    } else {
      certain_start <- start_max
      certain_end <- end_min

      if (start_is_single) {
        uncertain_left_start <- NA
        uncertain_left_end <- NA
      } else {
        uncertain_left_start <- start_min
        uncertain_left_end <- start_max
      }

      if (end_is_single) {
        uncertain_right_start <- NA
        uncertain_right_end <- NA
      } else {
        uncertain_right_start <- end_min
        uncertain_right_end <- end_max
      }
    }

    return(list(
      certain_start = certain_start,
      certain_end = certain_end,
      uncertain_left_start = uncertain_left_start,
      uncertain_left_end = uncertain_left_end,
      uncertain_right_start = uncertain_right_start,
      uncertain_right_end = uncertain_right_end,
      start_is_single = start_is_single,
      end_is_single = end_is_single
    ))
  }

  # Process missing regions
  message("Processing missing regions...")

  # Initialize new columns
  missing_casted$start_min <- NA_real_
  missing_casted$start_max <- NA_real_
  missing_casted$end_min <- NA_real_
  missing_casted$end_max <- NA_real_

  # Parse coordinates
  for (i in 1:nrow(missing_casted)) {
    start_result <- parse_coordinate_range(missing_casted$start[i])
    missing_casted$start_min[i] <- start_result$min
    missing_casted$start_max[i] <- start_result$max

    end_result <- parse_coordinate_range(missing_casted$end[i])
    missing_casted$end_min[i] <- end_result$min
    missing_casted$end_max[i] <- end_result$max
  }

  # Determine region types and calculate regions
  missing_casted$region_type <- character(nrow(missing_casted))
  missing_casted$certain_start <- NA_real_
  missing_casted$certain_end <- NA_real_
  missing_casted$uncertain_left_start <- NA_real_
  missing_casted$uncertain_left_end <- NA_real_
  missing_casted$uncertain_right_start <- NA_real_
  missing_casted$uncertain_right_end <- NA_real_
  missing_casted$start_is_single <- logical(nrow(missing_casted))
  missing_casted$end_is_single <- logical(nrow(missing_casted))

  for (i in 1:nrow(missing_casted)) {
    missing_casted$region_type[i] <- determine_region_type(
      missing_casted$start_min[i], missing_casted$start_max[i],
      missing_casted$end_min[i], missing_casted$end_max[i]
    )

    regions_result <- calculate_regions(
      missing_casted$start_min[i], missing_casted$start_max[i],
      missing_casted$end_min[i], missing_casted$end_max[i],
      missing_casted$region_type[i]
    )

    missing_casted$certain_start[i] <- regions_result$certain_start
    missing_casted$certain_end[i] <- regions_result$certain_end
    missing_casted$uncertain_left_start[i] <- regions_result$uncertain_left_start
    missing_casted$uncertain_left_end[i] <- regions_result$uncertain_left_end
    missing_casted$uncertain_right_start[i] <- regions_result$uncertain_right_start
    missing_casted$uncertain_right_end[i] <- regions_result$uncertain_right_end
    missing_casted$start_is_single[i] <- regions_result$start_is_single
    missing_casted$end_is_single[i] <- regions_result$end_is_single
  }

  # Add seqid
  missing_casted$seqid <- "NZ_CP128453.1"

  message("Missing data processing completed!")
} # missing_df 처리 종료

###############################################################################################################################################
# Function to parse mutation_data
parse_mutation_data <- function(data_frame, genome_id, chr_length) {
  if (nrow(data_frame) == 0) {
    return(data.frame(
      position = numeric(), ALE_version = character(), mutation = character(),
      annotation = character(), gene = character(), description = character(),
      label = character(), new = logical(), start = numeric(), end = numeric(),
      type = character(), seqid = character(), stringsAsFactors = FALSE
    ))
  }

  parsed_data <- dplyr::mutate(data_frame,
                               start = position,
                               end = NA_real_,
                               type = NA_character_,
                               seqid = NA_character_)

  if (all(is.na(parsed_data$seqid))) {
    if (nrow(chr_length) == 1) {
      parsed_data$seqid <- chr_length$seqid[1]
    } else {
      parsed_data$seqid <- chr_length$seqid[1]
    }
  }

  for (i in 1:nrow(parsed_data)) {
    mut <- as.character(parsed_data$mutation[i])
    pos <- as.numeric(parsed_data$start[i])

    if (is.na(mut) || is.na(pos)) {
      next
    }

    if (stringr::str_detect(mut, "^\\+[ATCG]+$")) {
      parsed_data$end[i] <- pos + 1
      parsed_data$type[i] <- "INS"
    } else if (stringr::str_detect(mut, "^\\d+\\s*bp→[ATGC]+$")) {
      num_bp <- as.numeric(stringr::str_extract(mut, "\\d+"))
      parsed_data$end[i] <- pos + num_bp - 1
      parsed_data$type[i] <- "REP"
    } else if (stringr::str_detect(mut, "^[ATGC]→[ATGC]$")) {
      parsed_data$end[i] <- pos
      parsed_data$type[i] <- "SNP"
    } else if (stringr::str_detect(mut, "^[∆Δ][0-9,]+(?:\\s*bp)?$")) {
      num_deleted_str <- stringr::str_extract(mut, "[0-9,]+")
      num_deleted <- as.numeric(gsub(",", "", num_deleted_str))

      if (!is.na(num_deleted)) {
        parsed_data$end[i] <- pos + num_deleted - 1
        parsed_data$type[i] <- "DEL"
      } else {
        parsed_data$end[i] <- NA
        parsed_data$type[i] <- "UNKNOWN"
      }
    } else {
      parsed_data$end[i] <- NA
      parsed_data$type[i] <- "UNKNOWN"
    }
  }
  return(parsed_data)
}

###############################################################################################################################################
# Enhanced breseqVisualization function with missing data integration
breseqVisualization <- function(fna_count = NULL, completeness = NULL, mutation_data = NULL, missing_data = NULL,
                                hide_intergenic_labels = HIDE_INTERGENIC_LABELS,
                                gc_colors = GC_COLORS) {
  gc_palette <- .breseq_one_shot_gc_colors(gc_colors)
  # Load required packages
  library(circlize)
  library(dplyr)
  library(ComplexHeatmap)
  library(stringr)

  # Function to get screen size
  get_screen_size <- function(default_width = 12, default_height = 7, dpi = 600) {
    screen_width_in <- default_width
    screen_height_in <- default_height
    tryCatch({
      display_info <- system("system_profiler SPDisplaysDataType | grep Resolution", intern = TRUE)
      if (length(display_info) > 0) {
        resolution <- strsplit(display_info[1], "x")[[1]]
        width_px <- as.numeric(gsub("[^0-9]", "", resolution[1]))
        height_px <- as.numeric(gsub("[^0-9]", "", resolution[2]))
        screen_width_in <- width_px / dpi
        screen_height_in <- height_px / dpi
        message("Detected screen resolution: ", width_px, "x", height_px, " pixels")
        message("Converted to inches (DPI=", dpi, "): ", round(screen_width_in, 2), "x", round(screen_height_in, 2))
        screen_width_in <- min(screen_width_in, 15)
        screen_height_in <- min(screen_height_in, 15 * (height_px / width_px))
      }
    }, error = function(e) {
      message("Failed to detect screen size: ", e$message)
      message("Using default PDF size: ", default_width, "x", default_height, " inches")
    })
    return(list(width = screen_width_in, height = screen_height_in))
  }

  # Set PDF dimensions
  if (!exists("pdf_width") || !exists("pdf_height") ||
      is.null(pdf_width) || is.null(pdf_height)) {
    screen_size <- get_screen_size()
    pdf_width <- screen_size$width*2
    pdf_height <- screen_size$height*2
  }
  message("Using PDF dimensions: ", pdf_width, "x", pdf_height, " inches")

  # Debug: Check if missing_data is provided
  message("Checking missing_data at function entry: ", ifelse(is.null(missing_data), "NULL", "Provided"))
  if (!is.null(missing_data)) {
    message("Missing data structure at entry:")
    print(str(missing_data))
  }

  # Debug: Check if mutation_data is provided
  message("Checking mutation_data at function entry: ", ifelse(is.null(mutation_data), "NULL", "Provided"))
  if (!is.null(mutation_data)) {
    message("Mutation data structure at entry:")
    print(str(mutation_data))
  }

  # Set fna_count to the number of rows in contig_length if it is NULL
  if (is.null(fna_count)) {
    if (!exists("contig_length", envir = .GlobalEnv)) {
      stop("Error: contig_length not found. Define contig_length or specify fna_count.")
    }
    fna_count <- nrow(contig_length)
  }

  # Check gff existence
  if (!exists("gff", envir = .GlobalEnv)) {
    stop("Error: gff object not found. Define gff data.")
  }

  for (i in 1:fna_count) {
    seqid <- paste((gff %>% filter(source == "Local" | type == "region") %>% arrange(desc(end)) %>% select(seqid))[i, ] %>% as.character())

    # Initialize factors for the chromosome i
    if (!exists(paste(seqid, "length", sep = "_"), envir = .GlobalEnv)) {
      message("Error: ", seqid, "_length not found. Skipping contig.")
      dev.off()
      next
    }

    # Clear previous plots and set up circos parameters
    circos.clear()
    col_text <- "grey40"
    circos.par("start.degree" = 90, "track.height" = 0.8, gap.degree = 0, cell.padding = c(0, 0, 0, 0))

    chr_length <- get(paste(seqid, "length", sep = "_"), envir = .GlobalEnv)
    message("Processing contig: ", seqid, ", length: start = ", chr_length$start, ", end = ", chr_length$end)

    # Parse and color mutation_data for each sample in sample_mutation_data_list
    mutation_data_list_parsed <- list()
    if (exists("mutation_list", envir = .GlobalEnv) && length(mutation_list) > 0) {
      for (sample_name in names(mutation_list)) {
        data_for_parsing <- mutation_list[[sample_name]]
        if (nrow(data_for_parsing) > 0) {
          parsed_sample_data <- parse_mutation_data(data_for_parsing, Genome_summary$VERSION[1], chr_length)

          # Add color column based on mutation type
          parsed_sample_data <- parsed_sample_data %>%
            dplyr::mutate(
              type_color = dplyr::case_when(
                type == "INS" ~ "#8C1A18",
                type == "DEL" ~ "#1D448F",
                type == "SNP" ~ "#24572A",
                type == "REP" ~ "#E5A544",
                TRUE ~ "gray"
              )
            )
          mutation_data_list_parsed[[sample_name]] <- parsed_sample_data
        } else {
          mutation_data_list_parsed[[sample_name]] <- data.frame()
        }
      }

      message("Parsed mutation_data_list_parsed (first sample head, if available):")
      if (length(mutation_data_list_parsed) > 0 && nrow(mutation_data_list_parsed[[1]]) > 0) {
        print(head(mutation_data_list_parsed[[1]]))
      } else {
        message("No parsed data in mutation_data_list_parsed or first sample is empty.")
      }

      # Aggregate all parsed mutations for the general points track and global legend
      mutation_data <- dplyr::bind_rows(mutation_data_list_parsed) %>% as.data.frame()
      cat("Total number of total identified SNVs across all ALE versions:", nrow(mutation_data), "\n")
      mutation_new_data <- mutation_data %>% dplyr::filter(new==TRUE)
      cat("Total number of unique SNVs identified across all ALE versions:", nrow(mutation_new_data), "\n")
    } else {
      message("No mutation data to parse for individual samples.")
    }

    # Integrate missing regions into existing mutation_new_data
    if (!is.null(missing_data) && nrow(missing_data) > 0 && exists("mutation_new_data") && !is.null(mutation_new_data)) {
      message("Integrating missing regions into mutation_new_data")

      # Filter missing data for current seqid
      if ("seqid" %in% colnames(missing_data)) {
        seqid_missing <- missing_data %>% dplyr::filter(trimws(.data$seqid) == trimws(seqid))
      } else {
        seqid_missing <- missing_data %>% dplyr::mutate(seqid = seqid)
      }

      if (nrow(seqid_missing) > 0) {
        # Create missing rectangles
        missing_rects_list <- list()
        region_colors <- list("MC_inversion" = "#FF6B6B", "MC_repeat" = "#4ECDC4")

        for (j in 1:nrow(seqid_missing)) {
          row <- seqid_missing[j, ]
          base_color <- region_colors[[row$region_type]]
          if (is.null(base_color)) base_color <- "#808080"

          # Certain region
          if (!is.na(row$certain_start) && !is.na(row$certain_end)) {
            certain_rect <- data.frame(
              seqid = seqid, start = row$certain_start, end = row$certain_end,
              type = "certain", region_type = row$region_type, alpha = 1.0,
              color = base_color, ALE_version = row$ALE_version,
              description = as.character(row$description), stringsAsFactors = FALSE
            )
            missing_rects_list <- append(missing_rects_list, list(certain_rect))
          }

          # Uncertain left region
          if (!is.na(row$uncertain_left_start) && !is.na(row$uncertain_left_end)) {
            uncertain_left_rect <- data.frame(
              seqid = seqid, start = row$uncertain_left_start, end = row$uncertain_left_end,
              type = "uncertain_left", region_type = row$region_type, alpha = 0.5,
              color = base_color, ALE_version = row$ALE_version,
              description = as.character(row$description), stringsAsFactors = FALSE
            )
            missing_rects_list <- append(missing_rects_list, list(uncertain_left_rect))
          }

          # Uncertain right region
          if (!is.na(row$uncertain_right_start) && !is.na(row$uncertain_right_end)) {
            uncertain_right_rect <- data.frame(
              seqid = seqid, start = row$uncertain_right_start, end = row$uncertain_right_end,
              type = "uncertain_right", region_type = row$region_type, alpha = 0.5,
              color = base_color, ALE_version = row$ALE_version,
              description = as.character(row$description), stringsAsFactors = FALSE
            )
            missing_rects_list <- append(missing_rects_list, list(uncertain_right_rect))
          }
        }

        # Convert to dataframe and filter
        if (length(missing_rects_list) > 0) {
          missing_rects <- do.call(rbind, missing_rects_list)
          missing_rects <- missing_rects %>%
            dplyr::filter(start >= chr_length$start & end <= chr_length$end & start <= end)

          # Ensure version_colors exists
          if (!exists("version_colors", envir = .GlobalEnv)) {
            ALE_versions <- unique(c(mutation_new_data$ALE_version, missing_rects$ALE_version))
            ALE_versions <- ALE_versions[gtools::mixedorder(ALE_versions)]
            n_versions <- length(ALE_versions)
            version_colors <- setNames(
              RColorBrewer::brewer.pal(min(n_versions, 8), "Set2")[1:n_versions],
              ALE_versions
            )
            assign("version_colors", version_colors, envir = .GlobalEnv)
          }

          # Convert missing regions to mutation format
          missing_as_mutations <- missing_rects %>%
            dplyr::mutate(
              position = (start + end) / 2,
              mutation = paste0("MISSING_", region_type, "_", type),
              annotation = description,
              gene = description,
              label = ifelse(type == "certain", paste0(region_type, " (certain)"), paste0(region_type, " (uncertain)")),
              new = TRUE,
              type_color = color,
              alpha_val = alpha,
              version_color = version_colors[as.character(ALE_version)]
            ) %>%
            dplyr::select(position, ALE_version, mutation, annotation, gene, description, label, new,
                          version_color, start, end, type, seqid, type_color, color, alpha_val)

          # Add alpha_val to mutation data
          if (!"alpha_val" %in% colnames(mutation_new_data)) {
            mutation_new_data$alpha_val <- 1.0
          }

          # Ensure all required columns exist
          required_cols <- c("position", "ALE_version", "mutation", "annotation", "gene", "description",
                             "label", "new", "version_color", "start", "end", "type", "seqid", "type_color", "color")

          for (col in required_cols) {
            if (!col %in% colnames(mutation_new_data)) {
              mutation_new_data[[col]] <- NA
            }
          }

          # Combine data
          mutation_new_data <- dplyr::bind_rows(mutation_new_data, missing_as_mutations)
          message("Enhanced mutation_new_data with missing regions: ", nrow(mutation_new_data), " total rows")
        }
      }
    }

    # Dynamically adjust xlim to cover all data points
    max_end <- chr_length$end

    # recording plot
    pdf_file <- paste0(seqid, "_map.pdf")
    grDevices::pdf(file = pdf_file, width = pdf_width*2, height = pdf_height*2)

    circos.initialize(factors = c(seqid), xlim = matrix(c(chr_length$start, max_end), ncol = 2))

    # Calculate genome length and determine step size and label units
    genome_length <- chr_length$end
    if (exists("nucleotide", envir = .GlobalEnv) && seqid %in% names(nucleotide)) {
      genome_length <- nchar(get(names(nucleotide)[names(nucleotide) == seqid], envir = .GlobalEnv))
      message("Using nucleotide sequence length: ", genome_length)
    } else {
      message("Warning: nucleotide object missing or does not contain ", seqid, ". Using chr_length$end: ", genome_length)
    }
    step_size <- genome_length / 12
    step_size <- ceiling(step_size / (10^(nchar(step_size) - 1))) * 10^(nchar(step_size) - 1)

    # Determine label unit (Mb for large genomes, kb for smaller genomes)
    if (genome_length >= 1e6) {
      label_unit <- "Mb"
      scale_factor <- 1e6
    } else {
      label_unit <- "kb"
      scale_factor <- 1e3
    }

    # Create dynamic break points for axis labeling
    brk <- seq(0, genome_length, by = step_size)

    # Process mutation data
    mut_data <- NULL
    if (!is.null(mutation_new_data)) {
      message("Mutation data provided. Structure:")
      print(str(mutation_new_data))

      # Debug: Check ALE_version distribution
      message("ALE_version distribution in mutation_new_data: ", paste(table(mutation_new_data$ALE_version), collapse = ", "))

      # Check if raw data (position, mutation) or pre-parsed (start, end, type)
      if (all(c("position", "mutation") %in% colnames(mutation_new_data))) {
        message("Detected raw mutation data. Parsing...")
        mut_data <- parse_mutation_data(mutation_new_data, seqid, chr_length)
        # Ensure version_color and color are preserved
        if (!"version_color" %in% colnames(mut_data)) {
          mut_data$version_color <- mutation_new_data$version_color[match(mut_data$position, mutation_new_data$position)]
        }
        if (!"color" %in% colnames(mut_data)) {
          mut_data$color <- "gray"
        }
      } else {
        # Validate pre-parsed data
        required_cols <- c("start", "end", "type", "color", "ALE_version")
        if (!all(required_cols %in% colnames(mutation_new_data))) {
          message(paste("Skipping mutation tracks for", seqid, "due to missing required columns in mutation_data:", paste(required_cols[!required_cols %in% colnames(mutation_new_data)], collapse = ", ")))
        } else {
          message("Unique seqid values in mutation_data: ", paste(unique(mutation_new_data$seqid), collapse = ", "))
          mut_data <- if ("seqid" %in% colnames(mutation_new_data)) {
            mutation_new_data %>% dplyr::filter(trimws(.data$seqid) == trimws(seqid))
          } else {
            mutation_new_data %>% dplyr::mutate(seqid = seqid)
          }
          message("Rows before filtering: ", nrow(mut_data))
          # Filter coordinates
          mut_data <- mut_data %>% dplyr::filter(start >= chr_length$start & end <= chr_length$end & start <= end)
          message("Rows after coordinate filtering: ", nrow(mut_data))
        }
      }

      # Debug: Check mut_data contents
      if (!is.null(mut_data) && nrow(mut_data) > 0) {
        message("mut_data columns: ", paste(colnames(mut_data), collapse = ", "))
        message("mut_data ALE_version distribution: ", paste(table(mut_data$ALE_version), collapse = ", "))
      }

      # Add labels for mutations with non-empty label column
      if (!is.null(mut_data) && "label" %in% colnames(mut_data) && any(!is.na(mut_data$label) & mut_data$label != "")) {
        label_data <- mut_data %>%
          dplyr::filter(!is.na(label) & label != "" & !is.na(start) & !is.na(end) & !is.na(version_color)) %>%
          dplyr::filter(start <= end)
        label_data <- .breseq_one_shot_filter_intergenic_labels(label_data, hide_intergenic_labels) %>%
          dplyr::select(seqid, start, end, label, version_color)

        if (nrow(label_data) > 0) {
          message("Adding labels for ", nrow(label_data), " mutations")

          # Additional validation
          label_data <- label_data %>%
            dplyr::filter(
              !is.na(start) & !is.na(end) &
                start > 0 & end > 0 &
                start <= end &
                !is.na(version_color) & version_color != ""
            )

          if (nrow(label_data) > 0) {
            tryCatch({
              circos.genomicLabels(
                bed = label_data,
                labels.column = 4,
                side = "outside",
                col = label_data$version_color,
                line_col = label_data$version_color,
                cex = 0.6,
                connection_height = mm_h(5),
                line_lwd = 0.5
              )
            }, error = function(e) {
              message("Error adding labels: ", e$message)
              message("Skipping label addition...")
            })
          }
        }
      }
    } else {
      message(paste("Skipping mutation tracks for", seqid, "due to missing mutation_data."))
    }

    # Plot cumulative mutation tracks and labels (INTEGRATED with missing regions)
    if (!is.null(mut_data) && nrow(mut_data) > 0 && "ALE_version" %in% colnames(mut_data)) {
      message("Using mut_data for tracks (includes missing regions)")
      mutation_data <- mut_data

      message("Mutation data before processing: ", nrow(mutation_data), " rows")
      message("mutation_data columns: ", paste(colnames(mutation_data), collapse = ", "))
      message("mutation_data ALE_version distribution: ", paste(table(mutation_data$ALE_version), collapse = ", "))

      # Plot cumulative tracks for each ALE_version
      groups <- samples
      message("Mutation groups (ALE_version): ", paste(groups, collapse = ", "))
      track_index <- get.current.track.index() + 1

      # Set transparency (alpha: 0 = fully transparent, 1 = fully opaque)
      alpha_bg <- 0.3    # Transparency for track background

      for (i in seq_along(groups)) {
        ALE_group <- groups[i]
        group_data <- mutation_data %>% dplyr::filter(ALE_version %in% groups[1:i])
        bg_color <- version_colors[as.character(ALE_group)]

        message("Plotting cumulative track for ALE_version <= ", ALE_group, ", rows: ", nrow(group_data))

        if (nrow(group_data) > 0) {
          # Apply transparency based on alpha_val if available
          if ("alpha_val" %in% colnames(group_data)) {
            group_data$version_color_trans <- mapply(function(col, alpha_val) {
              adjustcolor(col, alpha.f = alpha_val)
            }, group_data$version_color, group_data$alpha_val)
          } else {
            group_data$version_color_trans <- sapply(group_data$version_color, function(col) {
              adjustcolor(col, alpha.f = 1.0)
            })
          }

          # Apply transparency to background color
          bg_color_trans <- adjustcolor(bg_color, alpha.f = alpha_bg)

          circos.track(
            factors = as.character(group_data$seqid),
            panel.fun = function(region, value, ...) {
              for (j in 1:nrow(group_data)) {
                if ("alpha_val" %in% colnames(group_data) && grepl("MISSING_", group_data$type[j])) {
                  # Missing regions with border
                  circos.rect(
                    xleft = as.numeric(group_data$start[j]),
                    xright = as.numeric(group_data$end[j]),
                    ybottom = 0,
                    ytop = 1,
                    col = group_data$version_color_trans[j],
                    border = "black",
                    lwd = 0.3
                  )
                } else {
                  # Regular mutations
                  circos.rect(
                    xleft = as.numeric(group_data$start[j]),
                    xright = as.numeric(group_data$end[j]),
                    ybottom = 0,
                    ytop = 1,
                    col = group_data$version_color_trans[j],
                    border = NA
                  )
                }
              }
            },
            ylim = c(0, 1),
            track.index = track_index,
            bg.col = bg_color_trans,
            bg.border = FALSE,
            track.height = 0.03,
            track.margin = c(0.0, 0.0)
          )
          track_index <- track_index + 1
        } else {
          message("Skipping cumulative track for ALE_version <= ", ALE_group, " in ", seqid, " due to lack of data.")
        }
      }
    } else {
      message("Skipping mutation tracks for ", seqid, " due to no valid data in mutation_data.")
    }

    # Add genome x-axis with dynamic labels
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
      circos.axis(h = "top",
                  major.at = brk,
                  labels = paste(round(brk / scale_factor, 1), label_unit),
                  labels.cex = 0.8,
                  col = col_text,
                  labels.col = col_text,
                  lwd = 0.7,
                  labels.facing = "clockwise")
    }, bg.border = FALSE)

    # Plot genome base circle
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      xlim <- CELL_META$xlim
      ylim <- CELL_META$ylim
      circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = col_text, facing = "bending.inside", niceFacing = TRUE)
    }, bg.col = "grey70", bg.border = FALSE, track.height = 0.01)

    # Plot COG in 3' to 5' and 5' to 3' directions
    plot_cog <- function(seqid_suffix) {
      cog_data <- tryCatch(get(paste(seqid, seqid_suffix, sep = "_"), envir = .GlobalEnv), error = function(e) NULL)
      if (!is.null(cog_data) && nrow(cog_data) > 0) {
        if (nrow(cog_data) > 0) {
          circos.track(factors = cog_data$seqid,
                       panel.fun = function(region, value, ...) {
                         circos.genomicRect(
                           region = matrix(c(cog_data$start %>% as.numeric(),
                                             cog_data$end %>% as.numeric()), ncol = 2),
                           value = cog_data$COG_category_for_plot,
                           col = cog_data$COG_color,
                           border = NA
                         )
                       }, ylim = c(0, 1), track.index = get.current.track.index() + 1, bg.col = "grey90", bg.border = FALSE, track.height = 0.03)
        } else {
          message(paste("Skipping COG track for", seqid, "suffix", seqid_suffix, "due to no valid data after filtering."))
        }
      } else {
        message(paste("Skipping COG track for", seqid, "suffix", seqid_suffix, "due to lack of data."))
      }
    }

    plot_cog("3to5")
    plot_cog("5to3")

    # Plot RNA annotations
    if (!is.null(RNA) && any(RNA$seqid == seqid)) {
      RNA_data <- RNA %>% dplyr::filter(seqid == seqid)
      if (nrow(RNA_data) > 0) {
        # Filter coordinates
        RNA_data <- RNA_data %>% dplyr::filter(start >= chr_length$start & end <= chr_length$end & start <= end)
        if (nrow(RNA_data) > 0) {
          circos.track(factors = RNA_data$seqid %>% as.character(),
                       panel.fun = function(region, value, ...) {
                         circos.genomicRect(
                           region = matrix(c(RNA_data$start %>% as.numeric(), RNA_data$end %>% as.numeric()), ncol = 2),
                           value = RNA_data$type %>% as.character(),
                           col = RNA_data$color,
                           border = NA
                         )
                       }, ylim = c(0, 1), track.index = get.current.track.index() + 1, bg.col = "grey90", bg.border = FALSE, track.height = 0.03)
        } else {
          message(paste("Skipping RNA track for", seqid, "due to no valid data after filtering."))
        }
      } else {
        message(paste("Skipping RNA track for", seqid, "due to lack of data."))
      }
    }

    # Plot GC skew and GC ratio tracks if data exists
    plot_gc_track <- function(type) {
      gc_data <- get(paste("window", names(nucleotide)[i], sep = "_"), inherits = TRUE)
      if (!is.null(gc_data) && nrow(gc_data) > 0 && all(gc_data[[paste0(type, "_minus_average")]] != 0)) {
        if (any(!is.na(gc_data[[paste0(type, "_minus_average")]]))) {
          circlize::circos.track(factors = gc_data$seqid,
                                 panel.fun = function(x,y) {
                                   current_sector_data <- gc_data[gc_data$seqid == circlize::get.cell.meta.data("sector.index"),]
                                   if(nrow(current_sector_data) > 0){
                                     circlize::circos.barplot(
                                       value = current_sector_data[[paste0(type, "_minus_average")]] %>% as.numeric(),
                                       pos = (current_sector_data$start + current_sector_data$end) / 2,
                                       col = .breseq_one_shot_gc_track_colors(type,
                                         as.numeric(current_sector_data[[paste0(type, "_minus_average")]]), gc_palette),
                                       border = NA,
                                       bar_width = 10000
                                     )
                                   }
                                 }, ylim = c(-max(abs(gc_data[[paste0(type, "_minus_average")]] %>% as.numeric()), na.rm = TRUE),
                                             max(abs(gc_data[[paste0(type, "_minus_average")]] %>% as.numeric()), na.rm = TRUE)),
                                 track.index = circlize::get.current.track.index() + 1, bg.col = NA, bg.border = FALSE, track.height = 0.06)
        } else {
          message(paste("Skipping GC track for", seqid, "type", type, "due to lack of non-NA data."))
        }
      } else {
        message(paste("Skipping GC track for", seqid, "type", type, "due to lack of data or zero range."))
      }
    }

    # Check and plot GC skew and GC ratio within loop to prevent errors
    if (exists(paste("window", names(nucleotide)[i], sep = "_"))) {
      plot_gc_track("gc_skew")
      plot_gc_track("gc_ratio")
    } else {
      message(paste("Skipping GC skew and ratio tracks for", seqid, "due to missing data."))
    }

    # Add genome name and total base pairs at the center
    graphics::text(0, 0, paste(seqid,
                               paste(formatC(genome_length, format = "f", digits = 0, big.mark = ","), "bp"),
                               sep = "\n"), cex = 1.5)

    # Draw legends for COG, mutations, missing regions, and version_color
    lgd_list <- list()

    # COG legend
    if (exists("genbank_table", envir = .GlobalEnv)) {
      COG <- genbank_table %>% group_by(COG_legend) %>% slice_head() %>% select(COG_category_for_plot, COG_color, COG_legend) %>% na.omit()
      if (nrow(COG) > 0) {
        lgd_cog <- ComplexHeatmap::Legend(
          labels = COG$COG_legend,
          title = "COG category",
          type = "points",
          legend_gp = gpar(col = NA),
          ncol = 1,
          by_row = TRUE,
          direction = "vertical",
          background = COG$COG_color
        )
        lgd_list <- c(lgd_list, list(lgd_cog))
      }
    }

    # Mutation and Missing region legends (integrated)
    if (!is.null(mut_data) && nrow(mut_data) > 0) {
      # Separate mutation types from missing region types
      mut_legend <- mut_data %>%
        dplyr::filter(!grepl("MISSING_", type)) %>%
        group_by(type) %>%
        slice_head() %>%
        select(type, color) %>%
        na.omit()

      if (nrow(mut_legend) > 0) {
        lgd_mut <- ComplexHeatmap::Legend(
          labels = mut_legend$type,
          title = "Mutation Type",
          type = "points",
          legend_gp = gpar(col = NA),
          ncol = 1,
          by_row = TRUE,
          direction = "vertical",
          background = mut_legend$color
        )
        lgd_list <- c(lgd_list, list(lgd_mut))
      }

      # Missing region legend
      missing_legend <- mut_data %>%
        dplyr::filter(grepl("MISSING_", type)) %>%
        group_by(type) %>%
        slice_head() %>%
        select(type, color, alpha_val) %>%
        na.omit()

      if (nrow(missing_legend) > 0) {
        # Clean up type names for display
        display_labels <- gsub("MISSING_", "", missing_legend$type)
        display_labels <- gsub("_", " ", display_labels)

        # Apply transparency to colors for legend
        legend_colors <- mapply(function(color, alpha) {
          adjustcolor(color, alpha.f = alpha)
        }, missing_legend$color, missing_legend$alpha_val)

        lgd_missing <- ComplexHeatmap::Legend(
          labels = display_labels,
          title = "Missing Regions",
          type = "points",
          legend_gp = gpar(col = NA),
          ncol = 1,
          by_row = TRUE,
          direction = "vertical",
          background = legend_colors
        )
        lgd_list <- c(lgd_list, list(lgd_missing))
      }

      # Version_color legend (integrated data)
      version_legend <- mut_data %>% group_by(ALE_version) %>% slice_head() %>% select(ALE_version, version_color) %>% na.omit()
      if (nrow(version_legend) > 0) {
        lgd_version <- ComplexHeatmap::Legend(
          labels = version_legend$ALE_version,
          title = "ALE Version",
          type = "points",
          legend_gp = gpar(col = NA),
          ncol = 1,
          by_row = TRUE,
          direction = "vertical",
          background = version_legend$version_color
        )
        lgd_list <- c(lgd_list, list(lgd_version))
      }
    }

    # Pack and draw legends
    if (length(lgd_list) > 0) {
      pd = ComplexHeatmap::packLegend(list = lgd_list)
      ComplexHeatmap::draw(pd, x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
    }

    # Display plot and pause for inspection
    message("Displaying plot for ", seqid)
    dev.flush() # Ensure plot is rendered

    grDevices::dev.off()
    message("Saved PDF for ", seqid, " as ", pdf_file)
  }
  circos.clear()
}
