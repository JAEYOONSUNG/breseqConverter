#' Circos-style Genome Visualization
#'
#' Creates publication-quality circular genome maps using circlize package,
#' showing mutations with proper type-based coloring and cumulative ALE version tracks.
#'
#' @param mutation_data Data frame with mutation data (from pipeline or Excel)
#' @param missing_data Optional data frame with missing coverage regions
#' @param reference_file Path to reference GenBank file for genome info
#' @param output_file Path for output PDF file
#' @param samples Sample/stage names in order (auto-detected if NULL)
#' @param genome_length Optional genome length (auto-detected from reference if NULL)
#' @param seqid Sequence ID for the genome (auto-detected if NULL)
#' @return Path to the created PDF file
#' @export
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_circos_genome(
#'   mutation_data = mutations,
#'   reference_file = "reference.gbk",
#'   output_file = "genome_map.pdf"
#' )
#' }
plot_circos_genome <- function(mutation_data,
                               missing_data = NULL,
                               reference_file = NULL,
                               output_file = "genome_map.pdf",
                               samples = NULL,
                               genome_length = NULL,
                               seqid = NULL) {

  # Check for required packages
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' required. Install with: install.packages('circlize')")
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' required. Install with: BiocManager::install('ComplexHeatmap')")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' required. Install with: install.packages('stringr')")
  }

  library(circlize)
  library(dplyr)
  library(ComplexHeatmap)
  library(stringr)

  # Mutation type colors (사용자 지정 팔레트)
  type_colors <- c(
    "INS" = "#FF5C8D",   # Bright Pink - Insertion
    "DEL" = "#03045E",   # Dark Navy - Deletion
    "SNP" = "#139487",   # Teal - Single nucleotide polymorphism
    "REP" = "#E5C87B",   # Golden Yellow - Replacement
    "UNKNOWN" = "gray"   # Gray
  )

  # Parse reference genome info
  if (!is.null(reference_file) && file.exists(reference_file)) {
    ref_info <- parse_reference_info(reference_file)
    if (is.null(genome_length)) genome_length <- ref_info$length[1]
    if (is.null(seqid)) seqid <- ref_info$seq_id[1]
  } else {
    if (is.null(genome_length)) {
      genome_length <- max(as.numeric(gsub(",", "", mutation_data$position)), na.rm = TRUE) * 1.1
    }
    if (is.null(seqid)) seqid <- "Genome"
  }

  message("Creating Circos plot for ", seqid, " (", format(genome_length, big.mark = ","), " bp)")

  # Auto-detect samples from data
  if (is.null(samples)) {
    if ("ALE_version" %in% names(mutation_data)) {
      samples <- unique(mutation_data$ALE_version)
      if (requireNamespace("gtools", quietly = TRUE)) {
        samples <- samples[gtools::mixedorder(samples)]
      }
    } else {
      # Try to detect from column names
      sample_cols <- grep("\\.(45)?SJY\\d+|/output/index", names(mutation_data), value = TRUE, ignore.case = TRUE)
      if (length(sample_cols) > 0) {
        samples <- unique(na.omit(sapply(sample_cols, function(col) {
          # Extract sample name like "45SJY3" from column name
          m <- regmatches(col, regexpr("(45)?SJY\\d+", col, ignore.case = TRUE))
          if (length(m) > 0 && nchar(m) > 0) return(m)
          return(NA)
        })))
        if (requireNamespace("gtools", quietly = TRUE)) {
          samples <- samples[gtools::mixedorder(samples)]
        }
      }
    }
  }

  if (length(samples) == 0) {
    samples <- "Sample"
  }

  message("Samples (ALE versions): ", paste(samples, collapse = " -> "))

  # Sample/version colors for track backgrounds
  n_samples <- length(samples)
  sample_colors <- setNames(
    RColorBrewer::brewer.pal(max(n_samples, 3), "Set2")[1:n_samples],
    samples
  )

  # Parse and process mutation data
  mut_df <- parse_and_prepare_mutations(mutation_data, seqid, genome_length, samples, type_colors, sample_colors)

  if (nrow(mut_df) == 0) {
    message("Warning: No mutations to plot after parsing")
    return(NULL)
  }

  message("Total mutations parsed: ", nrow(mut_df))
  message("Mutation types: ", paste(table(mut_df$type), collapse = ", "))

  # Calculate cumulative "new" mutations per ALE version
  mut_df <- track_new_mutations(mut_df, samples)

  new_count <- sum(mut_df$new, na.rm = TRUE)
  message("New (unique) mutations: ", new_count)

  # Open PDF device
  pdf_width <- 14
  pdf_height <- 14
  grDevices::pdf(file = output_file, width = pdf_width, height = pdf_height)

  # Clear and initialize circos
  circos.clear()
  circos.par(
    "start.degree" = 90,
    "track.height" = 0.8,
    "gap.degree" = 0,
    "cell.padding" = c(0, 0, 0, 0)
  )

  # Initialize with genome coordinates
  circos.initialize(factors = seqid, xlim = c(1, genome_length))

  # Calculate axis labels
  step_size <- genome_length / 12
  step_size <- ceiling(step_size / (10^(nchar(as.character(as.integer(step_size))) - 1))) * 10^(nchar(as.character(as.integer(step_size))) - 1)
  brk <- seq(0, genome_length, by = step_size)

  if (genome_length >= 1e6) {
    label_unit <- "Mb"
    scale_factor <- 1e6
  } else {
    label_unit <- "kb"
    scale_factor <- 1e3
  }

  # Add outer labels for NEW mutations only
  if ("label" %in% names(mut_df)) {
    label_data <- mut_df %>%
      dplyr::filter(new == TRUE & !is.na(label) & label != "" & !is.na(start) & !is.na(end)) %>%
      dplyr::mutate(seqid = seqid) %>%
      dplyr::select(seqid, start, end, label, version_color) %>%
      dplyr::filter(start <= end)

    if (nrow(label_data) > 0) {
      message("Adding labels for ", nrow(label_data), " new mutations")
      tryCatch({
        circos.genomicLabels(
          bed = label_data,
          labels.column = 4,
          side = "outside",
          col = label_data$version_color,
          line_col = label_data$version_color,
          cex = 0.5,
          connection_height = mm_h(5),
          line_lwd = 0.5
        )
      }, error = function(e) {
        message("  Warning: Could not add labels: ", e$message)
      })
    }
  }

  # Create cumulative mutation tracks for each ALE version
  track_index <- get.current.track.index() + 1
  alpha_bg <- 0.3  # Background transparency

  for (i in seq_along(samples)) {
    sample_name <- samples[i]

    # Cumulative: show all NEW mutations from first sample up to current sample
    cumulative_samples <- samples[1:i]
    group_data <- mut_df %>%
      dplyr::filter(ALE_version %in% cumulative_samples & new == TRUE)

    bg_color <- adjustcolor(sample_colors[sample_name], alpha.f = alpha_bg)

    message("Track ", i, " (", sample_name, "): ", nrow(group_data), " cumulative new mutations")

    if (nrow(group_data) > 0) {
      group_data$seqid <- seqid

      tryCatch({
        circos.track(
          factors = seqid,
          ylim = c(0, 1),
          track.index = track_index,
          bg.col = bg_color,
          bg.border = FALSE,
          track.height = 0.04,
          track.margin = c(0.002, 0.002),
          panel.fun = function(x, y) {
            for (j in 1:nrow(group_data)) {
              # Color by mutation TYPE (not by version)
              circos.rect(
                xleft = group_data$start[j],
                xright = group_data$end[j],
                ybottom = 0.1,
                ytop = 0.9,
                col = group_data$type_color[j],
                border = NA
              )
            }
          }
        )
        track_index <- track_index + 1
      }, error = function(e) {
        message("  Warning: Could not create track for ", sample_name, ": ", e$message)
      })
    } else {
      # Empty track with just background
      tryCatch({
        circos.track(
          factors = seqid,
          ylim = c(0, 1),
          track.index = track_index,
          bg.col = bg_color,
          bg.border = FALSE,
          track.height = 0.04,
          track.margin = c(0.002, 0.002)
        )
        track_index <- track_index + 1
      }, error = function(e) {
        message("  Warning: Could not create empty track for ", sample_name, ": ", e$message)
      })
    }
  }

  # Axis track
  circos.track(
    track.index = get.current.track.index(),
    panel.fun = function(x, y) {
      circos.axis(
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
  circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      chr <- CELL_META$sector.index
      xlim <- CELL_META$xlim
      ylim <- CELL_META$ylim
      circos.text(mean(xlim), mean(ylim), chr, cex = 0.5, col = "grey40",
                  facing = "bending.inside", niceFacing = TRUE)
    },
    bg.col = "grey70",
    bg.border = FALSE,
    track.height = 0.01
  )

  # Center text
  graphics::text(0, 0, paste(seqid,
                             paste(formatC(genome_length, format = "f", digits = 0, big.mark = ","), "bp"),
                             sep = "\n"), cex = 1.2)

  # Build legends
  lgd_list <- list()

  # 1. Mutation type legend (by color)
  used_types <- unique(mut_df$type[mut_df$new == TRUE])
  used_types <- used_types[!is.na(used_types)]

  if (length(used_types) > 0) {
    type_legend_colors <- type_colors[used_types]
    lgd_mut <- ComplexHeatmap::Legend(
      labels = names(type_legend_colors),
      title = "Mutation Type",
      type = "points",
      legend_gp = grid::gpar(col = NA),
      ncol = 1,
      by_row = TRUE,
      direction = "vertical",
      background = unname(type_legend_colors)
    )
    lgd_list <- c(lgd_list, list(lgd_mut))
  }

  # 2. ALE version legend (track background colors)
  lgd_sample <- ComplexHeatmap::Legend(
    labels = samples,
    title = "ALE Version (Track)",
    type = "points",
    legend_gp = grid::gpar(col = NA),
    ncol = 1,
    by_row = TRUE,
    direction = "vertical",
    background = adjustcolor(sample_colors, alpha.f = alpha_bg)
  )
  lgd_list <- c(lgd_list, list(lgd_sample))

  # Draw legends
  if (length(lgd_list) > 0) {
    pd <- ComplexHeatmap::packLegend(list = lgd_list)
    ComplexHeatmap::draw(pd, x = grid::unit(4, "mm"), y = grid::unit(4, "mm"),
                         just = c("left", "bottom"))
  }

  grDevices::dev.off()
  circos.clear()

  message("Saved Circos plot: ", output_file)
  return(normalizePath(output_file))
}


#' Parse and Prepare Mutation Data for Circos
#'
#' Parses mutation strings to determine type and calculates proper end positions.
#'
#' @param mutation_data Raw mutation data frame
#' @param seqid Sequence ID
#' @param genome_length Genome length
#' @param samples Sample names in order
#' @param type_colors Named vector of colors by type
#' @param sample_colors Named vector of colors by sample
#' @return Parsed data frame ready for plotting
#' @keywords internal
parse_and_prepare_mutations <- function(mutation_data, seqid, genome_length, samples, type_colors, sample_colors) {

  df <- mutation_data

  # If already in long format with ALE_version, use directly
  if ("ALE_version" %in% names(df) && "position" %in% names(df) && "mutation" %in% names(df)) {
    message("Data is already in long format")

    # Ensure position is numeric
    df$position <- as.numeric(gsub(",", "", as.character(df$position)))
    df$start <- df$position

    # Parse mutation type and calculate end position
    for (i in 1:nrow(df)) {
      result <- parse_mutation_string(df$mutation[i], df$position[i])
      df$end[i] <- result$end
      df$type[i] <- result$type
    }

  } else {
    # Convert from wide format to long format
    message("Converting from wide to long format")

    # Find position column
    if (!"position" %in% names(df)) {
      stop("Data must have a 'position' column")
    }

    df$position <- as.numeric(gsub(",", "", as.character(df$position)))

    # Find mutation columns (e.g., mutation.45SJY3/output/index)
    mut_cols <- grep("^mutation\\.", names(df), value = TRUE, ignore.case = TRUE)
    ann_cols <- grep("^annotation\\.", names(df), value = TRUE, ignore.case = TRUE)
    gene_cols <- grep("^gene\\.", names(df), value = TRUE, ignore.case = TRUE)

    if (length(mut_cols) == 0) {
      # Try simple column names
      if ("mutation" %in% names(df)) {
        mut_cols <- "mutation"
      } else {
        stop("No mutation columns found")
      }
    }

    # Extract sample names from column names
    extract_sample <- function(col) {
      # Pattern 1: 45SJY followed by numbers
      m <- regmatches(col, regexpr("45SJY\\d+", col, ignore.case = TRUE))
      if (length(m) > 0 && nchar(m) > 0) return(m)

      # Pattern 2: SJY followed by numbers
      m <- regmatches(col, regexpr("SJY\\d+", col, ignore.case = TRUE))
      if (length(m) > 0 && nchar(m) > 0) return(m)

      # Fallback: extract between . and /
      parts <- strsplit(col, "\\.")[[1]]
      if (length(parts) >= 2) {
        tag <- sub("/.*", "", parts[2])
        if (nchar(tag) > 0) return(tag)
      }
      return(NA)
    }

    message("  Mutation columns found: ", paste(mut_cols, collapse = ", "))

    # Build long format data
    long_data <- data.frame()

    for (mut_col in mut_cols) {
      sample_name <- extract_sample(mut_col)
      if (is.na(sample_name)) next

      # Find corresponding annotation and gene columns
      ann_col <- ann_cols[grepl(sample_name, ann_cols, ignore.case = TRUE)]
      gene_col <- gene_cols[grepl(sample_name, gene_cols, ignore.case = TRUE)]

      for (i in 1:nrow(df)) {
        mut_val <- df[[mut_col]][i]
        if (is.na(mut_val) || mut_val == "") next

        pos <- df$position[i]
        ann_val <- if (length(ann_col) > 0) df[[ann_col[1]]][i] else NA
        gene_val <- if (length(gene_col) > 0) df[[gene_col[1]]][i] else NA

        # Parse mutation
        result <- parse_mutation_string(mut_val, pos)

        row_data <- data.frame(
          position = pos,
          start = pos,
          end = result$end,
          mutation = as.character(mut_val),
          annotation = as.character(ann_val),
          gene = as.character(gene_val),
          type = result$type,
          ALE_version = sample_name,
          stringsAsFactors = FALSE
        )

        long_data <- rbind(long_data, row_data)
      }
    }

    df <- long_data
  }

  # Remove rows with missing essential data
  df <- df[!is.na(df$position) & !is.na(df$start) & !is.na(df$end), ]

  # Add seqid
  df$seqid <- seqid

  # Filter by genome bounds
  df <- df[df$start >= 1 & df$end <= genome_length & df$start <= df$end, ]

  # Add type colors
  df$type_color <- sapply(df$type, function(t) {
    if (t %in% names(type_colors)) return(type_colors[t])
    return(type_colors["UNKNOWN"])
  })

  # Add version colors
  df$version_color <- sapply(df$ALE_version, function(v) {
    if (v %in% names(sample_colors)) return(sample_colors[v])
    return("#808080")
  })

  # Create label
  df$label <- paste(df$mutation, df$annotation, sep = " | ")
  df$label <- gsub("^NA \\| | \\| NA$|^NA \\| NA$", "", df$label)
  df$label <- trimws(df$label)
  df$label[df$label == "|" | df$label == "" | df$label == "NA"] <- NA

  return(df)
}


#' Parse Mutation String
#'
#' Parses a mutation string to determine type and end position.
#' Uses the same logic as the original one_shot.R parse_mutation_data().
#'
#' @param mut_str Mutation string (e.g., "+A", "A→G", "Δ100 bp")
#' @param position Start position
#' @return List with 'end' position and 'type'
#' @keywords internal
parse_mutation_string <- function(mut_str, position) {

  result <- list(end = position, type = "UNKNOWN")

  if (is.na(mut_str) || is.null(mut_str) || nchar(as.character(mut_str)) == 0) {
    return(result)
  }

  mut <- as.character(mut_str)
  pos <- as.numeric(position)

  if (is.na(pos)) {
    return(result)
  }

  # Insertion: +ATCG (adds bases after position)
  if (stringr::str_detect(mut, "^\\+[ATCG]+$")) {
    result$end <- pos + nchar(gsub("\\+", "", mut))
    result$type <- "INS"
    return(result)
  }

  # Homopolymer repeat change: (T)6→7, (A)8→9, etc.
  if (stringr::str_detect(mut, "^\\([ATCG]+\\)\\d+→\\d+$")) {
    # Extract the numbers before and after →
    nums <- as.numeric(stringr::str_extract_all(mut, "\\d+")[[1]])
    if (length(nums) >= 2) {
      old_len <- nums[1]
      new_len <- nums[2]
      if (new_len > old_len) {
        result$type <- "INS"  # Expansion
        result$end <- pos + new_len
      } else {
        result$type <- "DEL"  # Contraction
        result$end <- pos + old_len
      }
      return(result)
    }
  }

  # Replacement: 10 bp→ATCG (replaces N bp with new sequence)
  if (stringr::str_detect(mut, "^\\d+\\s*bp.*→[ATCG]+$")) {
    num_bp <- as.numeric(stringr::str_extract(mut, "\\d+"))
    result$end <- pos + num_bp - 1
    result$type <- "REP"
    return(result)
  }

  # SNP: A→G (single base change)
  if (stringr::str_detect(mut, "^[ATCG]→[ATCG]$")) {
    result$end <- pos
    result$type <- "SNP"
    return(result)
  }

  # Deletion: Δ100 bp or ∆100 bp or Δ2,340 bp (deletes N bases)
  if (stringr::str_detect(mut, "^[∆Δ][0-9,]+")) {
    num_deleted_str <- stringr::str_extract(mut, "[0-9,]+")
    num_deleted <- as.numeric(gsub(",", "", num_deleted_str))

    if (!is.na(num_deleted)) {
      result$end <- pos + num_deleted - 1
      result$type <- "DEL"
      return(result)
    }
  }

  # IS insertion or mobile element (typically shown with specific notation)
  if (stringr::str_detect(mut, "IS\\d+") || stringr::str_detect(mut, "\\+IS")) {
    result$end <- pos + 1000  # Approximate size for IS element
    result$type <- "INS"
    return(result)
  }

  # Default: point mutation
  result$end <- pos
  result$type <- "UNKNOWN"
  return(result)
}


#' Track New Mutations Across ALE Versions
#'
#' Identifies which mutations are "new" (first appearing) at each ALE version.
#'
#' @param df Data frame with mutations
#' @param samples Ordered sample names
#' @return Data frame with 'new' column added
#' @keywords internal
track_new_mutations <- function(df, samples) {

  if (nrow(df) == 0) return(df)

  # Sort samples by order
  if (requireNamespace("gtools", quietly = TRUE)) {
    samples <- samples[gtools::mixedorder(samples)]
  }

  # Track seen mutations (position:mutation combination)
  seen_mutations <- character()
  df$new <- FALSE

  for (sample in samples) {
    sample_rows <- which(df$ALE_version == sample)

    for (i in sample_rows) {
      # Create unique identifier for this mutation
      mut_id <- paste(df$position[i], df$mutation[i], sep = ":")

      if (!(mut_id %in% seen_mutations)) {
        df$new[i] <- TRUE
        seen_mutations <- c(seen_mutations, mut_id)
      }
    }
  }

  return(df)
}


#' Parse Reference GenBank File for Genome Info
#'
#' Extracts sequence ID and length from GenBank file.
#'
#' @param reference_file Path to GenBank file
#' @return Data frame with seq_id and length columns
#' @keywords internal
parse_reference_info <- function(reference_file) {
  lines <- readLines(reference_file, n = 100)

  # Find LOCUS line
  locus_line <- grep("^LOCUS", lines, value = TRUE)

  if (length(locus_line) == 0) {
    return(data.frame(seq_id = "Unknown", length = 1000000))
  }

  # Parse LOCUS line: LOCUS       NZ_CP128453.1         4628477 bp    DNA     circular BCT 05-APR-2024
  parts <- strsplit(locus_line[1], "\\s+")[[1]]

  seq_id <- parts[2]
  length_idx <- which(grepl("^\\d+$", parts))
  if (length(length_idx) > 0) {
    genome_length <- as.numeric(parts[length_idx[1]])
  } else {
    genome_length <- 1000000
  }

  return(data.frame(seq_id = seq_id, length = genome_length))
}
