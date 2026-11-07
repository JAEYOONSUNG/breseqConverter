# breseqConverter: one file, automatic extraction and plotting (2026-11-07)
# RStudio: edit the settings below, then Source. Terminal: Rscript --vanilla one_shot.R
# No breseqConverter installation, external R scripts, or GenomeDrawer globals required.
# Required R libraries are checked before reading data. Original inputs are read only.
# Set INPUT_DIR to your data folder. New results are saved outside the inputs.
INPUT_DIR <- getwd()  # Or set an absolute path to your breseq input folder
OUTPUT_DIR <- file.path(path.expand('~/Desktop'), paste0('breseqConverter_run_', format(Sys.time(), '%Y%m%d_%H%M%S')))
REFERENCE_FILE <- NULL  # NULL: find the reference GenBank/FASTA automatically
EXCEL_FILE <- NULL      # NULL: read breseq HTML; or supply a canonical exported workbook
GFF_FILE <- NULL        # NULL: auto-detect gene annotations
EGGNOG_FILE <- NULL     # NULL: auto-detect COG annotations
HIDE_INTERGENIC_LABELS <- FALSE  # TRUE hides labels only; marks remain
SAMPLE_ORDER <- NULL  # Or c('WT', 'sample1', 'sample2') in chronological order
SAMPLE_LABELS <- NULL # Optional named map, e.g. c(WT='Parent', sample1='Strain 1')
BASELINE_ID <- 'WT'  # Empty reference baseline when WT was not resequenced
SAMPLE_ALIASES <- NULL # Optional map from audited-evidence IDs to input sample IDs
DISPLAY_DELETIONS_FILE <- NULL  # Auto-find the existing audited coverage evidence TSV
GC_COLORS <- c(skew_positive='#00FF00', skew_negative='#FF0000', ratio_positive='#0000FF', ratio_negative='#FFFF00')
CREATE_GENBANK <- TRUE # Apply called mutations to the sequence, not just annotations
CREATE_DNA <- TRUE     # Convert changed GenBank files when SnapGene CLI is available
SNAPGENE_PATH <- NULL  # NULL: auto-detect; running SnapGene sessions are never closed
# Everything below is bundled implementation; no other code file is loaded.
.breseq_one_shot_config <- list(input_dir=INPUT_DIR, output_dir=OUTPUT_DIR, reference_file=REFERENCE_FILE, excel_file=EXCEL_FILE, gff_file=GFF_FILE, eggnog_file=EGGNOG_FILE, hide_intergenic_labels=HIDE_INTERGENIC_LABELS, sample_order=SAMPLE_ORDER, sample_labels=SAMPLE_LABELS, baseline_id=BASELINE_ID, sample_aliases=SAMPLE_ALIASES, display_deletions_file=DISPLAY_DELETIONS_FILE, gc_colors=GC_COLORS, create_genbank=CREATE_GENBANK, create_dna=CREATE_DNA, snapgene_path=SNAPGENE_PATH)
breseq_one_shot <- local({
# Bundled: breesq_drawer.R
#' breseqVisualization
#'
#' @param fna_count Number of contigs to process. If NULL, uses the number of rows in contig_length.
#' @param completeness Logical indicating if the genome is complete. If TRUE, saves output as PDF after interactive display.
#' @param mutation_data A data frame with mutation data (columns: start, end, type, color, ALE_version; or position, mutation, color, ALE_version for raw data).
#' @param hide_intergenic_labels Hide labels whose annotation starts with
#'   "intergenic" (default: FALSE). Mutation marks and data are retained.
#' @param gc_colors Optional named GC palette overrides: skew_positive,
#'   skew_negative, ratio_positive, and ratio_negative. Defaults are green/red
#'   for GC skew and blue/yellow for GC ratio. Signs refer to deviations from
#'   the existing average, and cached color columns are ignored.
#' @return A genome map displayed interactively in the plot window, with optional PDF output.
#' @export


breseqVisualization <- function(fna_count = NULL, completeness = NULL, mutation_data = NULL,
                                hide_intergenic_labels = FALSE, gc_colors = NULL) {
  gc_palette <- .breseq_gc_colors(gc_colors)
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
    if (length(mutation_list) > 0) {
      for (sample_name in names(mutation_list)) {
        data_for_parsing <- mutation_list[[sample_name]]
        if (nrow(data_for_parsing) > 0) {
          parsed_sample_data <- parse_mutation_data(data_for_parsing, Genome_summary$VERSION[1], chr_length) # Pass genome_id and chr_length

          # Add color column based on mutation type
          parsed_sample_data <- parsed_sample_data %>%
            dplyr::mutate(
              type_color = dplyr::case_when(
                type == "INS" ~ "#FF5C8D",
                type == "DEL" ~ "#03045E",
                type == "SNP" ~ "#139487",
                type == "REP" ~ "#E5C87B",
                TRUE ~ "gray" # Fallback for UNKNOWN
              )
            )
          mutation_data_list_parsed[[sample_name]] <- parsed_sample_data
        } else {
          mutation_data_list_parsed[[sample_name]] <- data.frame() # Keep empty if no data
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

    # Dynamically adjust xlim to cover all data points
    max_end <- chr_length$end

    # recording plot
    pdf_file <- paste0(seqid, "_map.pdf")
    grDevices::pdf(file = pdf_file, width = pdf_width*2, height = pdf_height*2)


    circos.initialize(factors = c(seqid), xlim = matrix(c(chr_length$start, max_end), ncol = 2))

    # Calculate genome length and determine step size and label units
    genome_length <- chr_length$end # Use chr_length$end as default
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
      if ("label" %in% colnames(mut_data) && any(!is.na(mut_data$label) & mut_data$label != "")) {
        label_data <- mut_data %>%
          dplyr::filter(!is.na(label) & label != "")
        label_data <- .breseq_filter_intergenic_labels(label_data, hide_intergenic_labels) %>%
          dplyr::select(seqid, start, end, label, version_color)
        if (nrow(label_data) > 0) {
          message("Adding labels for ", nrow(label_data), " mutations")
          # Debug: Check label_data contents
          message("label_data columns: ", paste(colnames(label_data), collapse = ", "))
          message("label_data version_color values: ", paste(head(label_data$version_color, 6), collapse = ", "))
          circos.genomicLabels(
            bed = label_data,
            labels.column = 4,
            side = "outside",
            col = label_data$version_color, # Use ALE_version gradient color
            line_col = label_data$version_color,
            cex = 0.6,
            connection_height = mm_h(5),
            line_lwd = 0.5
          )
        }
      }
    } else {
      message(paste("Skipping mutation tracks for", seqid, "due to missing mutation_data."))
    }

    # Plot cumulative mutation tracks and labels
    if (!is.null(mut_data) && nrow(mut_data) > 0 && "ALE_version" %in% colnames(mut_data)) {
      message("Using mut_data for tracks (placeholder until mutation_data is provided)")
      mutation_data <- mut_data  # Replace this with actual mutation_data when available

      message("Mutation data before processing: ", nrow(mutation_data), " rows")
      message("mutation_data columns: ", paste(colnames(mutation_data), collapse = ", "))
      message("mutation_data ALE_version distribution: ", paste(table(mutation_data$ALE_version), collapse = ", "))

      # Plot cumulative tracks for each ALE_version
      groups <- samples
      message("Mutation groups (ALE_version): ", paste(groups, collapse = ", "))
      track_index <- get.current.track.index() + 1

      # Set transparency (alpha: 0 = fully transparent, 1 = fully opaque)
      alpha_rect <- 1.0  # Transparency for mutation rectangles
      alpha_bg <- 0.3    # Transparency for track background

      for (i in seq_along(groups)) {
        ALE_group <- groups[i]
        group_data <- mutation_data %>% dplyr::filter(ALE_version <= ALE_group)
        bg_color <- version_colors[as.character(ALE_group)]

        message("Plotting cumulative track for ALE_version <= ", ALE_group, ", rows: ", nrow(group_data))

        if (nrow(group_data) > 0) {
          message("group_data columns: ", paste(colnames(group_data), collapse = ", "))
          message("group_data version_color values: ", paste(head(group_data$version_color, 6), collapse = ", "))
          message("group_data ALE_version distribution: ", paste(table(group_data$ALE_version), collapse = ", "))

          # Apply transparency to version_color for rectangles
          group_data$version_color_trans <- sapply(group_data$version_color, function(col) {
            adjustcolor(col, alpha.f = alpha_rect)
          })

          # Apply transparency to background color
          bg_color_trans <- adjustcolor(bg_color, alpha.f = alpha_bg)

          circos.track(
            factors = as.character(group_data$seqid),
            panel.fun = function(region, value, ...) {
              circos.genomicRect(
                region = matrix(c(as.numeric(group_data$start), as.numeric(group_data$end)), ncol = 2),
                value = as.character(group_data$type),
                col = group_data$version_color_trans, # Use transparent colors
                ybottom = 0,
                ytop = 1,
                border = NA
              )
            },
            ylim = c(0, 1),
            track.index = track_index,
            bg.col = bg_color_trans, # Use transparent background
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
      # Relies on 'window_nucleotideName_' existing
      gc_data <- get(paste("window", names(nucleotide)[i], sep = "_"), inherits = TRUE)
      if (!is.null(gc_data) && nrow(gc_data) > 0 && all(gc_data[[paste0(type, "_minus_average")]] != 0)) {
        if (any(!is.na(gc_data[[paste0(type, "_minus_average")]]))) {
          circlize::circos.track(factors = gc_data$seqid, # factors should match the sector
                                 panel.fun = function(x,y) { # x and y are implicit here for the current sector
                                   current_sector_data <- gc_data[gc_data$seqid == circlize::get.cell.meta.data("sector.index"),]
                                   if(nrow(current_sector_data) > 0){
                                     gc_values <- as.numeric(current_sector_data[[paste0(type, "_minus_average")]])
                                     palette_prefix <- if (type == "gc_skew") "skew" else "ratio"
                                     circlize::circos.barplot(
                                       value = gc_values,
                                       pos = (current_sector_data$start + current_sector_data$end) / 2,
                                       col = ifelse(gc_values >= 0,
                                         gc_palette[[paste0(palette_prefix, "_positive")]],
                                         gc_palette[[paste0(palette_prefix, "_negative")]]),
                                       border = NA,
                                       bar_width = 10000 # bar_width might need adjustment based on data scale
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

    # Draw legends for COG, mutations, and version_color
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

    # Mutation type legend
    if (!is.null(mut_data) && nrow(mut_data) > 0) {
      mut_legend <- mut_data %>% group_by(type) %>% slice_head() %>% select(type, color) %>% na.omit()
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

      # Version_color legend (NEW)
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
    #readline(prompt = "Press Enter to continue to the next contig...")

    grDevices::dev.off()
    message("Saved PDF for ", seqid, " as ", pdf_file)
  }
  circos.clear()
}

# Bundled: circos_visualization.R
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
#' @param hide_intergenic_labels Hide labels whose annotation starts with
#'   "intergenic" (default: FALSE). Mutation marks and data are retained.
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
                               seqid = NULL,
                               hide_intergenic_labels = FALSE) {

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
      stage_columns <- .breseq_stage_columns(mutation_data)
      samples <- attr(mutation_data, "breseq_sample_order", exact = TRUE)
      if (is.null(samples)) {
        samples <- unique(stage_columns$sample)
        if (length(samples) && requireNamespace("gtools", quietly = TRUE)) {
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
      dplyr::filter(new == TRUE & !is.na(label) & label != "" & !is.na(start) & !is.na(end))
    label_data <- .breseq_filter_intergenic_labels(label_data, hide_intergenic_labels) %>%
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


# Filter the label layer only; callers retain their complete mutation data.
.breseq_filter_intergenic_labels <- function(data, hide_intergenic_labels = FALSE) {
  if (!isTRUE(hide_intergenic_labels) || !"annotation" %in% names(data)) {
    return(data)
  }

  annotation <- trimws(gsub("\u00a0", " ", as.character(data$annotation), fixed = TRUE))
  intergenic <- !is.na(annotation) &
    grepl("^intergenic([[:space:]]|$)", annotation, ignore.case = TRUE)
  data[!intergenic, , drop = FALSE]
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

    # Match complete stage IDs, including suffixes such as 3w and embedded dots.
    stage_columns <- .breseq_stage_columns(df, sample_order = samples)
    mutation_columns <- stage_columns[stage_columns$field == "mutation", , drop = FALSE]
    mut_cols <- names(df)[mutation_columns$index]

    if (length(mut_cols) == 0) {
      # Try simple column names
      if ("mutation" %in% names(df)) {
        mut_cols <- "mutation"
      } else {
        stop("No mutation columns found")
      }
    }

    message("  Mutation columns found: ", paste(mut_cols, collapse = ", "))

    # Build long format data
    long_data <- data.frame()

    for (mut_col in mut_cols) {
      sample_name <- mutation_columns$sample[match(mut_col, names(df)[mutation_columns$index])]
      if (!length(sample_name)) sample_name <- NA_character_
      if (is.na(sample_name)) next

      # Find corresponding annotation and gene columns
      ann_col <- stage_columns$index[stage_columns$field == "annotation" & stage_columns$sample == sample_name]
      gene_col <- stage_columns$index[stage_columns$field == "gene" & stage_columns$sample == sample_name]

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

# Bundled: classify_mutations.R
#' Classify Mutation Types
#'
#' Analyzes mutation strings from breseq output and classifies them into
#' standardized categories: SNP, insertion, deletion, substitution, amplification, inversion, or mobile element.
#'
#' @param mutation_string Character string describing the mutation (e.g., "A→G", "Δ1 bp", "+GC")
#' @param annotation Optional annotation string for additional context
#' @return Character string indicating mutation type
#' @export
#' @examples
#' classify_mutation("A→G")           # "SNP"
#' classify_mutation("Δ500 bp")       # "large_deletion"
#' classify_mutation("+ATCG")         # "insertion"
#' classify_mutation("(G)8→9")        # "insertion"  (repeat expansion)
#' classify_mutation("(G)9→8")        # "deletion"   (repeat contraction)
classify_mutation <- function(mutation_string, annotation = NULL) {

  if (is.na(mutation_string) || is.null(mutation_string) ||
      tolower(mutation_string) %in% c("na", "nan", "unknown", "")) {
    return("unknown")
  }

  mut <- as.character(mutation_string)

  # Mobile element insertion (IS element)
  if (grepl("IS\\d+|mobile.element|transpos", mut, ignore.case = TRUE)) {
    return("mobile_element")
  }

  # Large deletion (Δ notation)
  if (grepl("^Δ|Δ\\d+", mut)) {
    # Extract size if available
    size_match <- regmatches(mut, regexpr("[0-9][0-9,]*", mut))
    if (length(size_match) > 0) {
      size <- as.numeric(gsub(",", "", size_match[1]))
      if (size >= 50) return("large_deletion")
    }
    return("deletion")
  }

  # Simple insertion (+bases)
  if (grepl("^\\+[ATCG]+$", mut, ignore.case = TRUE)) {
    inserted <- gsub("^\\+", "", mut)
    if (nchar(inserted) >= 50) return("large_insertion")
    return("insertion")
  }

  # Simple deletion (-bases)
  if (grepl("^-[ATCG]+$", mut, ignore.case = TRUE)) {
    deleted <- gsub("^-", "", mut)
    if (nchar(deleted) >= 50) return("large_deletion")
    return("deletion")
  }

  # Repeat expansion/contraction: (G)8→9 or (ATCG)5→3
  if (grepl("\\([ATCG]+\\)\\d+.*→.*\\d+", mut, ignore.case = TRUE)) {
    # Extract numbers
    nums <- as.numeric(regmatches(mut, gregexpr("\\d+", mut))[[1]])
    if (length(nums) >= 2) {
      if (nums[2] > nums[1]) {
        diff <- nums[2] - nums[1]
        if (diff >= 50) return("large_insertion")
        return("insertion")  # repeat expansion
      } else {
        diff <- nums[1] - nums[2]
        if (diff >= 50) return("large_deletion")
        return("deletion")   # repeat contraction
      }
    }
    return("indel")
  }

  # SNP: single base change A→G, T→C, etc.
  if (grepl("^[ATCG]→[ATCG]$", mut, ignore.case = TRUE)) {
    return("SNP")
  }

  # Complex substitution: multiple bases → multiple bases
  if (grepl("[ATCG]+→[ATCG]+", mut, ignore.case = TRUE)) {
    parts <- strsplit(mut, "→")[[1]]
    if (length(parts) == 2) {
      old_len <- nchar(gsub("[^ATCGatcg]", "", parts[1]))
      new_len <- nchar(gsub("[^ATCGatcg]", "", parts[2]))
      if (old_len == new_len) {
        if (old_len == 1) return("SNP")
        return("substitution")
      } else if (new_len > old_len) {
        return("insertion")
      } else {
        return("deletion")
      }
    }
  }

  # Amplification
  if (grepl("amplification|duplication|×\\d+", mut, ignore.case = TRUE)) {
    return("amplification")
  }

  # Inversion
  if (grepl("inversion|inv", mut, ignore.case = TRUE)) {
    return("inversion")
  }

  # Check annotation for additional clues
  if (!is.null(annotation) && !is.na(annotation)) {
    ann_lower <- tolower(annotation)
    if (grepl("intergenic|promoter|upstream|downstream", ann_lower)) {
      # Intergenic regions - type determined by mutation itself
    }
    if (grepl("pseudogene", ann_lower)) {
      # Still classify by mutation type

    }
  }

  # Missing coverage (from separate sheet)
  if (grepl("missing.coverage|missing_coverage", mut, ignore.case = TRUE)) {
    return("missing_coverage")
  }

  return("unknown")
}

#' Classify Multiple Mutations (Vectorized)
#'
#' @param mutations Character vector of mutation strings
#' @param annotations Optional character vector of annotations
#' @return Character vector of mutation types
#' @export
classify_mutations <- function(mutations, annotations = NULL) {
  if (is.null(annotations)) {
    annotations <- rep(NA, length(mutations))
  }
  mapply(classify_mutation, mutations, annotations, USE.NAMES = FALSE)
}

#' Get Mutation Type Color
#'
#' Returns a color for each mutation type for visualization.
#'
#' @param mutation_type Character string of mutation type
#' @return Hex color code
#' @export
get_mutation_color <- function(mutation_type) {
  colors <- c(
    "SNP"              = "#3498db",  # Blue
    "insertion"        = "#2ecc71",  # Green
    "deletion"         = "#e74c3c",  # Red
    "large_insertion"  = "#27ae60",  # Dark green
    "large_deletion"   = "#c0392b",  # Dark red
    "substitution"     = "#9b59b6",  # Purple
    "amplification"    = "#f39c12",  # Orange
    "inversion"        = "#1abc9c",  # Teal
    "mobile_element"   = "#e67e22",  # Dark orange
    "missing_coverage" = "#95a5a6",  # Gray
    "indel"            = "#f1c40f",  # Yellow
    "unknown"          = "#bdc3c7"   # Light gray
  )

  ifelse(mutation_type %in% names(colors),
         colors[mutation_type],
         colors["unknown"])
}

#' Add Mutation Type Column to Data Frame
#'
#' Adds a 'mutation_type' column to a data frame containing mutation data.
#' When the requested mutation column is absent, the first non-empty canonical
#' stage mutation column in each row supplies its classification. Source columns
#' are retained, including later allele changes.
#'
#' @param df Data frame with a 'mutation' column
#' @param mutation_col Name of the column containing mutation strings
#' @param annotation_col Optional name of annotation column
#' @return Data frame with added 'mutation_type' column
#' @export
add_mutation_types <- function(df, mutation_col = "mutation", annotation_col = "annotation") {

  mutation_cols <- mutation_col
  if (!mutation_col %in% names(df)) {
    # Try to find mutation column
    mut_cols <- grep("^(mutation|mut)\\.", names(df), value = TRUE, ignore.case = TRUE)
    if (length(mut_cols) > 0) {
      mutation_cols <- mut_cols
    } else {
      warning("No mutation column found. Returning original data frame.")
      return(df)
    }
  }

  mutations <- rep(NA_character_, nrow(df))
  annotations <- rep(NA_character_, nrow(df))
  for (col in mutation_cols) {
    values <- as.character(df[[col]])
    missing <- is.na(mutations) | !nzchar(trimws(mutations))
    available <- !is.na(values) & nzchar(trimws(values)) &
      !tolower(trimws(values)) %in% c("na", "nan")
    take <- missing & available
    mutations[take] <- values[take]
    matching_annotation <- sub("^(mutation|mut)([.]|$)", paste0(annotation_col, "\\2"), col)
    if (matching_annotation %in% names(df)) {
      annotations[take] <- as.character(df[[matching_annotation]])[take]
    } else if (annotation_col %in% names(df)) {
      annotations[take] <- as.character(df[[annotation_col]])[take]
    }
  }

  df$mutation_type <- classify_mutations(mutations, annotations)
  df$mutation_color <- sapply(df$mutation_type, get_mutation_color)

  return(df)
}

# Bundled: compare_multiple_tables.R
#' Sequentially Compare Multiple Lists of Tables
#'
#' This function sequentially compares multiple named lists of data frames. It is designed
#' to compare lists in a stepwise manner (e.g., A to B, B to C, etc.), storing the results
#' of each comparison.
#'
#' @param list_of_tables A list of named lists of data frames. Each named list should have identical names
#' across all entries to ensure accurate comparisons.
#' @param sample_order Optional character vector containing every input sample exactly once.
#' @return A named list where each entry corresponds to the comparison results of a pair of lists.
#' Each comparison result is a named list of data frames, as generated by `compare_named_tables`.
#' @export
#' @examples
#' # Example usage:
#' list_of_tables <- list(A = tables_A, B = tables_B, C = tables_C)
#' all_comparisons <- compare_multiple_tables(list_of_tables)
compare_multiple_tables <- function(list_of_tables, sample_order = NULL) {
  # Check if the input is a list of named lists
  if (!is.list(list_of_tables) || length(list_of_tables) < 2) {
    stop("The input must be a list containing at least two named lists of data frames.")
  }

  # Initialize an empty list to store all comparison results
  all_comparisons <- list()

  # Get sample names from list keys
  sample_names <- names(list_of_tables)
  if (is.null(sample_names) || anyNA(sample_names) || any(!nzchar(sample_names)) || anyDuplicated(sample_names))
    stop("Input tables must have distinct non-empty sample names.")
  if (!is.null(sample_order)) {
    if (anyDuplicated(sample_order) || !setequal(sample_order, sample_names))
      stop("sample_order must contain every input sample exactly once.")
    sample_names <- sample_order
    list_of_tables <- list_of_tables[sample_names]
  }

  # Sequentially compare each pair of tables
  for (i in 1:(length(list_of_tables) - 1)) {
    name_A <- sample_names[i]
    name_B <- sample_names[i + 1]

    cat("Comparing:", name_A, "to", name_B, "\n")

    # Extract the current and next lists to compare
    list_A <- list_of_tables[[i]]
    list_B <- list_of_tables[[i + 1]]

    # Perform the comparison with actual sample names as suffixes
    comparison_results <- compare_named_tables(
      list_A, list_B,
      name_A = paste0(name_A, "/output/index"),
      name_B = paste0(name_B, "/output/index")
    )

    # Store the comparison results with a meaningful name
    comparison_name <- paste0("Comparison_", name_A, "_to_", name_B)
    all_comparisons[[comparison_name]] <- comparison_results
  }

  return(all_comparisons)
}

# Bundled: compare_named_tables.R
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

# Bundled: excel_presentation.R
# Shared presentation helpers. Canonical data and sample IDs are never renamed.
.breseq_text <- function(x) trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))

.breseq_stage_columns <- function(data, sample_order = NULL) {
  fields <- c("Forward aligned reads", "Reverse aligned reads", "Read Coverage Depth",
              "aligned reads1", "aligned reads2", "reads (cov).1", "reads (cov)",
              "description", "annotation", "evidence", "mutation", "product",
              "seq id", "seq_id", "gene", "end", "size", "\u2190reads", "reads\u2192",
              "score", "skew", "freq")
  fields <- fields[order(nchar(fields), decreasing = TRUE)]
  found <- list()
  for (i in seq_along(names(data))) {
    nm <- .breseq_text(names(data)[i])
    if (nm %in% fields) next
    for (field in fields) {
      prefix <- paste0(field, ".")
      if (!startsWith(nm, prefix)) next
      sample <- substring(nm, nchar(prefix) + 1L)
      sample <- sub("/output(?:/index)?$", "", sample, perl = TRUE)
      if (nzchar(sample)) found[[length(found) + 1L]] <-
        data.frame(index = i, field = if (field == "seq_id") "seq id" else field,
                   sample = sample, stringsAsFactors = FALSE)
      break
    }
  }
  result <- if (length(found)) do.call(rbind, found) else
    data.frame(index = integer(), field = character(), sample = character())
  detected <- unique(result$sample)
  if (!is.null(sample_order)) {
    if (!is.character(sample_order) || anyNA(sample_order) || any(!nzchar(sample_order)) ||
        anyDuplicated(sample_order) || !all(detected %in% sample_order))
      stop("sample_order must contain each detected sample exactly once (additional empty baselines are allowed).")
  }
  if (anyDuplicated(paste(result$field, result$sample, sep = "\r")))
    stop("Ambiguous duplicate field/sample columns in the wide table.")
  attr(result, "sample_order") <- if (is.null(sample_order)) detected else sample_order
  result
}

.breseq_table_kind <- function(name) {
  name <- tolower(name)
  if (grepl("missing|coverage|^source mc$", name)) return("missing")
  if (grepl("junction", name)) return("junction")
  if (grepl("mutation", name)) return("mutation")
  "other"
}

.breseq_seq <- function(x) sub("\\.[0-9]+$", "", .breseq_text(x))
.breseq_coord <- function(x) gsub("[,[:space:]]", "", .breseq_text(x))
.breseq_allele <- function(x) gsub("\u2206", "\u0394", .breseq_coord(x), fixed = TRUE)

# Keys are kept as an attribute, so legacy exported columns remain unchanged.
.breseq_source_keys <- function(data, kind, default_seq = "") {
  n <- nrow(data)
  if (!n) return(character())
  if ("comparison_key" %in% names(data)) return(as.character(data$comparison_key))
  seq_col <- intersect(c("seq id", "seq_id", "seqid"), names(data))
  seqid <- if (length(seq_col)) .breseq_seq(data[[seq_col[1]]]) else rep(default_seq, n)
  seqid[is.na(seqid)] <- default_seq
  coordinate <- if (kind == "missing") "start" else "position"
  endpoint <- paste(seqid, .breseq_coord(data[[coordinate]]), sep = "::")
  if (kind != "junction") return(endpoint)
  if (n %% 2L) stop("Junction tables must contain two endpoint rows per event.")
  keys <- character(n)
  for (i in seq.int(1L, n, by = 2L)) {
    pair <- endpoint[i:(i + 1L)]
    ordering <- order(pair)
    roles <- order(ordering)
    keys[i:(i + 1L)] <- paste(paste(pair[ordering], collapse = "||"),
                             paste0("end", roles), sep = "||")
  }
  keys
}

.breseq_wide_seq <- function(data, columns) {
  result <- rep(NA_character_, nrow(data))
  globals <- intersect(c("seq_id", "seq id", "seqid"), names(data))
  indices <- unique(c(match(globals, names(data)), columns$index[columns$field == "seq id"]))
  for (i in indices) {
    x <- .breseq_seq(data[[i]])
    ok <- !is.na(x) & nzchar(x)
    conflict <- ok & !is.na(result) & result != x
    if (any(conflict)) stop("A merged row contains conflicting replicon IDs; rebuild comparisons with composite keys.")
    result[ok & is.na(result)] <- x[ok & is.na(result)]
  }
  result[is.na(result)] <- ""
  result
}

.breseq_field <- function(data, columns, field, sample) {
  at <- columns$index[columns$field == field & columns$sample == sample]
  if (!length(at)) return(rep(NA_character_, nrow(data)))
  .breseq_text(data[[at[1]]])
}
.breseq_present <- function(x) !is.na(x) & nzchar(x)
.breseq_bounds <- function(x) {
  if (!.breseq_present(x)) return(c(NA_real_, NA_real_))
  values <- as.numeric(regmatches(gsub(",", "", x, fixed = TRUE),
                                 gregexpr("[0-9]+", gsub(",", "", x, fixed = TRUE)))[[1]])
  if (!length(values) || anyNA(values)) c(NA_real_, NA_real_) else range(values)
}

# Return row-by-stage states without modifying data or using changing read metrics
# as event identity. MC states describe evidence, never newly acquired deletions.
.breseq_stage_status <- function(data, table_name, sample_order = NULL, mutation_data = NULL) {
  columns <- .breseq_stage_columns(data, sample_order)
  samples <- attr(columns, "sample_order")
  states <- matrix("not_reported", nrow(data), length(samples), dimnames = list(NULL, samples))
  kind <- .breseq_table_kind(table_name)
  if (!nrow(data) || !length(samples) || kind == "other") return(list(samples = samples, states = states))
  seqid <- .breseq_wide_seq(data, columns)
  if (kind == "mutation") {
    if (!"position" %in% names(data)) stop("Mutation tables require position for stage highlighting.")
    locus <- paste(seqid, .breseq_coord(data$position), sep = "::")
    if (anyDuplicated(locus)) stop("Duplicate mutation loci cannot be highlighted safely; rebuild the merged table.")
    history <- vector("list", nrow(data))
    previous <- rep(NA_character_, nrow(data))
    for (sample in samples) {
      allele <- .breseq_allele(.breseq_field(data, columns, "mutation", sample))
      for (i in which(.breseq_present(allele))) {
        states[i, sample] <- if (.breseq_present(previous[i])) {
          if (allele[i] == previous[i]) "retained" else "changed"
        } else if (allele[i] %in% history[[i]]) "reappeared" else "new"
        history[[i]] <- unique(c(history[[i]], allele[i]))
      }
      previous <- allele
    }
  } else if (kind == "junction") {
    keys <- if ("comparison_key" %in% names(data)) as.character(data$comparison_key) else attr(data, "breseq_keys")
    if (is.null(keys) || length(keys) != nrow(data) || anyNA(keys) || any(!grepl("\\|\\|", keys)))
      stop("Junction highlighting requires paired comparison keys from compare_multiple_tables().")
    event <- sub("\\|\\|end[12]$", "", keys)
    groups <- split(seq_len(nrow(data)), event)
    if (any(lengths(groups) != 2L)) stop("Every junction event must retain both endpoint rows.")
    previous <- character(); seen <- character()
    for (sample in samples) {
      present <- .breseq_present(.breseq_field(data, columns, "seq id", sample))
      if (!any(columns$field == "seq id" & columns$sample == sample))
        present <- .breseq_present(.breseq_field(data, columns, "reads (cov)", sample))
      current <- character()
      for (e in names(groups)) {
        at <- groups[[e]]
        if (any(present[at]) && !all(present[at])) stop("Only one endpoint is present for a junction sample/event.")
        if (!all(present[at])) next
        states[at, sample] <- if (e %in% previous) "retained" else if (e %in% seen) "reappeared" else "new"
        current <- c(current, e)
      }
      seen <- union(seen, current); previous <- current
    }
  } else if (kind == "missing") {
    if (!"start" %in% names(data)) stop("Missing-coverage tables require start for highlighting.")
    previous <- list(); history <- list()
    mcols <- if (!is.null(mutation_data)) .breseq_stage_columns(mutation_data, sample_order) else NULL
    mseq <- if (!is.null(mcols)) .breseq_wide_seq(mutation_data, mcols) else NULL
    for (j in seq_along(samples)) {
      sample <- samples[j]
      ends <- .breseq_field(data, columns, "end", sample)
      genes <- .breseq_field(data, columns, "gene", sample)
      current <- list()
      for (i in which(.breseq_present(ends))) {
        a <- .breseq_bounds(as.character(data$start[i])); b <- .breseq_bounds(ends[i])
        record <- list(seq = seqid[i], start = a[2], end = b[1],
                       raw = paste(.breseq_coord(data$start[i]), .breseq_coord(ends[i])), gene = genes[i])
        if (anyNA(c(record$start, record$end)) || record$start > record$end)
          stop("Cannot classify an invalid missing-coverage interval.")
        overlap <- function(other) {
          common_locus <- TRUE
          if (.breseq_present(record$gene) && .breseq_present(other$gene)) {
            loci <- function(x) regmatches(x, gregexpr("[A-Za-z0-9]+_RS[0-9]+", x))[[1]]
            left <- loci(record$gene); right <- loci(other$gene)
            common_locus <- if (length(left) && length(right)) length(intersect(left, right)) > 0L else
              identical(record$gene, other$gene)
          }
          identical(record$seq, other$seq) && common_locus &&
            max(record$start, other$start) <= min(record$end, other$end)
        }
        old <- Filter(overlap, previous)
        state <- if (length(old)) {
          if (any(vapply(old, function(o) identical(record$raw, o$raw), logical(1)))) "retained" else "updated"
        } else if (length(Filter(overlap, history))) "reappeared" else "new"
        if (state == "new" && j > 1L && !is.null(mcols)) {
          alleles <- .breseq_field(mutation_data, mcols, "mutation", samples[j - 1L])
          is_del <- !is.na(alleles) & grepl("^[\u0394\u2206][0-9,]+[[:space:]]*bp$", alleles)
          size <- suppressWarnings(as.numeric(gsub("[^0-9]", "", alleles)))
          pos <- suppressWarnings(as.numeric(gsub(",", "", mutation_data$position, fixed = TRUE)))
          covered <- is_del & mseq == record$seq & pos <= record$start & pos + size - 1 >= record$end
          if (any(covered, na.rm = TRUE)) state <- "prior_del_region"
        }
        states[i, sample] <- state
        current[[length(current) + 1L]] <- record
      }
      history <- c(history, current); previous <- current
    }
  }
  list(samples = samples, states = states)
}

.breseq_validate_labels <- function(labels, samples) {
  if (is.null(labels)) return(stats::setNames(samples, samples))
  if (!is.character(labels) || is.null(names(labels)) || anyNA(labels) || any(!nzchar(labels)) ||
      anyNA(names(labels)) || any(!nzchar(names(labels))) || anyDuplicated(names(labels)) ||
      !all(names(labels) %in% samples)) stop("sample_labels must be a named character map of known sample IDs.")
  result <- stats::setNames(samples, samples); result[names(labels)] <- labels
  if (anyDuplicated(unname(result))) stop("sample_labels must produce distinct display labels.")
  result
}

.breseq_write_presentation <- function(wb, sheet, data, table_name, highlight_new = FALSE,
                                     clean_headers = FALSE, sample_labels = NULL,
                                     sample_order = NULL, mutation_data = NULL) {
  sheet <- unname(as.character(sheet))
  if (!is.logical(highlight_new) || length(highlight_new) != 1L || is.na(highlight_new) ||
      !is.logical(clean_headers) || length(clean_headers) != 1L || is.na(clean_headers))
    stop("highlight_new and clean_headers must be single TRUE/FALSE values.")
  if (!highlight_new && !clean_headers) {
    openxlsx::writeData(wb, sheet, data)
    return(invisible(NULL))
  }
  openxlsx::modifyBaseFont(wb, fontSize = 10, fontColour = "#333333", fontName = "Arial")
  body_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10, fontColour = "#333333")
  header_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10,
    fontColour = "#333333", textDecoration = "bold", fgFill = "#FFFFFF")
  columns <- .breseq_stage_columns(data, sample_order)
  samples <- attr(columns, "sample_order")
  labels <- .breseq_validate_labels(sample_labels, samples)
  displayed <- data
  if (clean_headers) {
    base <- c(seq_id = "Replicon", position = "Position (bp)", start = "Start (bp)", comparison_key = "Comparison key")
    at <- match(names(displayed), names(base)); use <- !is.na(at)
    names(displayed)[use] <- unname(base[at[use]])
    for (k in seq_len(nrow(columns))) {
      field <- columns$field[k]
      field <- paste0(toupper(substring(field, 1, 1)), substring(field, 2))
      names(displayed)[columns$index[k]] <- paste(unname(labels[columns$sample[k]]), field, sep = " | ")
    }
    if (anyDuplicated(names(displayed))) stop("Clean headers would be ambiguous.")
    metadata <- data.frame(sheet = sheet, column = seq_along(data), canonical = names(data),
                           display = names(displayed), sample_order = paste(samples, collapse = "\t"),
                           stringsAsFactors = FALSE)
    meta_sheet <- "_breseq_columns"
    if (!meta_sheet %in% names(wb)) {
      openxlsx::addWorksheet(wb, meta_sheet, visible = FALSE)
      openxlsx::writeData(wb, meta_sheet, metadata)
      metadata_start <- 2L
    } else {
      old <- openxlsx::readWorkbook(wb, sheet = meta_sheet, check.names = FALSE)
      metadata_start <- nrow(old) + 2L
      openxlsx::writeData(wb, meta_sheet, metadata, startRow = metadata_start, colNames = FALSE)
    }
    openxlsx::addStyle(wb, meta_sheet, header_style, rows = 1L,
      cols = seq_len(ncol(metadata)), gridExpand = TRUE)
    if (nrow(metadata)) openxlsx::addStyle(wb, meta_sheet, body_style,
      rows = seq.int(metadata_start, length.out = nrow(metadata)),
      cols = seq_len(ncol(metadata)), gridExpand = TRUE)
  }
  openxlsx::writeData(wb, sheet, displayed)
  if (ncol(displayed)) {
    openxlsx::addStyle(wb, sheet, header_style, rows = 1L,
      cols = seq_len(ncol(displayed)), gridExpand = TRUE)
    if (nrow(displayed)) openxlsx::addStyle(wb, sheet, body_style,
      rows = seq_len(nrow(displayed)) + 1L, cols = seq_len(ncol(displayed)), gridExpand = TRUE)
  }
  if (highlight_new && nrow(data)) {
    status <- .breseq_stage_status(data, table_name, sample_order, mutation_data)
    # Color means observed in this sample, including retained calls. Bold marks
    # a first report, changed allele or reappearance; it is not a de novo claim.
    palette <- c("#EEF3F8", "#F0F5EF", "#F7F2EA", "#F3F0F6")
    is_wt <- toupper(status$samples) == "WT"
    fills <- stats::setNames(rep("#F3F3F3", length(status$samples)), status$samples)
    fills[!is_wt] <- rep(palette, length.out = sum(!is_wt))
    bold_style <- openxlsx::createStyle(fontName = "Arial", fontSize = 10,
      fontColour = "#333333", textDecoration = "bold")
    for (sample in status$samples) {
      rows <- which(status$states[, sample] != "not_reported")
      cols <- columns$index[columns$sample == sample]
      if (length(rows) && length(cols)) openxlsx::addStyle(wb, sheet,
        openxlsx::createStyle(fgFill = fills[[sample]]), rows = rows + 1L, cols = cols,
        gridExpand = TRUE, stack = TRUE)
      bold_rows <- which(status$states[, sample] %in% c("new", "changed", "reappeared"))
      if (length(bold_rows) && length(cols)) openxlsx::addStyle(wb, sheet,
        bold_style, rows = bold_rows + 1L, cols = cols, gridExpand = TRUE, stack = TRUE)
    }
  }
  invisible(NULL)
}

#' Restore Canonical Headers from a Presentation Export
#'
#' Internal loader used for clean-header workbook round trips. Files without
#' breseqConverter metadata are returned unchanged.
#' @keywords internal
restore_canonical_export_headers <- function(file, sheet, data) {
  metadata_sheet <- "_breseq_columns"
  if (!metadata_sheet %in% openxlsx::getSheetNames(file)) return(data)
  meta <- openxlsx::read.xlsx(file, sheet = metadata_sheet, check.names = FALSE)
  if (!all(c("sheet", "column", "canonical", "display", "sample_order") %in% names(meta)))
    stop("Invalid breseqConverter canonical-column metadata.")
  if (is.numeric(sheet)) sheet <- openxlsx::getSheetNames(file)[sheet]
  meta <- meta[meta$sheet == sheet, , drop = FALSE]
  if (!nrow(meta)) return(data)
  meta <- meta[order(as.integer(meta$column)), , drop = FALSE]
  if (nrow(meta) != ncol(data) || !identical(as.integer(meta$column), seq_len(ncol(data))))
    stop("Canonical-column metadata does not match the loaded worksheet; preserve empty columns when reading.")
  displayed_names <- as.character(meta$display)
  # read.xlsx(check.names=FALSE) still uses sep.names='.' by default.
  if (!identical(names(data), displayed_names) &&
      !identical(names(data), gsub(" ", ".", displayed_names, fixed = TRUE)))
    stop("Worksheet headers differ from the presentation metadata; refusing to assign incorrect canonical names.")
  names(data) <- as.character(meta$canonical)
  order <- unique(as.character(meta$sample_order))
  order <- order[!is.na(order) & nzchar(order)]
  if (length(order) == 1L) attr(data, "breseq_sample_order") <- strsplit(order, "\t", fixed = TRUE)[[1]]
  data
}

# Bundled: export_genbank.R
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
#' @keywords internal
.get_genbank_mutation_color <- function(mutation_type) {
  colors <- c(
    "INS" = "#FF5C8D",   # Bright Pink
    "DEL" = "#03045E",   # Dark Navy
    "SNP" = "#139487",   # Teal
    "REP" = "#E5C87B",   # Golden Yellow
    "UNKNOWN" = "#808080", # Gray
    "insertion" = "#FF5C8D",
    "large_insertion" = "#FF5C8D",
    "deletion" = "#03045E",
    "large_deletion" = "#03045E",
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
#' @keywords internal
.classify_genbank_mutation <- function(mutation) {
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
#' @param mutation_type Type of mutation (from .classify_genbank_mutation)
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
        restore_canonical_export_headers(mutation_data, sheet,
          openxlsx::read.xlsx(mutation_data, sheet = sheet, check.names = FALSE,
                             skipEmptyCols = FALSE)),
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
          mut_type <- if (!is.na(mutation)) .classify_genbank_mutation(mutation) else "unknown"
        }
      } else {
        mut_type <- if (!is.na(mutation)) .classify_genbank_mutation(mutation) else "unknown"
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
        color <- .get_genbank_mutation_color(mut_type)
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
      restore_canonical_export_headers(mutation_xlsx, sheet,
        openxlsx::read.xlsx(mutation_xlsx, sheet = sheet, check.names = FALSE,
                           skipEmptyCols = FALSE)),
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
      mut_type <- .classify_genbank_mutation(mutation)

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
        color <- .get_genbank_mutation_color(mut_type)
      } else if (!is.na(first_stage) && length(stage_names) > 0) {
        stage_idx <- match(first_stage, stage_names)
        if (is.na(stage_idx)) stage_idx <- 1
        color <- get_stage_color(stage_idx - 1, length(stage_names))
      } else {
        color <- .get_genbank_mutation_color(mut_type)
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

# Bundled: extract_named_tables.R
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
    tables <- html_data |> rvest::html_table(fill = TRUE)
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

# Bundled: find_snapgene.R
#' Find SnapGene Executable Path
#'
#' Automatically detects SnapGene installation across different operating systems.
#' Supports macOS, Windows, and Linux with common installation paths.
#'
#' @param custom_path Optional custom path to SnapGene executable
#' @return Character string with the path to SnapGene executable
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect
#' snapgene_path <- find_snapgene()
#'
#' # Use custom path
#' snapgene_path <- find_snapgene("/custom/path/to/snapgene")
#'
#' # Set via environment variable
#' Sys.setenv(SNAPGENE_PATH = "/path/to/snapgene")
#' snapgene_path <- find_snapgene()
#' }
find_snapgene <- function(custom_path = getOption("breseqConverter.snapgene_path")) {

  # 1. User-specified path takes priority
  if (!is.null(custom_path) && file.exists(custom_path)) {
    message("Using custom SnapGene path: ", custom_path)
    return(normalizePath(custom_path))
  }

  # 2. Check environment variable
  env_path <- Sys.getenv("SNAPGENE_PATH")
  if (nzchar(env_path) && file.exists(env_path)) {
    message("Using SnapGene from SNAPGENE_PATH: ", env_path)
    return(normalizePath(env_path))
  }

  # 3. OS-specific default paths
  os <- Sys.info()["sysname"]

  candidates <- switch(os,
    "Darwin" = c(
      "/Applications/SnapGene.app/Contents/MacOS/SnapGene",
      path.expand("~/Applications/SnapGene.app/Contents/MacOS/SnapGene"),
      "/Applications/SnapGene Viewer.app/Contents/MacOS/SnapGene Viewer",
      path.expand("~/Applications/SnapGene Viewer.app/Contents/MacOS/SnapGene Viewer")
    ),
    "Windows" = {
      prog_files <- Sys.getenv("PROGRAMFILES")
      prog_files_x86 <- Sys.getenv("PROGRAMFILES(X86)")
      local_app <- Sys.getenv("LOCALAPPDATA")
      c(
        file.path(prog_files, "SnapGene", "snapgene.exe"),
        file.path(prog_files_x86, "SnapGene", "snapgene.exe"),
        file.path(local_app, "Programs", "SnapGene", "snapgene.exe"),
        file.path(prog_files, "SnapGene Viewer", "snapgene.exe"),
        file.path(prog_files_x86, "SnapGene Viewer", "snapgene.exe")
      )
    },
    "Linux" = c(
      "/usr/bin/snapgene",
      "/usr/local/bin/snapgene",
      "/opt/snapgene/snapgene",
      "/opt/gslbiotech/snapgene/snapgene.sh",
      path.expand("~/snapgene/snapgene"),
      "/snap/bin/snapgene"
    ),
    character(0)
  )

  # 4. Check PATH
  path_result <- Sys.which("snapgene")
  if (nzchar(path_result)) {
    candidates <- c(unname(path_result), candidates)
  }

  # 5. Return first existing path

  for (path in candidates) {
    if (!is.na(path) && nzchar(path) && file.exists(path)) {
      message("SnapGene found: ", path)
      return(normalizePath(path))
    }
  }

  stop(
    "SnapGene not found. Please either:\n",
    "  1. Set environment variable: Sys.setenv(SNAPGENE_PATH = '/path/to/snapgene')\n",
    "  2. Set R option: options(breseqConverter.snapgene_path = '/path/to/snapgene')\n",
    "  3. Pass custom_path argument directly\n",
    "  4. Install SnapGene in the default location\n\n",
    "Searched locations:\n  ", paste(candidates, collapse = "\n  ")
  )
}

#' Check if SnapGene is available
#'
#' @return Logical indicating if SnapGene is found
#' @export
has_snapgene <- function() {
  tryCatch({
    find_snapgene()
    TRUE
  }, error = function(e) FALSE)
}

# Bundled: full_workflow.R
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
#' @param hide_intergenic_labels Hide intergenic annotation labels in Circos
#'   plots without removing mutation marks or data (default: FALSE).
#' @param highlight_new Shade every reported call by stage and bold new or
#'   changed calls in Excel exports (default: FALSE).
#' @param clean_headers Use concise Excel column headers (default: FALSE).
#' @param sample_labels Optional named vector mapping sample IDs to display labels
#'   in Excel exports.
#' @param sample_order Optional sample IDs in analysis order. When NULL, an order
#'   saved in Excel metadata is used, otherwise samples are naturally sorted.
#' @param gc_colors Optional named color overrides: skew_positive, skew_negative,
#'   ratio_positive, ratio_negative. Defaults are green/red (#00FF00/#FF0000) for
#'   GC skew and blue/yellow (#0000FF/#FFFF00) for GC ratio. Colors follow the sign
#'   of deviation from the window mean; zero uses the positive color.
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
                              create_genbank = TRUE,
                              hide_intergenic_labels = FALSE,
                              highlight_new = FALSE,
                              clean_headers = FALSE,
                              sample_labels = NULL,
                              sample_order = NULL,
                              gc_colors = NULL) {

  gc_colors <- .breseq_gc_colors(gc_colors)

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
      export_snapgene = FALSE,
      hide_intergenic_labels = hide_intergenic_labels,
      highlight_new = highlight_new,
      clean_headers = clean_headers,
      sample_labels = sample_labels,
      sample_order = sample_order
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
  processed_data <- process_mutation_data_for_viz(mutation_data, sample_order = sample_order)
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
        gc_data <- calculate_gc_content(reference_file, window_size = 10000,
                                        gc_colors = gc_colors)
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
        output_file = circos_file,
        hide_intergenic_labels = hide_intergenic_labels,
        gc_colors = gc_colors
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
      data <- openxlsx::read.xlsx(excel_file, sheet = sheet, check.names = FALSE, sep.names = " ")
      data_list[[sheet]] <- restore_canonical_export_headers(excel_file, sheet, data)
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
#' @param gc_colors Optional named color overrides for skew_positive,
#'   skew_negative, ratio_positive, ratio_negative. Defaults are #00FF00,
#'   #FF0000, #0000FF, and #FFFF00, respectively; zero deviation is positive.
#' @return Data frame with start, end, gc_skew, gc_ratio, gc_skew_minus_average, gc_ratio_minus_average, and colors
#' @keywords internal
calculate_gc_content <- function(reference_file, window_size = 10000, gc_colors = NULL) {

  gc_colors <- .breseq_gc_colors(gc_colors)

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

  # Color each metric by its own deviation from the window mean.
  gc_data$gc_skew_color <- ifelse(gc_data$gc_skew_minus_average >= 0,
    gc_colors[["skew_positive"]], gc_colors[["skew_negative"]])
  gc_data$gc_ratio_color <- ifelse(gc_data$gc_ratio_minus_average >= 0,
    gc_colors[["ratio_positive"]], gc_colors[["ratio_negative"]])

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
#' @param sample_order Optional sample IDs in analysis order. Saved Excel metadata
#'   supplies the order when this argument is NULL.
#' @return List with mutation_list, missing_data, and samples
#' @keywords internal
process_mutation_data_for_viz <- function(mutation_data, sample_order = NULL) {

  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package 'reshape2' required. Install with: install.packages('reshape2')")
  }
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' required. Install with: install.packages('gtools')")
  }

  mutations_df <- mutation_data[["Mutations"]]
  missing_df <- mutation_data[["Missing"]]
  requested_order <- sample_order
  if (is.null(requested_order)) {
    requested_order <- attr(mutations_df, "breseq_sample_order", exact = TRUE)
  }
  if (!is.null(requested_order)) {
    requested_order <- as.character(requested_order)
    if (anyNA(requested_order) || any(!nzchar(requested_order)) ||
        anyDuplicated(requested_order)) {
      stop("sample_order must contain unique, non-empty sample IDs")
    }
  }

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
    value_cols <- c("evidence", "seq id", "mutation", "annotation", "gene", "description")
    all_cols <- colnames(mutations_df)
    tag_cols <- setdiff(all_cols, c(base_cols, value_cols))
    stage_columns <- .breseq_stage_columns(mutations_df, sample_order = requested_order)
    sample_by_column <- stats::setNames(stage_columns$sample, all_cols[stage_columns$index])
    field_by_column <- stats::setNames(stage_columns$field, all_cols[stage_columns$index])

    # Extract sample names from column names
    get_sample_name <- function(col_name) {
      if (col_name %in% value_cols) return("None")
      if (col_name %in% names(sample_by_column)) return(unname(sample_by_column[[col_name]]))
      return("Unknown")
    }

    # Debug: Print column names and detected samples
    message("    Column names: ", paste(utils::head(all_cols, 10), collapse = ", "),
            if(length(all_cols) > 10) "..." else "")
    test_samples <- unique(sapply(tag_cols, get_sample_name))
    message("    Detected samples from columns: ", paste(test_samples, collapse = ", "))

    # Melt to long format
    melted <- reshape2::melt(
      transform(mutations_df, .breseq_source_row = seq_len(nrow(mutations_df))),
      id.vars = c(".breseq_source_row", "position"),
      measure.vars = setdiff(all_cols, "position"),
      variable.name = "variable",
      value.name = "value"
    )

    # Extract sample and variable type
    melted$ALE_version <- sapply(as.character(melted$variable), get_sample_name)
    melted$variable_type <- sapply(as.character(melted$variable), function(v) {
      if (v %in% value_cols) return(v)
      if (v %in% names(field_by_column)) return(unname(field_by_column[[v]]))
      return("other")
    })

    # Filter and cast
    melted <- melted[!is.na(melted$value) & melted$variable_type %in% value_cols, ]

    if (nrow(melted) > 0) {
      mutation_casted <- reshape2::dcast(
        melted,
        .breseq_source_row + position + ALE_version ~ variable_type,
        value.var = "value",
        fun.aggregate = function(x) paste(stats::na.omit(x), collapse = ";")
      )

      # Ensure required columns exist
      for (col in c("mutation", "annotation", "gene", "description")) {
        if (!col %in% names(mutation_casted)) mutation_casted[[col]] <- NA
      }
      source_seq_col <- intersect(c("seq id", "seq_id", "seqid"), names(mutations_df))
      if (!"seq id" %in% names(mutation_casted) && length(source_seq_col)) {
        mutation_casted[["seq id"]] <- mutations_df[[source_seq_col[1]]][mutation_casted$.breseq_source_row]
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
      if (is.null(requested_order)) {
        samples <- samples[gtools::mixedorder(samples)]
      } else {
        omitted <- setdiff(samples, requested_order)
        if (length(omitted)) {
          stop("sample_order omits observed samples: ", paste(omitted, collapse = ", "))
        }
        samples <- requested_order
      }

      # Create sample colors
      n_samples <- length(samples)
      if (n_samples > 0) {
        sample_colors <- setNames(
          RColorBrewer::brewer.pal(max(n_samples, 3), "Set2")[1:n_samples],
          samples
        )

        # Track new mutations per sample
        previous_mutations <- character()

        for (sample in samples) {
          sample_data <- mutation_casted[mutation_casted$ALE_version == sample, ]

          if (nrow(sample_data) > 0) {
            # Mark new mutations
            replicon <- if ("seq id" %in% names(sample_data)) as.character(sample_data[["seq id"]]) else ""
            mutation_keys <- paste(replicon, sample_data$position, sample_data$mutation, sep = "\r")
            sample_data$new <- !mutation_keys %in% previous_mutations

            sample_data$version_color <- sample_colors[sample]

            sample_data$.breseq_source_row <- NULL
            mutation_list[[sample]] <- sample_data
            previous_mutations <- unique(c(previous_mutations, mutation_keys))
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
  stage_columns <- .breseq_stage_columns(missing_df, sample_order = samples)
  sample_by_column <- stats::setNames(stage_columns$sample, all_cols[stage_columns$index])
  field_by_column <- stats::setNames(stage_columns$field, all_cols[stage_columns$index])

  # Get sample names
  get_sample_name <- function(col_name) {
    if (col_name %in% value_cols) return("None")

    if (col_name %in% names(sample_by_column)) return(unname(sample_by_column[[col_name]]))
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
    if (v %in% value_cols) return(v)
    if (v %in% names(field_by_column)) return(unname(field_by_column[[v]]))
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
#' @param gc_data GC content and skew data for inner tracks (optional)
#' @param output_file Output PDF file path
#' @param hide_intergenic_labels Hide labels whose annotation starts with
#'   "intergenic" (default: FALSE). Mutation marks and data are retained.
#' @param gc_colors Optional named overrides for skew_positive, skew_negative,
#'   ratio_positive, ratio_negative. Defaults are #00FF00/#FF0000 for GC skew
#'   and #0000FF/#FFFF00 for GC ratio. Colors are recomputed from metric deviations
#'   at render time; cached color columns are ignored. Zero is positive.
#' @param observed_tracks If TRUE, each mutation ring contains only calls observed
#'   in that sample. Default FALSE preserves cumulative tracks. The new flag
#'   controls labels only, not whether an observed call receives a mark.
#' @param reverse_tracks If TRUE (default), the last sample is the outermost ring.
#'   If FALSE, samples run outside to inside in the supplied order. Sample colors
#'   and the chronological order used for cumulative calls are unchanged.
#' @return Path to created file
#' @keywords internal
create_full_circos_plot <- function(mutation_list,
                                    missing_data = NULL,
                                    samples,
                                    genome_info = NULL,
                                    gff_data = NULL,
                                    eggnog_data = NULL,
                                    gc_data = NULL,
                                    output_file = "circos_map.pdf",
                                    hide_intergenic_labels = FALSE,
                                    gc_colors = NULL,
                                    observed_tracks = FALSE,
                                    reverse_tracks = TRUE) {

  gc_colors <- .breseq_gc_colors(gc_colors)
  if (!is.logical(observed_tracks) || length(observed_tracks) != 1L || is.na(observed_tracks)) {
    stop("observed_tracks must be a single TRUE/FALSE value.")
  }
  if (!is.logical(reverse_tracks) || length(reverse_tracks) != 1L || is.na(reverse_tracks)) {
    stop("reverse_tracks must be a single TRUE/FALSE value.")
  }

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
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(file = output_file, width = 14, height = 14, family = "Arial")
  } else {
    grDevices::pdf(file = output_file, width = 14, height = 14)
  }
  plot_device <- grDevices::dev.cur()
  on.exit({
    if (plot_device %in% grDevices::dev.list()) {
      grDevices::dev.set(plot_device)
      try(circlize::circos.clear(), silent = TRUE)
      grDevices::dev.off(plot_device)
    }
  }, add = TRUE)

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
    df$sample <- rep(sample, nrow(df))
    if (!nrow(df)) {
      df$type <- character()
      df$type_color <- character()
      df$start <- numeric()
      df$end <- numeric()
      if (!"new" %in% names(df)) df$new <- logical()
      if (!"label" %in% names(df)) df$label <- character()
      return(df)
    }
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
  if (is.null(all_mutations)) {
    all_mutations <- data.frame(sample = character(), type = character(),
      type_color = character(), start = numeric(), end = numeric(),
      new = logical(), label = character())
  }

  # Filter new mutations for plotting
  new_mutations <- all_mutations[all_mutations$new == TRUE, ]

  # Add labels for new mutations
  if (nrow(new_mutations) > 0 && any(!is.na(new_mutations$label) & new_mutations$label != "")) {
    label_data <- new_mutations[!is.na(new_mutations$label) & new_mutations$label != "", ]
    label_data <- .breseq_filter_intergenic_labels(label_data, hide_intergenic_labels)
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

  # Ring direction is independent of chronological sample order and colors.
  track_idx <- circlize::get.current.track.index() + 1
  track_samples <- if (reverse_tracks) rev(samples) else samples
  mutation_track_start <- track_idx

  for (i in seq_along(track_samples)) {
    sample <- track_samples[i]
    if (observed_tracks) {
      cum_data <- all_mutations[all_mutations$sample == sample, , drop = FALSE]
    } else {
      # Cumulative: all mutations up to this sample in original order
      orig_idx <- which(samples == sample)
      cumulative_samples <- samples[1:orig_idx]
      cum_data <- all_mutations[all_mutations$sample %in% cumulative_samples, , drop = FALSE]
    }

    bg_color <- adjustcolor(sample_colors[sample], alpha.f = 0.2)

    message("    Track ", i, " (", sample, "): ", nrow(cum_data),
            if (observed_tracks) " observed mutations" else " cumulative mutations")

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
    if (!is.finite(max_skew) || max_skew == 0) max_skew <- 1

    # Force evaluation of colors before closure (closure 문제 방지)
    skew_values <- gc_plot[[skew_col]]
    skew_positions <- (gc_plot$start + gc_plot$end) / 2
    skew_widths <- gc_plot$end - gc_plot$start + 1
    skew_colors <- ifelse(skew_values >= 0,
      gc_colors[["skew_positive"]], gc_colors[["skew_negative"]])

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
          bar_width = skew_widths
        )
      }
    )
    track_idx <- track_idx + 1

    message("    Adding GC ratio track (circos.barplot, GenomeDrawer style)...")

    max_ratio <- max(abs(gc_plot[[ratio_col]]), na.rm = TRUE)
    if (!is.finite(max_ratio) || max_ratio == 0) max_ratio <- 1

    # Force evaluation of colors before closure (closure 문제 방지)
    ratio_values <- gc_plot[[ratio_col]]
    ratio_positions <- (gc_plot$start + gc_plot$end) / 2
    ratio_widths <- gc_plot$end - gc_plot$start + 1
    ratio_colors <- ifelse(ratio_values >= 0,
      gc_colors[["ratio_positive"]], gc_colors[["ratio_negative"]])

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
          bar_width = ratio_widths
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
    labels = track_samples,  # Match outside-to-inside track order
    title = "ALE Version (Track)",
    legend_gp = grid::gpar(fill = adjustcolor(sample_colors[track_samples], alpha.f = 0.3)),
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

  circlize::circos.clear()
  grDevices::dev.off(plot_device)

  return(output_file)
}

# Bundled: gc_colors.R
# Resolve GC track colors independently of cached data-frame color columns.
# The sign is the deviation from the window mean, with zero using positive.
.breseq_gc_colors <- function(gc_colors = NULL) {
  palette <- c(skew_positive = "#00FF00", skew_negative = "#FF0000",
               ratio_positive = "#0000FF", ratio_negative = "#FFFF00")
  if (is.null(gc_colors)) return(palette)
  if (!is.character(gc_colors) || !length(gc_colors) || is.null(names(gc_colors)) ||
      anyNA(gc_colors) || anyNA(names(gc_colors)) || any(!nzchar(names(gc_colors))) ||
      anyDuplicated(names(gc_colors)) || any(!names(gc_colors) %in% names(palette))) {
    stop("gc_colors must be a named character vector with distinct names from: ",
         paste(names(palette), collapse = ", "), ".")
  }
  valid <- tryCatch({
    grDevices::col2rgb(unname(gc_colors))
    all(nzchar(gc_colors))
  }, error = function(e) FALSE)
  if (!valid) stop("gc_colors must contain valid R color names or hexadecimal colors.")
  palette[names(gc_colors)] <- gc_colors
  palette
}

# Bundled: get_screen_size.R
get_screen_size <- function(default_width = 12, default_height = 7, dpi = 96) {
  screen_width_in <- default_width
  screen_height_in <- default_height
  tryCatch({
    # Run system command to get display info
    display_info <- system("system_profiler SPDisplaysDataType | grep Resolution", intern = TRUE)
    if (length(display_info) > 0) {
      # Example output: "Resolution: 2560 x 1600 Retina"
      resolution <- strsplit(display_info[1], "x")[[1]]
      width_px <- as.numeric(gsub("[^0-9]", "", resolution[1]))
      height_px <- as.numeric(gsub("[^0-9]", "", resolution[2]))
      # Convert pixels to inches
      screen_width_in <- width_px / dpi
      screen_height_in <- height_px / dpi
      message("Detected screen resolution: ", width_px, "x", height_px, " pixels")
      message("Converted to inches (DPI=", dpi, "): ", round(screen_width_in, 2), "x", round(screen_height_in, 2))
      # Cap at reasonable size (e.g., 15 inches)
      screen_width_in <- min(screen_width_in, 15)
      screen_height_in <- min(screen_height_in, 15 * (height_px / width_px)) # Maintain aspect ratio
    }
  }, error = function(e) {
    message("Failed to detect screen size: ", e$message)
    message("Using default PDF size: ", default_width, "x", default_height, " inches")
  })
  return(list(width = screen_width_in, height = screen_height_in))
}

# Bundled: normalize_position.R
#' Normalize Genome Position for Circular Sequences
#'
#' Handles position normalization for circular genomes (plasmids, bacterial chromosomes).
#' Ensures positions wrap around correctly and handles edge cases near origin.
#'
#' @param position Numeric position(s) to normalize
#' @param genome_length Total length of the circular genome
#' @param origin Position of the origin (default = 1)
#' @return Normalized position(s) within 1:genome_length
#' @export
#' @examples
#' # Position beyond genome length wraps around
#' normalize_position(5100, genome_length = 5000)  # Returns 100
#'
#' # Negative positions wrap from end
#' normalize_position(-50, genome_length = 5000)   # Returns 4950
normalize_position <- function(position, genome_length, origin = 1) {

  if (is.na(position) || is.null(genome_length) || genome_length <= 0) {
    return(NA_integer_)
  }

  # Handle numeric conversion
  pos <- as.numeric(gsub(",", "", as.character(position)))

  if (is.na(pos)) return(NA_integer_)

  # Normalize to 1-based coordinates
  normalized <- ((pos - origin) %% genome_length) + 1

  # Ensure positive
  if (normalized <= 0) {
    normalized <- normalized + genome_length
  }

  return(as.integer(normalized))
}

#' Normalize Position Range for Circular Sequences
#'
#' Normalizes a start-end range for circular genomes. Handles cases where
#' the feature spans the origin.
#'
#' @param start Start position
#' @param end End position
#' @param genome_length Total length of the circular genome
#' @return List with normalized start, end, and spans_origin flag
#' @export
#' @examples
#' # Normal range
#' normalize_range(100, 500, genome_length = 5000)
#'
#' # Range spanning origin
#' normalize_range(4900, 100, genome_length = 5000)
normalize_range <- function(start, end, genome_length) {

  start_norm <- normalize_position(start, genome_length)
  end_norm <- normalize_position(end, genome_length)

  # Check if range spans the origin
  spans_origin <- FALSE
  if (!is.na(start_norm) && !is.na(end_norm)) {
    if (start_norm > end_norm) {
      spans_origin <- TRUE
    }
  }

  list(
    start = start_norm,
    end = end_norm,
    spans_origin = spans_origin,
    # For features spanning origin, provide both segments
    segments = if (spans_origin) {
      list(
        c(start_norm, genome_length),
        c(1, end_norm)
      )
    } else {
      list(c(start_norm, end_norm))
    }
  )
}

#' Parse Position String from breseq Output
#'
#' Parses various position formats from breseq HTML output.
#' Handles formats like: "1234", "1,234", "1234-5678", "1234–5678" (en-dash)
#' and junction endpoints such as "= 1234" or "1234 =". Orientation markers are
#' ignored only for the returned coordinates; the reported input is not changed.
#'
#' @param position_str Position string to parse
#' @return List with start and end positions, or single position if not a range
#' @export
parse_position <- function(position_str) {

  if (is.null(position_str) || !length(position_str) || is.na(position_str)) {
    return(list(start = NA_integer_, end = NA_integer_))
  }

  # Convert to string and clean
  pos_str <- as.character(position_str)
  pos_str <- gsub(",", "", pos_str)           # Remove commas
  pos_str <- gsub("–", "-", pos_str)          # Replace en-dash with hyphen
  pos_str <- gsub("—", "-", pos_str)          # Replace em-dash with hyphen
  pos_str <- gsub("\u00a0", " ", pos_str, fixed = TRUE)
  pos_str <- trimws(pos_str)
  pos_str <- gsub("^[[:space:]=<>]+|[[:space:]=<>]+$", "", pos_str)
  pos_str <- gsub("-$", "", pos_str)          # Remove trailing dash
  pos_str <- trimws(pos_str)

  # Check for range format
  if (grepl("^[0-9]+[[:space:]]*-[[:space:]]*[0-9]+$", pos_str)) {
    parts <- strsplit(pos_str, "-")[[1]]
    if (length(parts) == 2) {
      start <- as.integer(parts[1])
      end <- as.integer(parts[2])
      return(list(start = start, end = end))
    }
  }

  # Single position
  pos <- if (grepl("^[0-9]+$", pos_str)) as.integer(pos_str) else NA_integer_
  return(list(start = pos, end = pos))
}

#' Calculate Distance Between Positions on Circular Genome
#'
#' Calculates the shortest distance between two positions on a circular genome.
#'
#' @param pos1 First position
#' @param pos2 Second position
#' @param genome_length Total genome length
#' @return Shortest distance (always positive)
#' @export
circular_distance <- function(pos1, pos2, genome_length) {

  if (any(is.na(c(pos1, pos2, genome_length)))) return(NA_real_)

  direct <- abs(pos2 - pos1)
  wraparound <- genome_length - direct

  min(direct, wraparound)
}

#' Merge Overlapping Mutations on Circular Genome
#'
#' Identifies and optionally merges mutations that are at the same position
#' or overlapping, accounting for circular genome topology.
#'
#' @param df Data frame with position column
#' @param position_col Name of position column
#' @param genome_length Optional genome length for circular handling
#' @param merge_distance Maximum distance to consider as same position (default = 0)
#' @param seq_id_col Optional column for sequence ID (for multi-replicon genomes)
#' @return Data frame with added 'position_group' column for identifying duplicates
#' @export
identify_position_duplicates <- function(df,
                                         position_col = "position",
                                         genome_length = NULL,
                                         merge_distance = 0,
                                         seq_id_col = NULL) {

  if (!position_col %in% names(df)) {
    warning("Position column '", position_col, "' not found.")
    df$position_group <- seq_len(nrow(df))
    return(df)
  }

  # Parse positions
  parsed <- lapply(df[[position_col]], parse_position)
  df$parsed_start <- sapply(parsed, `[[`, "start")

  # If sequence ID column exists, group by it first
  if (!is.null(seq_id_col) && seq_id_col %in% names(df)) {
    df$seq_key <- paste(df[[seq_id_col]], df$parsed_start, sep = "_")
  } else {
    df$seq_key <- as.character(df$parsed_start)
  }

  # Assign position groups
  if (is.null(genome_length) || merge_distance == 0) {
    # Simple exact matching
    df$position_group <- as.integer(factor(df$seq_key))
  } else {
    # Consider circular distance for merging
    positions <- df$parsed_start
    n <- length(positions)
    groups <- rep(0, n)
    current_group <- 0

    for (i in seq_len(n)) {
      if (groups[i] == 0) {
        current_group <- current_group + 1
        groups[i] <- current_group

        for (j in seq_len(n)[-i]) {
          if (groups[j] == 0) {
            dist <- circular_distance(positions[i], positions[j], genome_length)
            if (!is.na(dist) && dist <= merge_distance) {
              groups[j] <- current_group
            }
          }
        }
      }
    }
    df$position_group <- groups
  }

  # Clean up temporary columns
  df$parsed_start <- NULL
  df$seq_key <- NULL

  # Add duplicate flag
  dup_counts <- table(df$position_group)
  df$is_duplicate <- df$position_group %in% names(dup_counts[dup_counts > 1])

  return(df)
}

#' Normalize Positions in Data Frame
#'
#' Applies position normalization to an entire data frame, handling
#' multiple replicons (chromosome + plasmids) if present.
#' Can auto-detect genome lengths from reference file.
#'
#' @param df Data frame with position data
#' @param position_col Name of position column(s)
#' @param genome_lengths Named vector of genome lengths (names = seq_id), OR
#' @param reference_file Path to reference file (.gbk, .fasta, .dna) to auto-parse lengths
#' @param seq_id_col Column containing sequence ID (auto-detected if NULL)
#' @return Data frame with normalized positions
#' @export
#' @examples
#' \dontrun{
#' # Auto-detect from GenBank (recommended)
#' df <- normalize_positions_df(df, reference_file = "reference.gbk")
#'
#' # Or provide lengths manually
#' df <- normalize_positions_df(df, genome_lengths = c(chromosome = 4600000, pBR322 = 4361))
#' }
normalize_positions_df <- function(df,
                                   position_col = "position",
                                   genome_lengths = NULL,
                                   reference_file = NULL,
                                   seq_id_col = NULL) {

  # Auto-parse genome lengths from reference file
  if (is.null(genome_lengths) && !is.null(reference_file)) {
    genome_lengths <- get_genome_lengths(reference_file)
  }

  if (is.null(genome_lengths)) {
    message("No genome lengths provided. Positions not normalized for circularity.")
    return(df)
  }

  # Auto-detect seq_id column
  if (is.null(seq_id_col)) {
    seq_id_col <- find_seq_id_column(df)
  }

  # Normalize each position based on its replicon
  for (col in position_col) {
    if (!col %in% names(df)) next

    new_col <- paste0(col, "_normalized")

    if (!is.null(seq_id_col) && seq_id_col %in% names(df)) {
      df[[new_col]] <- mapply(function(pos, sid) {
        glen <- genome_lengths[sid]
        if (is.na(glen)) return(pos)
        normalize_position(pos, glen)
      }, df[[col]], df[[seq_id_col]])
    } else if (length(genome_lengths) == 1) {
      df[[new_col]] <- sapply(df[[col]], normalize_position,
                               genome_length = genome_lengths[1])
    } else {
      warning("Multiple replicons detected but no seq_id column found. ",
              "Using first genome length for all positions.")
      df[[new_col]] <- sapply(df[[col]], normalize_position,
                               genome_length = genome_lengths[1])
    }
  }

  return(df)
}

# Bundled: parse_reference.R
#' Parse Reference Genome Information
#'
#' Extracts replicon names, lengths, and topology (linear/circular) from
#' GenBank, FASTA, or SnapGene .dna files. Supports multi-replicon genomes
#' (chromosome + multiple plasmids).
#'
#' @param reference_file Path to reference file (.gbk, .gb, .gbff, .fasta, .fa, .fna, .dna)
#' @return Data frame with columns: seq_id, length, topology, description
#' @export
#' @examples
#' \dontrun{
#' # Get info from GenBank file
#' info <- parse_reference_info("reference.gbk")
#' # Returns:
#' #   seq_id    length topology description
#' # 1 chromosome 4600000 circular E. coli K-12
#' # 2 pBR322       4361 circular Cloning vector
#' # 3 pUC19        2686 circular Cloning vector
#'
#' # Use for normalization
#' genome_lengths <- setNames(info$length, info$seq_id)
#' }
parse_reference_info <- function(reference_file) {

  if (!file.exists(reference_file)) {
    stop("Reference file not found: ", reference_file)
  }

  ext <- tolower(tools::file_ext(reference_file))

  result <- switch(ext,
    "gbk" = , "gb" = , "gbff" = , "genbank" = parse_genbank_info(reference_file),
    "fasta" = , "fa" = , "fna" = , "faa" = parse_fasta_info(reference_file),
    "dna" = parse_dna_info(reference_file),
    stop("Unsupported file format: ", ext,
         "\nSupported: .gbk, .gb, .gbff, .fasta, .fa, .fna, .dna")
  )

  message("Parsed ", nrow(result), " replicon(s) from ", basename(reference_file))
  for (i in seq_len(nrow(result))) {
    message("  ", result$seq_id[i], ": ",
            format(result$length[i], big.mark = ","), " bp (",
            result$topology[i], ")")
  }

  return(result)
}

#' Parse GenBank File for Replicon Information
#'
#' Reads GenBank file(s) and extracts LOCUS information for each replicon.
#'
#' @param gbk_file Path to GenBank file
#' @return Data frame with seq_id, length, topology, description
#' @export
parse_genbank_info <- function(gbk_file) {

  lines <- readLines(gbk_file, warn = FALSE)

  # Find all LOCUS lines
  locus_idx <- grep("^LOCUS", lines)

  if (length(locus_idx) == 0) {
    stop("No LOCUS line found in GenBank file")
  }

  results <- lapply(locus_idx, function(idx) {
    locus_line <- lines[idx]

    # Parse LOCUS line
    # Format: LOCUS       name    length bp    mol_type    topology    division date
    # Example: LOCUS       NC_000913  4641652 bp    DNA     circular BCT 17-SEP-2019

    # Extract sequence ID (first field after LOCUS)
    parts <- strsplit(trimws(locus_line), "\\s+")[[1]]
    seq_id <- parts[2]

    # Extract length (number followed by "bp")
    bp_match <- regmatches(locus_line, regexpr("\\d+\\s*bp", locus_line))
    length <- if (length(bp_match) > 0) {
      as.integer(gsub("\\s*bp", "", bp_match))
    } else NA_integer_

    # Extract topology (circular or linear)
    topology <- if (grepl("circular", locus_line, ignore.case = TRUE)) {
      "circular"
    } else if (grepl("linear", locus_line, ignore.case = TRUE)) {
      "linear"
    } else {
      "unknown"
    }

    # Find DEFINITION line for description
    def_idx <- idx + which(grepl("^DEFINITION", lines[(idx+1):min(idx+10, length(lines))]))[1]
    description <- if (!is.na(def_idx) && def_idx <= length(lines)) {
      sub("^DEFINITION\\s+", "", lines[def_idx])
    } else { "" }

    data.frame(
      seq_id = seq_id,
      length = length,
      topology = topology,
      description = description,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Parse FASTA File for Sequence Information
#'
#' Reads FASTA file and calculates sequence lengths.
#' Topology is inferred from header if possible, otherwise assumed circular for bacteria.
#'
#' @param fasta_file Path to FASTA file
#' @param default_topology Default topology if not specified ("circular" or "linear")
#' @return Data frame with seq_id, length, topology, description
#' @export
parse_fasta_info <- function(fasta_file, default_topology = "circular") {

  lines <- readLines(fasta_file, warn = FALSE)

  # Find header lines
  header_idx <- grep("^>", lines)

  if (length(header_idx) == 0) {
    stop("No FASTA headers found")
  }

  # Add end position for last sequence
  seq_ends <- c(header_idx[-1] - 1, length(lines))

  results <- lapply(seq_along(header_idx), function(i) {
    header <- lines[header_idx[i]]
    seq_lines <- lines[(header_idx[i] + 1):seq_ends[i]]

    # Parse header: >seq_id description
    header_clean <- sub("^>", "", header)
    parts <- strsplit(header_clean, "\\s+", perl = TRUE)[[1]]
    seq_id <- parts[1]
    description <- if (length(parts) > 1) paste(parts[-1], collapse = " ") else ""

    # Calculate sequence length
    seq_length <- sum(nchar(gsub("\\s", "", seq_lines)))

    # Infer topology from description
    topology <- if (grepl("circular|plasmid", description, ignore.case = TRUE)) {
      "circular"
    } else if (grepl("linear|chromosome", description, ignore.case = TRUE)) {
      # Note: bacterial chromosomes are usually circular
      if (grepl("plasmid", description, ignore.case = TRUE)) "circular" else default_topology
    } else {
      default_topology
    }

    data.frame(
      seq_id = seq_id,
      length = seq_length,
      topology = topology,
      description = description,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Parse SnapGene .dna File for Sequence Information
#'
#' Converts .dna to GenBank using SnapGene CLI, then parses the result.
#' Falls back to binary parsing if SnapGene is not available.
#'
#' @param dna_file Path to .dna file
#' @return Data frame with seq_id, length, topology, description
#' @export
parse_dna_info <- function(dna_file) {

  # Try to convert to GenBank first
  if (has_snapgene()) {
    temp_gbk <- tempfile(fileext = ".gbk")
    on.exit(unlink(temp_gbk), add = TRUE)

    tryCatch({
      convert_dna_to_genbank(dna_file, temp_gbk)
      return(parse_genbank_info(temp_gbk))
    }, error = function(e) {
      message("SnapGene conversion failed, trying binary parse: ", e$message)
    })
  }

  # Fallback: parse .dna binary format directly
  # SnapGene .dna format is documented at:
  # https://www.snapgene.com/resources/file-format
  parse_dna_binary(dna_file)
}

#' Parse SnapGene .dna Binary Format
#'
#' Directly parses the SnapGene binary format without needing SnapGene installed.
#'
#' @param dna_file Path to .dna file
#' @return Data frame with seq_id, length, topology, description
#' @keywords internal
parse_dna_binary <- function(dna_file) {

  # Read binary file
  con <- file(dna_file, "rb")
  on.exit(close(con), add = TRUE)

  # Read header
  header <- readBin(con, "raw", n = 1)

  # Check magic byte (0x09 for SnapGene)
  if (as.integer(header) != 9) {
    warning("File may not be a valid SnapGene .dna file")
  }

  # Read file structure
  # SnapGene uses a chunk-based format

  seq_length <- NA_integer_
  topology <- "circular"  # Default assumption
  seq_id <- tools::file_path_sans_ext(basename(dna_file))
  description <- ""

  # Try to read sequence length from file
  tryCatch({
    # Seek to find DNA sequence chunk
    seek(con, 0)
    raw_data <- readBin(con, "raw", n = file.info(dna_file)$size)

    # Look for sequence data (usually after header chunks)
    # The DNA sequence is stored as ASCII in the file
    raw_text <- rawToChar(raw_data[raw_data >= 0x20 & raw_data <= 0x7E])

    # Try to find sequence (ATCG only regions)
    dna_matches <- gregexpr("[ATCGatcg]{100,}", raw_text)[[1]]
    if (dna_matches[1] > 0) {
      # Find the longest match (likely the main sequence)
      match_lengths <- attr(dna_matches, "match.length")
      best_match <- which.max(match_lengths)
      seq_length <- match_lengths[best_match]
    }

    # Check for circular indicator
    if (grepl("circular|\\bcircular\\b", raw_text, ignore.case = TRUE)) {
      topology <- "circular"
    } else if (grepl("linear", raw_text, ignore.case = TRUE)) {
      topology <- "linear"
    }

  }, error = function(e) {
    warning("Could not parse .dna file structure: ", e$message)
  })

  data.frame(
    seq_id = seq_id,
    length = seq_length,
    topology = topology,
    description = description,
    stringsAsFactors = FALSE
  )
}

#' Get Genome Lengths as Named Vector
#'
#' Convenience function to get genome lengths in a format ready for use
#' with normalization functions.
#'
#' @param reference_file Path to reference file
#' @return Named numeric vector (names = seq_id, values = length)
#' @export
#' @examples
#' \dontrun{
#' genome_lengths <- get_genome_lengths("reference.gbk")
#' # chromosome      pBR322       pUC19
#' #    4600000        4361        2686
#' }
get_genome_lengths <- function(reference_file) {
  info <- parse_reference_info(reference_file)
  setNames(info$length, info$seq_id)
}

#' Get Replicon Topology as Named Vector
#'
#' @param reference_file Path to reference file
#' @return Named character vector (names = seq_id, values = topology)
#' @export
get_replicon_topology <- function(reference_file) {
  info <- parse_reference_info(reference_file)
  setNames(info$topology, info$seq_id)
}

#' Check if Replicon is Circular
#'
#' @param seq_id Sequence ID to check
#' @param reference_file Path to reference file
#' @return Logical
#' @export
is_circular <- function(seq_id, reference_file) {
  topology <- get_replicon_topology(reference_file)
  if (seq_id %in% names(topology)) {
    return(topology[seq_id] == "circular")
  }
  # Default: assume circular for bacterial genomes
  TRUE
}

#' Auto-detect seq_id Column in breseq Data
#'
#' breseq output includes a "seq id" column that identifies which replicon
#' (chromosome or plasmid) each mutation is on.
#'
#' @param df Data frame from breseq
#' @return Name of the seq_id column, or NULL if not found
#' @keywords internal
find_seq_id_column <- function(df) {
  candidates <- c("seq id", "seq_id", "seqid", "seq.id",
                  "replicon", "contig", "chromosome")

  for (col in candidates) {
    if (col %in% names(df)) return(col)
  }

  # Try case-insensitive match
  col_lower <- tolower(names(df))
  for (cand in candidates) {
    idx <- which(col_lower == tolower(cand))
    if (length(idx) > 0) return(names(df)[idx[1]])
  }

  NULL
}

# Bundled: plot_mutations.R
#' Plot Mutation Landscape
#'
#' Creates a visual representation of mutations across the genome,
#' showing mutation types, positions, and evolution stages.
#'
#' @param mutation_data Data frame with mutation data or path to Excel file
#' @param genome_length Total genome length (for circular plot option)
#' @param plot_type Type of plot: "linear", "circular", or "timeline"
#' @param color_by Color mutations by "type" or "stage"
#' @param title Plot title
#' @param show_labels Show mutation labels (TRUE/FALSE)
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' plot_mutations(comparison_results$`Mutation predictions`,
#'                genome_length = 4600000,
#'                plot_type = "linear")
#' }
plot_mutations <- function(mutation_data,
                           genome_length = NULL,
                           plot_type = c("linear", "circular", "timeline"),
                           color_by = c("type", "stage"),
                           title = "Mutation Landscape",
                           show_labels = FALSE) {

  plot_type <- match.arg(plot_type)
  color_by <- match.arg(color_by)

  # Check for ggplot2
 if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required. Install with: install.packages('ggplot2')")
  }

  # Load data if path provided
  if (is.character(mutation_data) && file.exists(mutation_data)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' required for Excel files")
    }
    mutation_data <- openxlsx::read.xlsx(mutation_data, sheet = "Mutations")
  }

  df <- as.data.frame(mutation_data)

  # Find position column
  pos_col <- if ("position" %in% names(df)) "position" else
             if ("start" %in% names(df)) "start" else
             grep("position", names(df), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(pos_col) || !pos_col %in% names(df)) {
    stop("No position column found in data")
  }

  # Parse positions
  df$pos_numeric <- sapply(df[[pos_col]], function(x) {
    parsed <- parse_position(x)
    parsed$start
  })

  # Add mutation type if not present
  if (!"mutation_type" %in% names(df)) {
    mut_col <- grep("^mutation", names(df), value = TRUE)[1]
    if (!is.na(mut_col)) {
      df$mutation_type <- classify_mutations(df[[mut_col]])
    } else {
      df$mutation_type <- "unknown"
    }
  }

  # Add colors
  df$color <- sapply(df$mutation_type, get_mutation_color)

  # Detect stages for timeline plot
  stage_cols <- grep("\\.", names(df), value = TRUE)
  if (length(stage_cols) > 0) {
    stage_names <- unique(gsub(".*\\.([^/]+)/.*", "\\1", stage_cols))
    stage_names <- stage_names[!grepl("^Unnamed|^position|^mutation|^annotation", stage_names)]

    # Find first detection stage for each mutation
    df$first_stage <- sapply(seq_len(nrow(df)), function(i) {
      for (stage in stage_names) {
        cols <- grep(stage, names(df), value = TRUE)
        if (any(!is.na(df[i, cols]))) {
          return(stage)
        }
      }
      return(NA_character_)
    })

    df$stage_order <- match(df$first_stage, stage_names)
  } else {
    df$first_stage <- NA
    df$stage_order <- 1
    stage_names <- character(0)
  }

  # Remove NA positions
  df <- df[!is.na(df$pos_numeric), ]

  if (nrow(df) == 0) {
    message("No valid positions to plot")
    return(invisible(NULL))
  }

  # Create plot based on type
  if (plot_type == "linear") {
    p <- create_linear_plot(df, pos_col, color_by, title, show_labels, genome_length)
  } else if (plot_type == "circular") {
    p <- create_circular_plot(df, pos_col, color_by, title, genome_length)
  } else if (plot_type == "timeline") {
    p <- create_timeline_plot(df, pos_col, color_by, title, stage_names)
  }

  return(p)
}

#' Create Linear Mutation Plot
#' @keywords internal
create_linear_plot <- function(df, pos_col, color_by, title, show_labels, genome_length) {

  # Set up color mapping
  if (color_by == "type") {
    color_var <- "mutation_type"
    color_values <- c(
      "SNP" = "#3498db", "insertion" = "#2ecc71", "deletion" = "#e74c3c",
      "large_insertion" = "#27ae60", "large_deletion" = "#c0392b",
      "substitution" = "#9b59b6", "amplification" = "#f39c12",
      "inversion" = "#1abc9c", "mobile_element" = "#e67e22",
      "missing_coverage" = "#95a5a6", "indel" = "#f1c40f", "unknown" = "#bdc3c7"
    )
  } else {
    color_var <- "first_stage"
    n_stages <- length(unique(df$first_stage[!is.na(df$first_stage)]))
    color_values <- NULL  # Use default colors
  }

  # Create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pos_numeric, y = 1)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = pos_numeric, xend = pos_numeric, y = 0.8, yend = 1.2,
                   color = .data[[color_var]]),
      linewidth = 1.5, alpha = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_var]]),
      size = 3, alpha = 0.8
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = title,
      x = "Genome Position (bp)",
      y = "",
      color = if (color_by == "type") "Mutation Type" else "Evolution Stage"
    ) +
    ggplot2::ylim(0.5, 1.5)

  if (color_by == "type" && !is.null(color_values)) {
    p <- p + ggplot2::scale_color_manual(values = color_values, na.value = "#bdc3c7")
  }

  if (!is.null(genome_length)) {
    p <- p + ggplot2::xlim(0, genome_length)
  }

  if (show_labels && "annotation" %in% names(df)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        ggplot2::aes(label = annotation),
        size = 2.5, max.overlaps = 20
      )
    }
  }

  return(p)
}

#' Create Circular Genome Plot
#' @keywords internal
create_circular_plot <- function(df, pos_col, color_by, title, genome_length) {

  if (is.null(genome_length)) {
    genome_length <- max(df$pos_numeric, na.rm = TRUE) * 1.1
  }

  # Convert positions to angles (radians)
  df$angle <- (df$pos_numeric / genome_length) * 2 * pi - pi/2
  df$x <- cos(df$angle)
  df$y <- sin(df$angle)

  # Outer ring coordinates
  df$x_outer <- 1.15 * cos(df$angle)
  df$y_outer <- 1.15 * sin(df$angle)

  if (color_by == "type") {
    color_var <- "mutation_type"
    color_values <- c(
      "SNP" = "#3498db", "insertion" = "#2ecc71", "deletion" = "#e74c3c",
      "large_insertion" = "#27ae60", "large_deletion" = "#c0392b",
      "substitution" = "#9b59b6", "amplification" = "#f39c12",
      "inversion" = "#1abc9c", "mobile_element" = "#e67e22",
      "missing_coverage" = "#95a5a6", "indel" = "#f1c40f", "unknown" = "#bdc3c7"
    )
  } else {
    color_var <- "first_stage"
    color_values <- NULL
  }

  # Create circular reference line
  circle_df <- data.frame(
    angle = seq(0, 2*pi, length.out = 100)
  )
  circle_df$x <- cos(circle_df$angle)
  circle_df$y <- sin(circle_df$angle)

  p <- ggplot2::ggplot() +
    # Reference circle
    ggplot2::geom_path(data = circle_df, ggplot2::aes(x = x, y = y),
                       color = "gray70", linewidth = 2) +
    # Mutation markers
    ggplot2::geom_segment(
      data = df,
      ggplot2::aes(x = x, y = y, xend = x_outer, yend = y_outer,
                   color = .data[[color_var]]),
      linewidth = 1.5, alpha = 0.8
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = x_outer, y = y_outer, color = .data[[color_var]]),
      size = 3
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(
      title = title,
      color = if (color_by == "type") "Mutation Type" else "Evolution Stage"
    )

  if (color_by == "type" && !is.null(color_values)) {
    p <- p + ggplot2::scale_color_manual(values = color_values, na.value = "#bdc3c7")
  }

  # Add position labels
  label_positions <- c(0, 0.25, 0.5, 0.75) * genome_length
  label_angles <- (label_positions / genome_length) * 2 * pi - pi/2
  label_df <- data.frame(
    x = 1.3 * cos(label_angles),
    y = 1.3 * sin(label_angles),
    label = format(label_positions, big.mark = ",", scientific = FALSE)
  )

  p <- p + ggplot2::geom_text(data = label_df, ggplot2::aes(x = x, y = y, label = label),
                              size = 3, color = "gray40")

  return(p)
}

#' Create Timeline Plot Showing Mutation Accumulation
#' @keywords internal
create_timeline_plot <- function(df, pos_col, color_by, title, stage_names) {

  if (length(stage_names) == 0) {
    message("No stage information available for timeline plot")
    return(create_linear_plot(df, pos_col, color_by, title, FALSE, NULL))
  }

  # Ensure stage order factor
  df$first_stage <- factor(df$first_stage, levels = stage_names)

  if (color_by == "type") {
    color_var <- "mutation_type"
    color_values <- c(
      "SNP" = "#3498db", "insertion" = "#2ecc71", "deletion" = "#e74c3c",
      "large_insertion" = "#27ae60", "large_deletion" = "#c0392b",
      "substitution" = "#9b59b6", "amplification" = "#f39c12",
      "inversion" = "#1abc9c", "mobile_element" = "#e67e22",
      "missing_coverage" = "#95a5a6", "indel" = "#f1c40f", "unknown" = "#bdc3c7"
    )
  } else {
    color_var <- "first_stage"
    color_values <- NULL
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = first_stage, y = pos_numeric)) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data[[color_var]]),
      width = 0.2, height = 0, size = 3, alpha = 0.7
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = title,
      x = "Evolution Stage",
      y = "Genome Position (bp)",
      color = if (color_by == "type") "Mutation Type" else "Stage"
    )

  if (color_by == "type" && !is.null(color_values)) {
    p <- p + ggplot2::scale_color_manual(values = color_values, na.value = "#bdc3c7")
  }

  return(p)
}

#' Plot Mutation Type Summary
#'
#' Creates a bar chart showing the count of each mutation type.
#'
#' @param mutation_data Data frame with mutation data
#' @param title Plot title
#' @return ggplot object
#' @export
plot_mutation_summary <- function(mutation_data, title = "Mutation Type Distribution") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  df <- as.data.frame(mutation_data)

  # Add mutation type if needed
  if (!"mutation_type" %in% names(df)) {
    mut_col <- grep("^mutation", names(df), value = TRUE)[1]
    if (!is.na(mut_col)) {
      df$mutation_type <- classify_mutations(df[[mut_col]])
    } else {
      df$mutation_type <- "unknown"
    }
  }

  # Count by type
  type_counts <- as.data.frame(table(df$mutation_type))
  names(type_counts) <- c("type", "count")

  color_values <- c(
    "SNP" = "#3498db", "insertion" = "#2ecc71", "deletion" = "#e74c3c",
    "large_insertion" = "#27ae60", "large_deletion" = "#c0392b",
    "substitution" = "#9b59b6", "amplification" = "#f39c12",
    "inversion" = "#1abc9c", "mobile_element" = "#e67e22",
    "missing_coverage" = "#95a5a6", "indel" = "#f1c40f", "unknown" = "#bdc3c7"
  )

  p <- ggplot2::ggplot(type_counts, ggplot2::aes(x = reorder(type, -count), y = count, fill = type)) +
    ggplot2::geom_col(alpha = 0.8) +
    ggplot2::scale_fill_manual(values = color_values, na.value = "#bdc3c7") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = title,
      x = "Mutation Type",
      y = "Count"
    ) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.5, size = 3)

  return(p)
}

#' Plot Mutation Comparison Heatmap
#'
#' Creates a heatmap showing mutation presence/absence across evolution stages.
#'
#' @param mutation_data Data frame with mutation comparison data
#' @param title Plot title
#' @return ggplot object
#' @export
plot_mutation_heatmap <- function(mutation_data, title = "Mutation Presence Across Stages") {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required")
  }

  df <- as.data.frame(mutation_data)

  # Find stage columns
  stage_cols <- grep("\\.", names(df), value = TRUE)
  stage_names <- unique(gsub(".*\\.([^/]+)/.*", "\\1", stage_cols))
  stage_names <- stage_names[!grepl("^Unnamed|^position|^mutation|^annotation", stage_names)]

  if (length(stage_names) == 0) {
    message("No stage information found for heatmap")
    return(invisible(NULL))
  }

  # Find position column
  pos_col <- if ("position" %in% names(df)) "position" else "start"

  # Create presence/absence matrix
  heat_data <- data.frame(position = df[[pos_col]])

  for (stage in stage_names) {
    stage_col <- grep(paste0("mutation\\.", stage), names(df), value = TRUE)[1]
    if (!is.na(stage_col) && stage_col %in% names(df)) {
      heat_data[[stage]] <- as.integer(!is.na(df[[stage_col]]))
    }
  }

  # Reshape for ggplot
  heat_long <- stats::reshape(
    heat_data,
    direction = "long",
    varying = stage_names,
    v.names = "present",
    times = stage_names,
    timevar = "stage"
  )

  heat_long$position <- factor(heat_long$position)
  heat_long$stage <- factor(heat_long$stage, levels = stage_names)

  p <- ggplot2::ggplot(heat_long, ggplot2::aes(x = stage, y = position, fill = factor(present))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c("0" = "#f5f5f5", "1" = "#e74c3c"),
      labels = c("Absent", "Present"),
      name = "Mutation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 6)
    ) +
    ggplot2::labs(title = title, x = "Evolution Stage", y = "Position")

  return(p)
}

# Bundled: run_pipeline.R
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
#' @param hide_intergenic_labels Hide only intergenic Circos labels, keeping all mutation marks.
#' @param highlight_new Shade reported Excel calls by stage and bold new or changed calls; default FALSE.
#' @param clean_headers Export readable Excel headers and hidden canonical-column metadata.
#' @param sample_labels Named sample-ID to display-label map, used only with clean_headers=TRUE.
#' @param sample_order Optional exact chronological order of all detected sample IDs.
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
                                color_by = c("type", "stage"),
                                hide_intergenic_labels = FALSE,
                                highlight_new = FALSE,
                                clean_headers = FALSE,
                                sample_labels = NULL,
                                sample_order = NULL) {

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

    # Preserve the complete identifier, including suffixes such as 3w.
    if (!identical(tolower(filename), "index")) {
      sample_name <- sub("_index$", "", sub("^index_", "", filename, ignore.case = TRUE), ignore.case = TRUE)
    } else {
      dir_name <- basename(dirname(html_file))
      sample_name <- if (dir_name == "output") basename(dirname(dirname(html_file))) else dir_name
    }
    if (sample_name %in% names(extracted_tables)) stop("Duplicate sample ID: ", sample_name)

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
  if (!is.null(sample_order)) {
    if (!is.character(sample_order) || anyNA(sample_order) || anyDuplicated(sample_order) ||
        !setequal(sample_order, sample_names)) stop("sample_order must contain every detected sample exactly once.")
    sample_names <- sample_order
  } else if (requireNamespace("gtools", quietly = TRUE)) {
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
  if (interactive() && is.null(sample_order)) {
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
      # A shared junction endpoint is not a duplicate paired event.
      if (.breseq_table_kind(name) == "junction") next
      pos_col <- if ("position" %in% names(df)) "position" else "start"
      row_seqids <- .breseq_wide_seq(df, .breseq_stage_columns(df))
      seq_id_col <- NULL
      if (any(nzchar(row_seqids))) {
        seq_id_col <- ".breseq_normalization_seq"
        df[[seq_id_col]] <- row_seqids
      }
      comparison_keys <- attr(df, "breseq_keys")

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
      if (!is.null(seq_id_col)) merged_results[[name]][[seq_id_col]] <- NULL
      attr(merged_results[[name]], "breseq_keys") <- comparison_keys
    }
  }

  results$merged_results <- merged_results

  # Step 4: Save to Excel
  message("\n[5/7] Saving results to Excel...")

  excel_file <- file.path(output_dir, "merged_mutation_results.xlsx")
  save_merged_results_to_excel(merged_results, excel_file,
    highlight_new = highlight_new, clean_headers = clean_headers,
    sample_labels = sample_labels, sample_order = sample_names)
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
          row_seqids <- .breseq_wide_seq(mut_df, .breseq_stage_columns(mut_df))
          for (rep_name in names(genome_lengths)) {
            tryCatch({
              rep_df <- mut_df[row_seqids == .breseq_seq(rep_name), , drop = FALSE]
              if (nrow(rep_df) > 0) {
                fname <- paste0("circos_", gsub("[^a-zA-Z0-9]", "_", rep_name), ".pdf")
                plot_circos_genome(
                  mutation_data = rep_df,
                  reference_file = reference_file,
                  samples = sample_names, seqid = rep_name,
                  genome_length = unname(genome_lengths[[rep_name]]),
                  hide_intergenic_labels = hide_intergenic_labels,
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
              samples = sample_names,
              seqid = if (is.null(names(genome_lengths))) NULL else names(genome_lengths)[1],
              genome_length = unname(genome_lengths[1]),
              hide_intergenic_labels = hide_intergenic_labels,
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
      keys <- attr(merged_tbl, "breseq_keys")
      key_aware <- !is.null(keys)
      if (key_aware) {
        positions <- stats::setNames(merged_tbl[[pos_col]], keys)
        merged_tbl[[pos_col]] <- keys
      }
      for (i in seq_along(tbl_list)[-1]) {
        # Get columns that already exist in merged_tbl (excluding position)
        existing_cols <- setdiff(names(merged_tbl), pos_col)
        new_cols <- setdiff(names(tbl_list[[i]]), pos_col)

        # Only add columns that don't already exist
        cols_to_add <- setdiff(new_cols, existing_cols)
        if (length(cols_to_add) > 0) {
          add_df <- tbl_list[[i]][, c(pos_col, cols_to_add), drop = FALSE]
          if (key_aware) {
            next_keys <- attr(tbl_list[[i]], "breseq_keys")
            if (is.null(next_keys) || length(next_keys) != nrow(add_df))
              stop("Cannot mix keyed and legacy comparisons; rebuild the comparisons together.")
            positions <- c(positions, stats::setNames(add_df[[pos_col]], next_keys))
            add_df[[pos_col]] <- next_keys
          }
          merged_tbl <- merge(merged_tbl, add_df, by = pos_col, all = TRUE)
        }
      }

      # Remove status columns and duplicates
      merged_tbl <- merged_tbl[, !grepl("^status", names(merged_tbl))]
      if (key_aware) {
        keys <- merged_tbl[[pos_col]]
        merged_tbl[[pos_col]] <- unname(positions[match(keys, names(positions))])
        attr(merged_tbl, "breseq_keys") <- keys
      }
      attr(merged_tbl, "breseq_sample_order") <- attr(.breseq_stage_columns(merged_tbl), "sample_order")

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
#' @param highlight_new Shade every reported call by stage; bold new or changed calls.
#' @param clean_headers Use readable export headers with hidden canonical-column metadata.
#' @param sample_labels Named sample-ID to display-label map, used only with clean_headers=TRUE.
#' @param sample_order Optional exact chronological sample order.
#' @keywords internal
save_merged_results_to_excel <- function(merged_results, output_file,
                                         highlight_new = FALSE, clean_headers = FALSE,
                                         sample_labels = NULL, sample_order = NULL) {

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
    mutation_name <- names(merged_results)[vapply(names(merged_results), function(n) .breseq_table_kind(n) == "mutation", logical(1))]
    mutation_data <- if (length(mutation_name)) merged_results[[mutation_name[1]]] else NULL
    .breseq_write_presentation(wb, sheet_name, merged_results[[name]], name,
      highlight_new, clean_headers, sample_labels, sample_order, mutation_data)

    # Add conditional formatting for mutation types if present
    if (!highlight_new && !clean_headers && "mutation_type" %in% names(merged_results[[name]])) {
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

# Bundled: save_comparison_results.R
#' Save Comparison Results to Excel
#'
#' This function saves a list of data frames (comparison results) into an Excel workbook,
#' with each data frame written to a separate sheet. Sheet names are simplified for clarity.
#'
#' @param comparison_results A named list of data frames, typically generated by `compare_named_tables`.
#' @param file_name A character string specifying the output Excel file name.
#' @param highlight_new Shade every reported call by stage and bold new or
#'   changed calls. Default FALSE; this does not imply de novo mutations.
#' @param clean_headers Use readable export headers and hidden canonical-column metadata.
#' @param sample_labels Optional named character map from sample IDs to display labels;
#'   used only when clean_headers is TRUE. Internal IDs remain unchanged.
#' @param sample_order Optional exact sample IDs in chronological order.
#' @return None. The function saves the results directly to the specified Excel file.
#' @import openxlsx
#' @export
#' @examples
#' # Example usage:
#' save_comparison_results_to_excel(comparison_results, "output.xlsx")
save_comparison_results_to_excel <- function(comparison_results, file_name,
                                             highlight_new = FALSE, clean_headers = FALSE,
                                             sample_labels = NULL, sample_order = NULL) {

  # Create a new workbook
  wb <- openxlsx::createWorkbook()

  # Define a mapping for simplified sheet names
  sheet_names <- c("Mutation predictions" = "Mutations",
                   "Unassigned missing coverage evidence" = "Missing",
                   "Unassigned new junction evidence" = "New junction")

  # Loop through each data frame in `comparison_results` and add it as a new sheet
  for (name in names(comparison_results)) {
    # Use the simplified name from the mapping or default to the original name
    simplified_name <- sheet_names[name]
    if (is.na(simplified_name)) simplified_name <- name  # Fallback to original name if not mapped

    # Ensure sheet name is within the Excel limit (31 characters)
    simplified_name <- substr(simplified_name, 1, 31)

    # Add a new sheet and write the data frame to it
    openxlsx::addWorksheet(wb, simplified_name)
    mutation_names <- names(comparison_results)[vapply(names(comparison_results), function(n) .breseq_table_kind(n) == "mutation", logical(1))]
    mutation_data <- if (length(mutation_names)) comparison_results[[mutation_names[1]]] else NULL
    .breseq_write_presentation(wb, simplified_name, comparison_results[[name]], name,
      highlight_new, clean_headers, sample_labels, sample_order, mutation_data)
  }

  # Save the workbook to the specified file
  openxlsx::saveWorkbook(wb, file_name, overwrite = TRUE)

  # Confirm save location
  cat("Excel file saved to:", file_name, "\n")
}

# Bundled: save_results_all.R
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

# Bundled: sequence_edits.R
# Reconstruct one sample/replicon from finalized breseq mutation calls.
# Coordinate rules: https://breseq.barricklab.org/latest/genomediff-file-format/
# HTML repeat notation retains the underlying INS/DEL position, rather than
# changing it to the start of the repeat tract (breseq src/breseq/output.cpp;
# MutationPredictor::normalize_and_annotate_tandem_repeat_mutations).
# No unassigned MC evidence or coverage-derived plotting intervals belong here.
.breseq_apply_sequence_edits <- function(sequence, calls, strict = TRUE) {
  if (!is.character(sequence) || length(sequence) != 1L || is.na(sequence) ||
      !grepl("^[ACGTRYSWKMBDHVN]*$", toupper(sequence))) {
    stop("sequence must be one DNA string containing only IUPAC nucleotide symbols.")
  }
  if (!is.data.frame(calls)) stop("calls must be a data frame.")
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("strict must be a single TRUE/FALSE value.")
  }
  sequence <- toupper(sequence)
  reference_length <- nchar(sequence)
  empty_edits <- data.frame(position = numeric(), mutation = character(), type = character(),
    ref_start = numeric(), ref_end = numeric(), ref_length = numeric(), ref_sequence = character(),
    alt_length = numeric(), alt_sequence = character(), new_start = numeric(), new_end = numeric(),
    source_rows = character(), stringsAsFactors = FALSE)
  report <- data.frame(row = seq_len(nrow(calls)), position = rep(NA_character_, nrow(calls)),
    mutation = rep(NA_character_, nrow(calls)), status = rep("pending", nrow(calls)),
    message = rep("", nrow(calls)), edit_index = rep(NA_integer_, nrow(calls)),
    stringsAsFactors = FALSE)
  if (!nrow(calls)) return(list(sequence = sequence, edits = empty_edits, report = report,
                              complete = TRUE, reference_length = reference_length))
  if (!all(c("position", "mutation") %in% names(calls))) {
    stop("calls must contain position and mutation columns.")
  }
  seq_columns <- which(gsub("[[:space:]\u00a0_]", "", names(calls), perl = TRUE) == "seqid")
  if (length(seq_columns)) {
    ids <- unlist(lapply(calls[seq_columns], as.character), use.names = FALSE)
    ids <- unique(sub("\\.[0-9]+$", "", ids[!is.na(ids) & nzchar(ids)]))
    if (length(ids) > 1L) stop("calls contain multiple replicons; apply edits to one replicon at a time.")
  }
  report$position <- as.character(calls$position)
  report$mutation <- as.character(calls$mutation)
  # Nonbreaking/zero-width spaces come from breseq's HTML display.
  normalized <- gsub("[[:space:]\u00a0\u200b\ufeff]", "", report$mutation, perl = TRUE)
  normalized <- gsub("\u2192", "->", normalized, fixed = TRUE)
  normalized <- gsub("\u2212", "-", normalized, fixed = TRUE)
  normalized <- gsub("\u2206", "\u0394", normalized, fixed = TRUE)
  positions <- gsub(",", "", trimws(report$position), fixed = TRUE)
  matches <- function(pattern, text) regmatches(text, regexec(pattern, text, perl = TRUE))[[1]]
  fail <- function(status, message) {
    stop(structure(list(message = message, call = NULL, status = status),
                   class = c("breseq_edit_parse_error", "error", "condition")))
  }
  parsed <- vector("list", nrow(calls))
  for (i in seq_len(nrow(calls))) {
    parsed[i] <- list(tryCatch({
      if (is.na(positions[i]) || !grepl("^[0-9]+$", positions[i])) {
        fail("invalid_position", "Expected an integer reference coordinate, without ranges or junction symbols.")
      }
      position <- as.numeric(positions[i])
      if (!is.finite(position) || position < 1 || position > reference_length) {
        fail("invalid_position", "Position is outside the reference sequence.")
      }
      mutation <- normalized[i]
      if (is.na(mutation) || !nzchar(mutation)) fail("unsupported", "Mutation is missing.")
      mutation <- toupper(mutation)
      ref_start <- position
      ref_length <- 0
      alt <- ""
      expected <- NULL
      type <- NULL
      repeat_unit <- NULL
      old_copies <- NULL
      m <- matches("^([ACGT])->([ACGT])$", mutation)
      if (length(m)) {
        type <- "SNP"; ref_length <- 1; expected <- m[2]; alt <- m[3]
      } else if (grepl("^\\+[ACGT]+$", mutation)) {
        type <- "INS"; ref_start <- position + 1; alt <- substring(mutation, 2)
      } else if (length(m <- matches("^[-\u0394]([0-9][0-9,]*)(?:BP)?$", mutation))) {
        type <- "DEL"; ref_length <- as.numeric(gsub(",", "", m[2], fixed = TRUE))
      } else if (length(m <- matches("^-([ACGT]+)$", mutation))) {
        type <- "DEL"; expected <- m[2]; ref_length <- nchar(expected)
      } else if (length(m <- matches("^([0-9][0-9,]*)BP->([ACGT]+)$", mutation))) {
        type <- "REP"; ref_length <- as.numeric(gsub(",", "", m[2], fixed = TRUE)); alt <- m[3]
      } else if (length(m <- matches("^\\(([ACGT]+)\\)([0-9]+)->([0-9]+)$", mutation))) {
        repeat_unit <- m[2]; old_copies <- as.numeric(m[3]); new_copies <- as.numeric(m[4])
        difference <- new_copies - old_copies
        if (!is.finite(old_copies) || old_copies < 1 || !is.finite(new_copies) || difference == 0) {
          fail("invalid_repeat", "Repeat copy numbers must specify a change from a positive reference count.")
        }
        if (difference > 0) {
          type <- "INS"; ref_start <- position + 1
          alt <- strrep(repeat_unit, difference)
          repeat_end <- position
        } else {
          type <- "DEL"; ref_length <- -difference * nchar(repeat_unit)
          repeat_end <- position + ref_length - 1
          expected <- strrep(repeat_unit, -difference)
        }
        repeat_start <- repeat_end - old_copies * nchar(repeat_unit) + 1
        if (repeat_start < 1 || repeat_end > reference_length ||
            substr(sequence, repeat_start, repeat_end) != strrep(repeat_unit, old_copies)) {
          fail("reference_mismatch", "The reported repeat sequence/count does not match the reference at the breseq INS/DEL coordinate.")
        }
        # breseq counts complete units backwards at this phase, stopping before
        # reference position 1. Preserve that convention for origin-edge calls.
        preceding_start <- repeat_start - nchar(repeat_unit)
        if (preceding_start > 1 && substr(sequence, preceding_start, repeat_start - 1) == repeat_unit) {
          fail("reference_mismatch", "The reference contains more preceding repeat copies than reported.")
        }
      } else {
        fail("unsupported", "Unsupported or sequence-abbreviated mutation; supply an exact finalized call/GenomeDiff record.")
      }
      if (!is.finite(ref_length) || ref_length < 0 ||
          (type %in% c("DEL", "REP") && ref_length == 0)) {
        fail("invalid_length", "Deletion/replacement length must be a positive integer.")
      }
      ref_end <- ref_start + ref_length - 1
      if (ref_end > reference_length) fail("invalid_position", "Edit crosses the reference end; an origin-spanning edit needs explicit resolution.")
      reference <- if (ref_length) substr(sequence, ref_start, ref_end) else ""
      if (!is.null(expected) && !identical(reference, expected)) {
        fail("reference_mismatch", paste0("Reference allele mismatch: expected ", expected, ", found ", reference, "."))
      }
      if (ref_length && identical(reference, alt)) fail("no_change", "The call does not change the reference sequence.")
      report$status[i] <- "validated"
      data.frame(position = position, mutation = report$mutation[i], type = type,
        ref_start = ref_start, ref_end = ref_end, ref_length = ref_length, ref_sequence = reference,
        alt_length = nchar(alt), alt_sequence = alt, new_start = NA_real_, new_end = NA_real_,
        source_rows = as.character(i), stringsAsFactors = FALSE)
    }, breseq_edit_parse_error = function(e) {
      report$status[i] <<- e$status
      report$message[i] <<- conditionMessage(e)
      NULL
    }))
  }
  edits <- if (any(lengths(parsed))) do.call(rbind, parsed) else empty_edits
  # Equivalent exact edits are safe to collapse; every source row remains in the
  # report, including duplicate evidence rows and equivalent display spellings.
  if (nrow(edits)) {
    key <- paste(edits$ref_start, edits$ref_end, edits$alt_sequence, sep = "\r")
    unique_key <- unique(key)
    kept <- match(unique_key, key)
    sources <- lapply(unique_key, function(k) as.integer(edits$source_rows[key == k]))
    edits <- edits[kept, , drop = FALSE]
    edits$source_rows <- vapply(sources, paste, collapse = ",", character(1))
    edits <- edits[order(edits$ref_start, edits$ref_end), , drop = FALSE]
    rownames(edits) <- NULL
    sources <- lapply(strsplit(edits$source_rows, ",", fixed = TRUE), as.integer)
    for (j in seq_len(nrow(edits))) {
      report$edit_index[sources[[j]]] <- j
      if (length(sources[[j]]) > 1L) {
        duplicates <- sources[[j]][-1]
        report$status[duplicates] <- "duplicate"
        report$message[duplicates] <- paste0("Same exact edit as source row ", sources[[j]][1], ".")
      }
    }
    conflict <- logical(nrow(edits))
    if (nrow(edits) > 1L) for (a in seq_len(nrow(edits) - 1L)) for (b in seq.int(a + 1L, nrow(edits))) {
      ia <- edits$ref_length[a] == 0; ib <- edits$ref_length[b] == 0
      overlap <- if (ia && ib) {
        edits$ref_end[a] == edits$ref_end[b]
      } else if (!ia && !ib) {
        max(edits$ref_start[a], edits$ref_start[b]) <= min(edits$ref_end[a], edits$ref_end[b])
      } else {
        ins <- if (ia) a else b; other <- if (ia) b else a
        edits$ref_end[ins] >= edits$ref_start[other] && edits$ref_end[ins] < edits$ref_end[other]
      }
      if (overlap) conflict[c(a, b)] <- TRUE
    }
    for (j in which(conflict)) {
      report$status[sources[[j]]] <- "overlap"
      report$message[sources[[j]]] <- "Overlapping edits or differently ordered insertions require explicit resolution; no sequence was generated."
    }
  }
  invalid <- !report$status %in% c("validated", "duplicate")
  if (any(invalid)) {
    result <- list(sequence = NULL, edits = edits, report = report, complete = FALSE,
                   reference_length = reference_length)
    if (!strict) return(result)
    details <- paste0("row ", report$row[invalid], " (", report$position[invalid], ", ",
                      report$mutation[invalid], "): ", report$message[invalid])
    message <- paste("Cannot reconstruct sequence:", paste(utils::head(details, 8), collapse = "; "))
    if (sum(invalid) > 8L) message <- paste0(message, "; see $report for all ", sum(invalid), " unresolved rows.")
    stop(structure(c(list(message = message, call = NULL), result),
                   class = c("breseq_sequence_edit_error", "error", "condition")))
  }
  pieces <- character(2 * nrow(edits) + 1L)
  cursor <- 1
  shift <- 0
  for (j in seq_len(nrow(edits))) {
    pieces[2 * j - 1L] <- if (edits$ref_start[j] > cursor) substr(sequence, cursor, edits$ref_start[j] - 1) else ""
    pieces[2 * j] <- edits$alt_sequence[j]
    edits$new_start[j] <- edits$ref_start[j] + shift
    edits$new_end[j] <- edits$new_start[j] + edits$alt_length[j] - 1
    cursor <- edits$ref_end[j] + 1
    shift <- shift + edits$alt_length[j] - edits$ref_length[j]
  }
  pieces[length(pieces)] <- if (cursor <= reference_length) substr(sequence, cursor, reference_length) else ""
  mutated <- paste0(pieces, collapse = "")
  stopifnot(nchar(mutated) == reference_length + sum(edits$alt_length - edits$ref_length))
  report$status[report$status == "validated"] <- "applied"
  list(sequence = mutated, edits = edits, report = report, complete = TRUE,
       reference_length = reference_length)
}

# Bundled: sequence_genbank.R
# GenBank reconstruction uses edited bases, independently of display-only marks.
.breseq_gb_split_location <- function(text) {
  chars <- strsplit(text, "", fixed = TRUE)[[1]]
  depth <- 0L; starts <- 1L; result <- character()
  for (i in seq_along(chars)) {
    if (chars[i] == "(") depth <- depth + 1L
    if (chars[i] == ")") depth <- depth - 1L
    if (chars[i] == "," && depth == 0L) {
      result <- c(result, substr(text, starts, i - 1L)); starts <- i + 1L
    }
  }
  c(result, substr(text, starts, nchar(text)))
}

.breseq_gb_location <- function(text) {
  text <- gsub("[[:space:]]", "", text)
  if (grepl("^(complement|join|order)\\(", text)) {
    operator <- sub("\\(.*", "", text)
    if (!endsWith(text, ")")) stop("Unbalanced feature location")
    inner <- substr(text, nchar(operator) + 2L, nchar(text) - 1L)
    parts <- .breseq_gb_split_location(inner)
    if (operator == "complement" && length(parts) != 1L) stop("Invalid complement")
    return(list(type = operator, children = lapply(parts, .breseq_gb_location)))
  }
  if (grepl("^[<>]?[0-9]+(\\.\\.[<>]?[0-9]+)?$", text)) {
    parts <- strsplit(text, "..", fixed = TRUE)[[1]]
    if (length(parts) == 1L) parts <- rep(parts, 2L)
    return(list(type = "range", start = as.numeric(gsub("[<>]", "", parts[1])),
      end = as.numeric(gsub("[<>]", "", parts[2])),
      left = sub("[0-9].*", "", parts[1]), right = sub("[0-9].*", "", parts[2]),
      point = !grepl("..", text, fixed = TRUE)))
  }
  if (grepl("^[0-9]+\\^[0-9]+$", text)) {
    parts <- as.numeric(strsplit(text, "^", fixed = TRUE)[[1]])
    return(list(type = "between", start = parts[1], end = parts[2]))
  }
  stop("Unsupported feature location: ", text)
}

.breseq_gb_location_text <- function(node) {
  if (node$type == "range") {
    left <- paste0(node$left, node$start)
    if (isTRUE(node$point) && node$start == node$end) return(left)
    return(paste0(left, "..", node$right, node$end))
  }
  if (node$type == "between") return(paste0(node$start, "^", node$end))
  paste0(node$type, "(", paste(vapply(node$children, .breseq_gb_location_text, ""), collapse = ","), ")")
}

.breseq_gb_map_base <- function(position, edits) {
  position + sum((edits$alt_length - edits$ref_length)[edits$ref_end < position])
}

.breseq_gb_lift_location <- function(node, edits, reference_length, new_length) {
  if (!node$type %in% c("range", "between")) {
    lifted <- lapply(node$children, .breseq_gb_lift_location, edits = edits,
      reference_length = reference_length, new_length = new_length)
    children <- lapply(lifted, `[[`, "node")
    children <- Filter(Negate(is.null), children)
    changed <- any(vapply(lifted, `[[`, logical(1), "changed"))
    partial <- length(children) != length(node$children) || any(vapply(lifted, `[[`, logical(1), "partial"))
    if (!length(children)) return(list(node = NULL, changed = TRUE, partial = TRUE))
    node$children <- children
    if (node$type %in% c("join", "order") && length(children) == 1L) node <- children[[1]]
    return(list(node = node, changed = changed, partial = partial))
  }
  a <- node$start; b <- node$end
  if (node$type == "between") {
    if (!(b == a + 1 || (a == reference_length && b == 1))) stop("Invalid between-base location")
    inside <- edits$ref_length > 0 & edits$ref_start <= a & edits$ref_end > a
    if (any(inside)) return(list(node = NULL, changed = TRUE, partial = TRUE))
    boundary <- a + sum((edits$alt_length - edits$ref_length)[edits$ref_end <= a])
    if (a == reference_length && b == 1) {
      node$start <- new_length; node$end <- 1
    } else {
      if (boundary < 1 || boundary >= new_length) return(list(node = NULL, changed = TRUE, partial = TRUE))
      node$start <- boundary; node$end <- boundary + 1
    }
    return(list(node = node, changed = any(edits$ref_length == 0 & edits$ref_end == a), partial = FALSE))
  }
  if (a < 1 || b < a || b > reference_length) stop("Out-of-reference feature location")
  relevant <- which((edits$ref_length > 0 & edits$ref_start <= b & edits$ref_end >= a) |
    (edits$ref_length == 0 & edits$ref_end >= a & edits$ref_end < b))
  if (!length(relevant)) {
    node$start <- .breseq_gb_map_base(a, edits); node$end <- .breseq_gb_map_base(b, edits)
    return(list(node = node, changed = FALSE, partial = FALSE))
  }
  segments <- list(); cursor <- a; partial_left <- FALSE; partial_right <- FALSE
  add <- function(left, right) { if (left <= right) segments[[length(segments) + 1L]] <<- c(left, right) }
  for (i in relevant) {
    e <- edits[i, ]
    if (e$ref_length == 0) {
      add(e$new_start, e$new_end)
      next
    }
    before <- min(b, e$ref_start - 1)
    if (cursor <= before) add(.breseq_gb_map_base(cursor, edits), .breseq_gb_map_base(before, edits))
    if (e$alt_length > 0) {
      if (e$ref_length == e$alt_length) {
        add(e$new_start + max(a, e$ref_start) - e$ref_start,
            e$new_start + min(b, e$ref_end) - e$ref_start)
      } else if (a <= e$ref_start && b >= e$ref_end) add(e$new_start, e$new_end)
    }
    if (e$ref_length != e$alt_length && e$ref_start <= a && e$ref_end >= a) partial_left <- TRUE
    if (e$ref_length != e$alt_length && e$ref_start <= b && e$ref_end >= b) partial_right <- TRUE
    cursor <- max(cursor, e$ref_end + 1)
  }
  if (cursor <= b) add(.breseq_gb_map_base(cursor, edits), .breseq_gb_map_base(b, edits))
  if (!length(segments)) return(list(node = NULL, changed = TRUE, partial = TRUE))
  node$start <- min(vapply(segments, `[`, 0, 1)); node$end <- max(vapply(segments, `[`, 0, 2))
  if (partial_left) node$left <- "<"
  if (partial_right) node$right <- ">"
  if (node$start != node$end) node$point <- FALSE
  list(node = node, changed = TRUE, partial = partial_left || partial_right)
}

.breseq_gb_reverse_complement <- function(sequence) {
  paste(rev(strsplit(chartr("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN", sequence), "", fixed = TRUE)[[1]]), collapse = "")
}

.breseq_gb_feature_sequence <- function(node, sequence) {
  if (node$type == "range") return(substr(sequence, node$start, node$end))
  if (node$type == "between") return("")
  if (node$type == "order") stop("Ordered segments cannot be translated as a joined CDS")
  joined <- paste(vapply(node$children, .breseq_gb_feature_sequence, "", sequence = sequence), collapse = "")
  if (node$type == "complement") .breseq_gb_reverse_complement(joined) else joined
}

.breseq_gb_translation <- function(node, sequence, code = 1L, codon_start = 1L) {
  if (!code %in% c(1L, 11L)) stop("unsupported translation table")
  if (grepl("[<>]", .breseq_gb_location_text(node)) || codon_start != 1L) stop("partial CDS")
  dna <- .breseq_gb_feature_sequence(node, sequence)
  if (!nchar(dna) || nchar(dna) %% 3L != 0L) stop("CDS length is not a multiple of three")
  if (grepl("[^ACGT]", dna)) stop("ambiguous coding bases")
  starts <- seq.int(1L, nchar(dna), by = 3L)
  codons <- substring(dna, starts, starts + 2L)
  bases <- c("T", "C", "A", "G")
  # Explicit NCBI standard-code order: first base slowest, third base fastest.
  keys <- unlist(lapply(bases, function(a) unlist(lapply(bases, function(b) paste0(a, b, bases)))))
  amino <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "", fixed = TRUE)[[1]]
  protein <- unname(stats::setNames(amino, keys)[codons])
  allowed_starts <- if (code == 11L) c("TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG") else c("TTG", "CTG", "ATG")
  if (!codons[1] %in% allowed_starts) stop("start codon is not valid for this translation table")
  if (tail(protein, 1) != "*" || any(head(protein, -1L) == "*")) stop("missing terminal stop or internal stop codon")
  protein[1] <- "M"
  paste(head(protein, -1L), collapse = "")
}

.breseq_gb_qualifier <- function(key, value) {
  value <- gsub('"', "'", as.character(value), fixed = TRUE)
  chunks <- character(); remaining <- value
  repeat {
    width <- if (!length(chunks)) max(1L, 76L - 21L - nchar(key)) else 58L
    cut <- min(width, nchar(remaining))
    if (key != "translation" && nchar(remaining) > cut) {
      spaces <- gregexpr(" ", substr(remaining, 1L, cut), fixed = TRUE)[[1]]
      if (any(spaces > 0)) cut <- max(spaces)
    }
    chunks <- c(chunks, substr(remaining, 1L, cut))
    remaining <- substring(remaining, cut + 1L)
    if (key != "translation") remaining <- sub("^ +", "", remaining)
    if (!nzchar(remaining)) break
  }
  lines <- c(paste0(strrep(" ", 21), "/", key, '="', chunks[1]),
             if (length(chunks) > 1L) paste0(strrep(" ", 21), chunks[-1]) else character())
  lines[length(lines)] <- paste0(lines[length(lines)], '"')
  lines
}

.breseq_gb_feature_lines <- function(key, location) {
  starts <- seq.int(1L, nchar(location), by = 59L)
  chunks <- substring(location, starts, starts + 58L)
  c(paste0("     ", sprintf("%-16s", key), chunks[1]),
    if (length(chunks) > 1L) paste0(strrep(" ", 21), chunks[-1]) else character())
}

.breseq_gb_features <- function(lines) {
  # Some complete reference records have CONTIG/BASE COUNT footer fields even
  # when ORIGIN bases are present. They are not feature qualifier continuations.
  footer <- grep("^[^ ]", lines)
  if (length(footer)) lines <- lines[seq_len(footer[1] - 1L)]
  starts <- grep("^ {5}[^ ]", lines)
  if (!length(starts)) return(list())
  ends <- c(starts[-1] - 1L, length(lines))
  lapply(seq_along(starts), function(i) {
    block <- lines[starts[i]:ends[i]]
    q <- grep("^ {21}/", block)
    location_end <- if (length(q)) q[1] - 1L else length(block)
    qualifiers <- if (length(q)) lapply(seq_along(q), function(j) {
      last <- if (j < length(q)) q[j + 1L] - 1L else length(block)
      block[q[j]:last]
    }) else list()
    list(key = trimws(substr(block[1], 6L, 20L)),
      location = paste(trimws(substring(block[seq_len(location_end)], 22L)), collapse = ""),
      qualifiers = qualifiers)
  })
}

.breseq_gb_qualifier_key <- function(block) sub("[=[:space:]].*$", "", sub("^ */", "", block[1]))
.breseq_gb_qualifier_value <- function(block) {
  value <- paste(trimws(block), collapse = "")
  sub('"$', "", sub('^[^=]*="?', "", value))
}

.breseq_gb_read_records <- function(reference_file) {
  lines <- readLines(reference_file, warn = FALSE)
  if (any(grepl("^LOCUS", lines))) {
    starts <- grep("^LOCUS", lines); ends <- grep("^//", lines)
    if (length(starts) != length(ends)) stop("Incomplete GenBank reference records")
    return(lapply(seq_along(starts), function(i) {
      record <- lines[starts[i]:ends[i]]
      origin <- grep("^ORIGIN", record)[1]; features <- grep("^FEATURES", record)[1]
      if (is.na(origin)) stop("Reference record has no ORIGIN sequence")
      id <- strsplit(trimws(record[1]), "[[:space:]]+")[[1]][2]
      version <- grep("^VERSION[[:space:]]", record, value = TRUE)
      accession <- grep("^ACCESSION[[:space:]]", record, value = TRUE)
      if (length(version)) id <- strsplit(trimws(version[1]), "[[:space:]]+")[[1]][2] else
        if (length(accession)) id <- strsplit(trimws(accession[1]), "[[:space:]]+")[[1]][2]
      sequence <- toupper(gsub("[^A-Za-z]", "", paste(record[(origin + 1L):(length(record) - 1L)], collapse = "")))
      if (!nzchar(sequence) || grepl("[^ACGTRYMKSWBDHVN]", sequence)) stop("Invalid reference DNA bases")
      list(seq_id = .breseq_seq(id), sequence = sequence,
        topology = if (grepl("circular", record[1], fixed = TRUE)) "circular" else "linear",
        header = record[seq_len(if (is.na(features)) origin - 1L else features - 1L)],
        features = if (is.na(features)) list() else .breseq_gb_features(record[(features + 1L):(origin - 1L)]))
    }))
  }
  starts <- grep("^>", lines)
  if (!length(starts)) stop("Reference must be a GenBank or FASTA file")
  ends <- c(starts[-1] - 1L, length(lines))
  lapply(seq_along(starts), function(i) {
    id <- strsplit(sub("^>", "", lines[starts[i]]), "[[:space:]]+")[[1]][1]
    sequence <- toupper(gsub("[[:space:]]", "", paste(lines[(starts[i] + 1L):ends[i]], collapse = "")))
    if (!nzchar(sequence) || grepl("[^ACGTRYMKSWBDHVN]", sequence)) stop("Invalid reference DNA bases")
    topology <- if (grepl("linear", lines[starts[i]], ignore.case = TRUE)) "linear" else "circular"
    list(seq_id = .breseq_seq(id), sequence = sequence, topology = topology,
         header = c(paste0("DEFINITION  ", sub("^>", "", lines[starts[i]]))), features = list())
  })
}

.breseq_gb_lift_features <- function(record, edited_sequence, edits) {
  output <- character(); audit <- list(); length_new <- nchar(edited_sequence)
  features <- record$features
  if (!any(vapply(features, function(x) x$key == "source", logical(1))))
    features <- c(list(list(key = "source", location = paste0("1..", nchar(record$sequence)), qualifiers = list())), features)
  for (i in seq_along(features)) {
    f <- features[[i]]; q <- f$qualifiers
    keys <- vapply(q, .breseq_gb_qualifier_key, "")
    locus <- which(keys %in% c("locus_tag", "gene"))
    label <- if (length(locus)) .breseq_gb_qualifier_value(q[[locus[1]]]) else ""
    old_location <- f$location
    result <- tryCatch({
      node <- .breseq_gb_location(old_location)
      if (f$key == "source") {
        node <- .breseq_gb_location(paste0("1..", length_new))
        list(node = node, changed = FALSE, partial = FALSE)
      } else .breseq_gb_lift_location(node, edits, nchar(record$sequence), length_new)
    }, error = function(e) list(node = NULL, changed = FALSE, partial = FALSE, error = conditionMessage(e)))
    status <- "preserved"; note <- ""; location <- ""
    if (is.null(result$node)) {
      status <- if (!is.null(result$error)) "dropped_unsupported_location" else "dropped_fully_deleted"
      note <- if (!is.null(result$error)) result$error else "All supported location segments were deleted"
    } else {
      location <- .breseq_gb_location_text(result$node)
      if (location != old_location) status <- "lifted"
      if (result$changed) status <- if (result$partial) "partial_sequence_change" else "sequence_changed"
      # Positional qualifiers must not retain old coordinates after liftover.
      positional_failure <- FALSE
      for (j in which(keys %in% c("transl_except", "anticodon", "rpt_unit_range"))) {
        value <- .breseq_gb_qualifier_value(q[[j]])
        old <- if (keys[j] == "rpt_unit_range") value else sub(",aa:.*$", "", sub("^\\(pos:", "", value))
        lifted <- tryCatch(.breseq_gb_lift_location(.breseq_gb_location(old), edits,
          nchar(record$sequence), length_new), error = function(e) NULL)
        if (is.null(lifted) || is.null(lifted$node) || lifted$changed) {
          q[[j]] <- character(); positional_failure <- TRUE
        } else {
          replacement <- .breseq_gb_location_text(lifted$node)
          q[[j]] <- .breseq_gb_qualifier(keys[j], sub(old, replacement, value, fixed = TRUE))
        }
      }
      if (f$key == "CDS" && (result$changed || positional_failure)) {
        code_at <- which(keys == "transl_table"); start_at <- which(keys == "codon_start")
        code <- if (length(code_at)) as.integer(.breseq_gb_qualifier_value(f$qualifiers[[code_at[1]]])) else 1L
        offset <- if (length(start_at)) as.integer(.breseq_gb_qualifier_value(f$qualifiers[[start_at[1]]])) else 1L
        translated <- tryCatch({
          if (result$partial || positional_failure || any(keys %in% c("pseudo", "pseudogene", "transl_except"))) stop("partial or exceptional CDS requires reannotation")
          .breseq_gb_translation(result$node, edited_sequence, code, offset)
        }, error = function(e) e)
        q <- q[!keys %in% c("translation", "protein_id", "inference")]
        if (inherits(translated, "error")) {
          if (!any(keys %in% c("pseudo", "pseudogene"))) q <- c(q, list(paste0(strrep(" ", 21), "/pseudo")))
          note <- paste0("Altered CDS; translation not asserted: ", conditionMessage(translated),
                         ". Reference annotation retained conservatively; functional status unverified.")
          status <- "CDS_requires_reannotation"
        } else {
          q <- c(q, list(.breseq_gb_qualifier("translation", translated)))
          note <- "CDS translation recomputed from the called-variant sequence; reference protein_id removed."
          status <- "CDS_translation_recomputed"
        }
      } else if (result$changed || positional_failure) {
        note <- "Reference feature lifted across called sequence changes; altered or truncated annotation requires review."
      }
      if (nzchar(note)) q <- c(q, list(.breseq_gb_qualifier("note", note)))
      output <- c(output, .breseq_gb_feature_lines(f$key, location), unlist(q, use.names = FALSE))
    }
    audit[[i]] <- data.frame(feature_key = f$key, locus = label, original_location = old_location,
      new_location = location, status = status, note = note, stringsAsFactors = FALSE)
  }
  list(lines = output, audit = do.call(rbind, audit))
}

.breseq_gb_header <- function(record, id, sample, length_new) {
  header <- record$header
  # Reference accession IDs must not identify a newly reconstructed sequence.
  starts <- grep("^[A-Z]", header)
  ends <- c(starts[-1] - 1L, length(header))
  keep <- rep(TRUE, length(header))
  for (i in seq_along(starts)) if (trimws(substr(header[starts[i]], 1, 12)) %in%
      c("LOCUS", "ACCESSION", "VERSION", "BASE COUNT", "CONTIG")) keep[starts[i]:ends[i]] <- FALSE
  header <- header[keep]
  date <- Sys.Date(); months <- c("JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC")
  date <- paste(format(date, "%d"), months[as.integer(format(date, "%m"))], format(date, "%Y"), sep = "-")
  locus <- sprintf("LOCUS       %-16s %11d bp    DNA     %-8s BCT %s", id, length_new, record$topology, date)
  # Insert local identifiers after DEFINITION, ahead of source/reference metadata.
  at <- grep("^[A-Z]", header)
  after_definition <- if (length(at) > 1L) at[2] - 1L else length(header)
  header <- append(header, c(paste0("ACCESSION   ", id), paste0("VERSION     ", id)), after_definition)
  header <- sub(paste0("(bases 1 to ", nchar(record$sequence), ")"),
                paste0("(bases 1 to ", length_new, ")"), header, fixed = TRUE)
  c(locus, header,
    paste0("COMMENT     Called-variant sequence for sample ", sample, "; reference ", record$seq_id, "."),
    "            Uncertain missing-coverage intervals are not applied as exact edits.")
}

.breseq_gb_origin <- function(sequence) {
  starts <- seq.int(1L, nchar(sequence), by = 60L)
  vapply(starts, function(start) {
    part <- tolower(substr(sequence, start, start + 59L))
    chunks <- substring(part, seq.int(1L, nchar(part), by = 10L), seq.int(1L, nchar(part), by = 10L) + 9L)
    paste0(sprintf("%9d", start), " ", paste(chunks, collapse = " "))
  }, "")
}

.breseq_gb_edit_features <- function(edits, length_new, sample) {
  output <- character()
  for (i in seq_len(nrow(edits))) {
    edit <- edits[i, ]
    if (edit$alt_length > 0) {
      location <- if (edit$new_start == edit$new_end) as.character(edit$new_start) else
        paste0(edit$new_start, "..", edit$new_end)
      note <- if (edit$ref_length == 0) paste0("Insertion after reference position ", edit$ref_end, ".") else
        paste0("Reference interval ", edit$ref_start, "..", edit$ref_end, ".")
    } else {
      location <- as.character(max(1, min(length_new, edit$new_start - 1)))
      note <- paste0("Deletion junction marker on an adjacent remaining base; deleted reference interval ",
        edit$ref_start, "..", edit$ref_end, ". No deleted bases remain at this marker.")
    }
    note <- paste0(note, " Sample ", sample, "; applied ", edit$mutation,
                   "; reference length ", edit$ref_length, "; alternate length ", edit$alt_length, ".")
    color <- unname(.get_genbank_mutation_color(edit$type))
    output <- c(output, .breseq_gb_feature_lines("misc_feature", location),
      .breseq_gb_qualifier("label", paste0("[", edit$type, ":", sample, "] ", edit$mutation)),
      .breseq_gb_qualifier("note", note), .breseq_gb_qualifier("ApEinfo_fwdcolor", color),
      .breseq_gb_qualifier("ApEinfo_revcolor", color))
  }
  output
}

#' Export GenBank Files Containing Actually Reconstructed Sample Sequences
#'
#' Applies only final mutation calls, lifts reference features through sequence
#' edits, and recomputes safe table 1/11 CDS translations. Display-only coverage
#' evidence is excluded. Full plasmid losses have manifest entries without files.
#' @param reference_file GenBank or FASTA reference, including multiple records.
#' @param mutation_list Named processed per-sample final mutation tables.
#' @param output_dir Directory for per-sample, per-replicon GenBank and FASTA files.
#' @param samples Sample IDs in output order; a missing WT entry is an empty baseline.
#' @param sample_labels Optional named display labels for sample IDs.
#' @param replicon_losses Optional validated whole-replicon coverage-loss records
#'   with sample, seq_id, evidence_origin, and depth_ratio columns.
#' @return List with manifest, feature_audit, edits, and call_report tables.
#' @export
export_changed_sequences_to_genbank <- function(reference_file, mutation_list, output_dir,
                                                samples = names(mutation_list), sample_labels = NULL,
                                                replicon_losses = NULL) {
  records <- .breseq_gb_read_records(reference_file)
  ids <- vapply(records, `[[`, "", "seq_id")
  if (anyDuplicated(ids)) stop("Reference replicon IDs must be unique")
  if (!is.list(mutation_list) || !is.character(samples) || !length(samples) || anyNA(samples) ||
      any(!nzchar(samples)) || anyDuplicated(samples) || any(!samples %in% c(names(mutation_list), "WT")))
    stop("samples must identify unique mutation_list stages (or an empty WT baseline)")
  if (!is.null(replicon_losses)) {
    required <- c("sample", "seq_id", "evidence_origin", "depth_ratio")
    if (!all(required %in% names(replicon_losses))) stop("replicon_losses is missing required evidence columns")
    ratio <- suppressWarnings(as.numeric(replicon_losses$depth_ratio))
    if (any(!replicon_losses$sample %in% samples) || any(!.breseq_seq(replicon_losses$seq_id) %in% ids) ||
        anyNA(replicon_losses$evidence_origin) || any(replicon_losses$evidence_origin != "COVERAGE_RATIO") ||
        any(!is.finite(ratio) | ratio < 0 | ratio >= 0.1) ||
        anyDuplicated(paste(replicon_losses$sample, .breseq_seq(replicon_losses$seq_id))))
      stop("replicon_losses must contain unique, valid whole-replicon COVERAGE_RATIO evidence below 0.1")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  manifests <- list(); audits <- list(); all_edits <- list(); reports <- list()
  for (s in seq_along(samples)) for (r in seq_along(records)) {
    sample <- samples[s]; record <- records[[r]]; data <- mutation_list[[sample]]
    if (is.null(data)) data <- data.frame()
    calls <- data
    if (nrow(calls)) {
      names(calls) <- .breseq_text(names(calls))
      seq_col <- intersect(c("seq id", "seq_id", "seqid"), names(calls))
      if (!length(seq_col)) stop("Final mutation rows need a replicon ID column")
      if (any(!.breseq_seq(calls[[seq_col[1]]]) %in% ids)) stop("Mutation rows contain an unknown replicon ID")
      keep <- .breseq_seq(calls[[seq_col[1]]]) == record$seq_id
      if ("annotation" %in% names(calls)) keep <- keep &
        !(tolower(.breseq_text(calls$annotation)) %in% "coverage-derived del display")
      if ("evidence" %in% names(calls)) keep <- keep & !(calls$evidence %in% "COVERAGE_RATIO")
      calls <- calls[keep, , drop = FALSE]
    }
    loss <- !is.null(replicon_losses) && any(replicon_losses$sample == sample &
      .breseq_seq(replicon_losses$seq_id) == record$seq_id)
    result <- if (loss) NULL else tryCatch(.breseq_apply_sequence_edits(record$sequence, calls), error = function(e) e)
    status <- if (loss) "lost_by_coverage" else if (inherits(result, "error")) "unresolved_calls" else
      if (!nchar(result$sequence)) "lost_confirmed" else if (!nrow(calls)) "reference_baseline" else "exported_called_variants"
    filename <- fasta <- checksum <- ""; new_length <- if (loss) 0 else if (inherits(result, "error")) NA_real_ else nchar(result$sequence)
    label <- if (!is.null(sample_labels) && sample %in% names(sample_labels)) unname(sample_labels[sample]) else sample
    if (!loss && !inherits(result, "error")) {
      edits <- result$edits
      if (nrow(edits)) {
        if ("source_rows" %in% names(edits) && is.list(edits$source_rows)) edits$source_rows <- vapply(edits$source_rows, paste, "", collapse = ",")
        all_edits[[length(all_edits) + 1L]] <- cbind(sample = sample, seq_id = record$seq_id, edits, stringsAsFactors = FALSE)
      }
      if (!is.null(result$report) && nrow(result$report)) reports[[length(reports) + 1L]] <-
        cbind(sample = sample, seq_id = record$seq_id, result$report, stringsAsFactors = FALSE)
      if (new_length > 0) {
        lifted <- .breseq_gb_lift_features(record, result$sequence, edits)
        audits[[length(audits) + 1L]] <- cbind(sample = sample, seq_id = record$seq_id, lifted$audit, stringsAsFactors = FALSE)
        stem <- paste0(gsub("[^A-Za-z0-9_.-]", "_", sample), "_", record$seq_id, "_called_variants")
        filename <- file.path(output_dir, paste0(stem, ".gbk")); fasta <- file.path(output_dir, paste0(stem, ".fasta"))
        local_id <- sprintf("BRS_%03d_%03d", s, r)
        writeLines(c(.breseq_gb_header(record, local_id, label, new_length), "FEATURES             Location/Qualifiers",
          lifted$lines, .breseq_gb_edit_features(edits, new_length, sample),
          "ORIGIN", .breseq_gb_origin(result$sequence), "//"), filename, useBytes = TRUE)
        starts <- seq.int(1L, new_length, by = 80L)
        writeLines(c(paste0(">", local_id, " ", sample, " ", record$seq_id, " called variants"),
          substring(result$sequence, starts, starts + 79L)), fasta, useBytes = TRUE)
        hash_file <- tempfile(); writeChar(result$sequence, hash_file, eos = NULL, useBytes = TRUE)
        checksum <- unname(tools::md5sum(hash_file)); unlink(hash_file)
      }
    } else if (inherits(result, "error") && !is.null(result$report) && nrow(result$report))
      reports[[length(reports) + 1L]] <- cbind(sample = sample, seq_id = record$seq_id, result$report, stringsAsFactors = FALSE)
    manifests[[length(manifests) + 1L]] <- data.frame(sample = sample, sample_label = label, seq_id = record$seq_id,
      reference_length = nchar(record$sequence), new_length = new_length, sequence_length = new_length,
      final_call_rows = nrow(calls), status = status, file = filename, genbank_file = filename, fasta_file = fasta,
      sequence_md5 = checksum, message = if (inherits(result, "error")) conditionMessage(result) else "", stringsAsFactors = FALSE)
  }
  bind <- function(x) if (length(x)) do.call(rbind, x) else data.frame()
  list(manifest = bind(manifests), feature_audit = bind(audits), edits = bind(all_edits), call_report = bind(reports))
}

# Bundled: snapgene_export.R
# SnapGene's documented CLI requires that the application is not already open:
# https://support.snapgene.com/hc/en-us/articles/10384393330836-Command-Line-Converting-File-Formats
# Never close an existing session or terminate a delayed conversion process.
.breseq_snapgene_is_running <- function() {
  tryCatch({
    if (.Platform$OS.type == "windows") {
      result <- suppressWarnings(system2("tasklist", c("/FO", "CSV", "/NH"),
                                        stdout = TRUE, stderr = FALSE))
      pattern <- '^"SnapGene( Viewer)?[.]exe"'
    } else {
      result <- suppressWarnings(system2("ps", c("-A", "-o", "comm="),
                                        stdout = TRUE, stderr = FALSE))
      pattern <- "(^|/)SnapGene( Viewer)?([.]exe)?$"
    }
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) return(NA)
    any(grepl(pattern, trimws(result), ignore.case = TRUE))
  }, error = function(e) NA)
}

.breseq_snapgene_run <- function(executable, args, timeout) {
  if (!requireNamespace("processx", quietly = TRUE)) {
    return(list(unavailable = TRUE, pending = FALSE, status = NA_integer_,
                message = "Install the optional processx package for SnapGene conversion."))
  }
  # cleanup=FALSE is intentional: a timeout must not kill SnapGene or a prompt.
  process <- processx::process$new(executable, args, stdout = "|", stderr = "|",
                                  cleanup = FALSE, cleanup_tree = FALSE)
  process$wait(timeout = timeout * 1000)
  if (process$is_alive()) {
    return(list(unavailable = FALSE, pending = TRUE, status = NA_integer_,
                message = paste0("SnapGene did not exit within ", timeout,
                  " seconds; its process was left running and was not terminated.")))
  }
  list(unavailable = FALSE, pending = FALSE, status = process$get_exit_status(),
       message = trimws(paste(process$read_all_output(), process$read_all_error())))
}

.breseq_snapgene_dna_header <- function(file) {
  if (!file.exists(file) || file.info(file)$size <= 19L) return(FALSE)
  connection <- file(file, "rb")
  on.exit(close(connection), add = TRUE)
  header <- readBin(connection, "raw", n = 19L)
  # Cookie packet: type 9, size 14, signature SnapGene, sequence type DNA (1).
  length(header) == 19L &&
    identical(as.integer(header[1:5]), c(9L, 0L, 0L, 0L, 14L)) &&
    identical(header[6:13], charToRaw("SnapGene")) &&
    identical(as.integer(header[14:15]), c(0L, 1L))
}

# Return an explicit conversion status; never fabricate a .dna by renaming text.
# A single edited GenBank replicon is expected.  The caller retains that GenBank
# whether native SnapGene conversion succeeds, is unavailable, or is deferred.
.breseq_convert_genbank_to_dna <- function(input, output, snapgene_path = NULL,
                                          timeout = 60) {
  valid_path <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  if (!valid_path(input) || !file.exists(input) || dir.exists(input))
    stop("input must name an existing GenBank file.")
  if (!valid_path(output) || !grepl("[.]dna$", output, ignore.case = TRUE))
    stop("output must name a SnapGene .dna file.")
  if (!is.numeric(timeout) || length(timeout) != 1L || !is.finite(timeout) || timeout <= 0)
    stop("timeout must be a positive number of seconds.")
  input <- normalizePath(input, mustWork = TRUE)
  lines <- readLines(input, warn = FALSE)
  if (sum(grepl("^LOCUS[[:space:]]", lines)) != 1L ||
      sum(grepl("^ORIGIN([[:space:]]|$)", lines)) != 1L ||
      sum(grepl("^//[[:space:]]*$", lines)) != 1L)
    stop("SnapGene conversion requires one complete GenBank sequence record.")
  result <- function(status, message, success = FALSE, path = NULL, exit_status = NA_integer_) {
    list(success = success, status = status, output_file = path,
         message = message, exit_status = exit_status)
  }
  if (file.exists(output))
    return(result("output_exists", "The output already exists and was not overwritten."))
  if (!is.null(snapgene_path) &&
      (!valid_path(snapgene_path) || !file.exists(snapgene_path) || dir.exists(snapgene_path)))
    return(result("unavailable", "The configured SnapGene executable is unavailable."))
  executable <- tryCatch(suppressMessages(find_snapgene(snapgene_path)), error = function(e) NULL)
  if (is.null(executable) || file.access(executable, mode = 1L) != 0L)
    return(result("unavailable", "A usable SnapGene executable was not found; the GenBank file is retained."))
  running <- .breseq_snapgene_is_running()
  if (is.na(running))
    return(result("process_check_failed", "Could not determine whether SnapGene is running; conversion was not started."))
  if (running)
    return(result("application_running", "SnapGene is already running; conversion was deferred without closing it."))

  directory <- dirname(path.expand(output))
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
    stop("Could not create output directory: ", directory)
  output <- file.path(normalizePath(directory, mustWork = TRUE), basename(output))
  temporary <- tempfile(".snapgene-conversion-", tmpdir = dirname(output), fileext = ".dna")
  remove_temporary <- TRUE
  on.exit(if (remove_temporary && file.exists(temporary)) unlink(temporary), add = TRUE)
  args <- c("--quiet", "--convert", "SnapGene", "--input", input, "--output", temporary)
  converted <- tryCatch(.breseq_snapgene_run(executable, args, timeout),
    error = function(e) list(unavailable = FALSE, pending = FALSE,
                             status = NA_integer_, message = conditionMessage(e)))
  if (isTRUE(converted$unavailable)) return(result("unavailable", converted$message))
  if (isTRUE(converted$pending)) {
    remove_temporary <- FALSE
    pending <- result("pending", converted$message)
    pending$pending_file <- temporary
    return(pending)
  }
  if (is.na(converted$status) || converted$status != 0L)
    return(result("conversion_failed", converted$message, exit_status = converted$status))
  if (!.breseq_snapgene_dna_header(temporary))
    return(result("invalid_output", "SnapGene did not produce a valid DNA-format header.",
                  exit_status = converted$status))
  if (file.exists(output) || !file.copy(temporary, output, overwrite = FALSE))
    return(result("output_exists", "The destination could not be created without overwriting an existing file."))
  result("converted", "Native SnapGene DNA conversion completed and its binary header was verified.",
         success = TRUE, path = normalizePath(output), exit_status = converted$status)
}

# Bundled: standalone_driver.R
# Internal driver embedded into the generated single-file entry point.
.breseq_standalone_pick <- function(input_dir, pattern, explicit = NULL, optional = FALSE) {
  if (!is.null(explicit)) {
    if (!file.exists(explicit)) stop("Input file does not exist: ", explicit)
    return(normalizePath(explicit))
  }
  candidates <- list.files(input_dir, pattern = pattern, recursive = TRUE,
                           full.names = TRUE, ignore.case = TRUE)
  candidates <- candidates[!grepl("/provenance/|/\\.git/", candidates)]
  if (!length(candidates)) {
    if (optional) return(NULL)
    stop("No matching input file found in ", input_dir, ": ", pattern)
  }
  if (length(candidates) > 1L) {
    hashes <- unname(tools::md5sum(candidates))
    if (length(unique(hashes)) != 1L)
      stop("Multiple different input files match; set the corresponding path explicitly:\n",
           paste(candidates, collapse = "\n"))
    candidates <- candidates[order(nchar(candidates))]
  }
  normalizePath(candidates[1L])
}

.breseq_standalone_sample_id <- function(file) {
  stem <- tools::file_path_sans_ext(basename(file))
  if (tolower(stem) != "index")
    return(sub("_index$", "", sub("^index_", "", stem, ignore.case = TRUE), ignore.case = TRUE))
  parent <- dirname(file)
  if (basename(parent) == "output") parent <- dirname(parent)
  basename(parent)
}

.breseq_standalone_unresolved_coverage <- function(data, samples) {
  empty <- data.frame(sample = character(), seq_id = character(),
    source_start = character(), source_end = character(), gene = character(),
    interpretation = character(), stringsAsFactors = FALSE)
  if (is.null(data) || !nrow(data)) return(empty)
  columns <- .breseq_stage_columns(data, samples)
  seqid <- .breseq_wide_seq(data, columns)
  records <- list()
  for (sample in samples) {
    ends <- .breseq_field(data, columns, "end", sample)
    gene <- .breseq_field(data, columns, "gene", sample)
    present <- which(.breseq_present(ends))
    if (length(present)) records[[length(records) + 1L]] <- data.frame(
      sample = sample, seq_id = seqid[present], source_start = as.character(data$start[present]),
      source_end = ends[present], gene = gene[present],
      interpretation = "Unassigned missing coverage: exact deletion breakpoints are not established; not applied to reconstructed sequence.",
      stringsAsFactors = FALSE)
  }
  if (length(records)) do.call(rbind, records) else empty
}

.breseq_run_standalone <- function(config) {
  required <- c("dplyr", "rvest", "openxlsx", "circlize", "ComplexHeatmap",
                "RColorBrewer", "gtools", "reshape2", "stringr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". Install these packages, then run this same file again.")
  suppressPackageStartupMessages(library(dplyr))
  if (getRversion() < "4.1.0") stop("This script requires R 4.1 or newer.")
  input <- normalizePath(path.expand(config$input_dir), mustWork = TRUE)
  output <- path.expand(config$output_dir)
  if (!grepl("^/", output)) output <- file.path(getwd(), output)
  check_output <- file.path(normalizePath(dirname(output), mustWork = FALSE), basename(output))
  if (check_output == input || startsWith(check_output, paste0(input, "/")))
    stop("output_dir must be outside input_dir so original inputs remain unchanged.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  if (output == input || startsWith(output, paste0(input, "/")))
    stop("output_dir must be outside input_dir so original inputs remain unchanged.")
  message("Input: ", input, "\nOutput: ", output)
  work <- file.path(output, "provenance")
  dir.create(work, showWarnings = FALSE)

  reference <- .breseq_standalone_pick(input,
    "\\.(gbff|gbk|gb|genbank|fasta|fna|fa)$", config$reference_file)
  refs <- .breseq_prepare_standalone_reference(reference, file.path(work, "reference"))
  gff <- .breseq_standalone_pick(input, "\\.gff3?$", config$gff_file, optional = TRUE)
  eggnog <- .breseq_standalone_pick(input, "emapper[.]annotations([.](xlsx|tsv))?$",
                                   config$eggnog_file, optional = TRUE)
  gff_data <- if (is.null(gff)) NULL else parse_gff_file(gff)
  eggnog_data <- if (is.null(eggnog)) NULL else parse_eggnog_file(eggnog)
  html <- list.files(input, pattern = "index.*[.]html$", full.names = TRUE,
                     recursive = TRUE, ignore.case = TRUE)
  html <- html[!grepl("/provenance/|/\\.git/", html)]
  source_files <- c(reference, gff, eggnog)
  baseline <- config$baseline_id
  order <- config$sample_order
  if (is.null(config$excel_file) && length(html)) {
    ids <- vapply(html, .breseq_standalone_sample_id, character(1))
    if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("HTML sample IDs are empty or duplicated.")
    tables <- setNames(lapply(html, extract_single_html), ids)
    if (is.null(order)) {
      order <- ids[gtools::mixedorder(ids)]
      if (!is.null(baseline) && !baseline %in% order) order <- c(baseline, order)
    }
    if (anyNA(order) || anyDuplicated(order) || !all(ids %in% order) ||
        any(!order %in% c(ids, baseline)))
      stop("sample_order must contain each HTML sample once and, optionally, baseline_id.")
    if (!is.null(baseline) && baseline %in% order && !baseline %in% ids)
      tables[[baseline]] <- lapply(tables[[1L]], function(d) d[FALSE, , drop = FALSE])
    tables <- tables[order]
    if (length(tables) < 2L) stop("A single sample requires a distinct baseline_id.")
    merged <- merge_comparison_results(compare_multiple_tables(tables))
    mut_name <- names(merged)[vapply(names(merged), function(n) .breseq_table_kind(n) == "mutation", logical(1))]
    if (length(mut_name)) merged[[mut_name[1L]]] <- add_mutation_types(merged[[mut_name[1L]]])
    mutation_data <- setNames(merged, c("Mutations", "Missing", "New junction")[match(names(merged),
      c("Mutation predictions", "Unassigned missing coverage evidence", "Unassigned new junction evidence"))])
    source_files <- c(source_files, html)
  } else {
    excel <- .breseq_standalone_pick(input, "(merged_mutation_results|mutant_results)[.]xlsx$",
                                    config$excel_file)
    mutation_data <- load_mutation_excel(excel)
    if (is.null(mutation_data$Mutations)) stop("No canonical Mutations worksheet in ", excel)
    if (is.null(order)) {
      order <- attr(mutation_data$Mutations, "breseq_sample_order", exact = TRUE)
      if (is.null(order)) order <- attr(.breseq_stage_columns(mutation_data$Mutations), "sample_order")
    }
    merged <- setNames(mutation_data,
      c("Mutation predictions", "Unassigned missing coverage evidence", "Unassigned new junction evidence")[match(names(mutation_data), c("Mutations", "Missing", "New junction"))])
    source_files <- c(source_files, excel)
  }
  labels <- .breseq_validate_labels(config$sample_labels, order)
  excel_out <- file.path(output, "merged_mutation_results.xlsx")
  wb <- openxlsx::createWorkbook()
  sheet_map <- c("Mutation predictions" = "Mutations",
    "Unassigned missing coverage evidence" = "Missing",
    "Unassigned new junction evidence" = "New junction")
  for (name in names(merged)) {
    sheet <- unname(sheet_map[name])
    if (is.na(sheet)) sheet <- substr(name, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    highlight <- .breseq_table_kind(name) != "junction" ||
      !is.null(attr(merged[[name]], "breseq_keys")) || "comparison_key" %in% names(merged[[name]])
    if (!highlight) message("Junction pair identities are unavailable in the imported workbook; preserving its values without inferring stage event highlights.")
    .breseq_write_presentation(wb, sheet, merged[[name]], name,
      highlight_new = highlight, clean_headers = TRUE, sample_labels = labels,
      sample_order = order, mutation_data = mutation_data$Mutations)
  }
  openxlsx::saveWorkbook(wb, excel_out, overwrite = TRUE)
  processed <- process_mutation_data_for_viz(mutation_data, sample_order = order)
  # A workbook may retain region annotations in a sample with no variant call.
  # Annotation/replicon metadata alone must not become an UNKNOWN mutation.
  calls <- lapply(processed$mutation_list, function(d) {
    d[.breseq_present(.breseq_text(d$mutation)), , drop = FALSE]
  })
  raw_calls <- calls
  raw_counts <- setNames(vapply(order, function(s) if (is.null(calls[[s]])) 0L else nrow(calls[[s]]), integer(1)), order)

  display_file <- config$display_deletions_file
  if (is.null(display_file)) {
    candidates <- list.files(input, pattern = "^coverage_derived_DEL_rows[.]tsv$",
                             recursive = TRUE, full.names = TRUE)
    if (length(candidates) > 1L && length(unique(unname(tools::md5sum(candidates)))) != 1L)
      stop("Several different coverage display tables exist; set display_deletions_file.")
    if (length(candidates)) display_file <- candidates[order(nchar(candidates))][1L]
  }
  display_rows <- NULL
  if (!is.null(display_file)) {
    aliases <- config$sample_aliases
    if (is.null(aliases)) {
      short_ids <- sub("^[0-9]+(SJY)?", "", names(calls), ignore.case = TRUE, perl = TRUE)
      if (anyDuplicated(short_ids))
        stop("Coverage evidence IDs are ambiguous; set SAMPLE_ALIASES explicitly.")
      aliases <- stats::setNames(names(calls), short_ids)
    }
    display <- .breseq_standalone_display_deletions(calls, display_file,
      aliases, refs$info, merged_missing = mutation_data$Missing)
    calls <- display$mutation_list
    display_rows <- display$display_rows
    utils::write.table(display_rows, file.path(work, "display_deletion_evidence.tsv"),
                       sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
    source_files <- c(source_files, display_file)
  }
  plot_counts <- setNames(vapply(order, function(s) if (is.null(calls[[s]])) 0L else nrow(calls[[s]]), integer(1)), order)
  counts <- data.frame(sample = order, label = unname(labels[order]),
    raw_calls = unname(raw_counts), plot_calls = unname(plot_counts))
  utils::write.table(counts, file.path(output, "sample_call_counts.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  genome_export <- NULL
  create_genbank <- !identical(config$create_genbank, FALSE)
  if (create_genbank) {
    losses <- NULL
    if (!is.null(display_rows)) {
      loss_at <- display_rows$evidence_origin == "COVERAGE_RATIO" & display_rows$appended
      losses <- display_rows[loss_at, , drop = FALSE]
      if (nrow(losses)) losses$sample <- losses$sample_id
    }
    genome_export <- export_changed_sequences_to_genbank(reference, raw_calls,
      file.path(output, "changed_sequences"), samples = order, sample_labels = labels,
      replicon_losses = losses)
    for (artifact in c("feature_audit", "edits", "call_report")) {
      table <- genome_export[[artifact]]
      if (is.data.frame(table) && ncol(table)) utils::write.table(table,
        file.path(output, "changed_sequences", paste0(artifact, ".tsv")),
        sep = "\t", row.names = FALSE, quote = TRUE, na = "NA")
    }
    if (any(genome_export$manifest$status == "unresolved_calls")) {
      utils::write.table(genome_export$manifest,
        file.path(output, "changed_sequences", "sequence_manifest.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
      stop("Some finalized calls could not be applied exactly. See changed_sequences/call_report.tsv and sequence_manifest.tsv; DNA conversion was not started.")
    }
    unresolved <- .breseq_standalone_unresolved_coverage(mutation_data$Missing, order)
    if (nrow(unresolved)) {
        utils::write.table(unresolved, file.path(output, "changed_sequences", "unresolved_coverage_regions.tsv"),
          sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
        key <- paste(unresolved$sample, .breseq_seq(unresolved$seq_id), sep = "\r")
        manifest_key <- paste(genome_export$manifest$sample, .breseq_seq(genome_export$manifest$seq_id), sep = "\r")
        genome_export$manifest$unresolved_coverage <- manifest_key %in% key
    }
    if (!"unresolved_coverage" %in% names(genome_export$manifest))
      genome_export$manifest$unresolved_coverage <- FALSE
    genome_export$manifest$sequence_basis <- ifelse(genome_export$manifest$unresolved_coverage,
      "Finalized calls only; unresolved missing coverage was not applied.",
      "Finalized calls relative to the supplied reference.")
    genome_export$manifest$dna_file <- NA_character_
    genome_export$manifest$dna_status <- "not_requested"
    if (!identical(config$create_dna, FALSE)) {
      for (i in which(!is.na(genome_export$manifest$genbank_file) & nzchar(genome_export$manifest$genbank_file))) {
        gbk <- genome_export$manifest$genbank_file[i]
        dna <- sub("[.]gb[k]?$", ".dna", gbk, ignore.case = TRUE)
        converted <- .breseq_convert_genbank_to_dna(gbk, dna, config$snapgene_path)
        genome_export$manifest$dna_status[i] <- converted$status
        if (isTRUE(converted$success)) genome_export$manifest$dna_file[i] <- converted$output_file
        else message("SnapGene export ", converted$status, ": ", converted$message)
      }
    }
    utils::write.table(genome_export$manifest,
      file.path(output, "changed_sequences", "sequence_manifest.tsv"), sep = "\t",
      row.names = FALSE, quote = FALSE, na = "NA")
  }
  plots <- character()
  plot_dir <- file.path(output, "plots")
  dir.create(plot_dir, showWarnings = FALSE)
  ref_ids <- .breseq_seq(refs$info$seq_id)
  known <- unique(unlist(lapply(calls, function(d) if (nrow(d)) .breseq_seq(d[["seq id"]]) else character())))
  if (any(!known %in% ref_ids)) stop("Mutation replicons missing from the reference: ", paste(setdiff(known, ref_ids), collapse = ", "))
  for (rep in refs$replicons) {
    seqid <- .breseq_seq(rep$genome_info$seq_id[1L])
    per_calls <- lapply(calls, function(d) d[.breseq_seq(d[["seq id"]]) == seqid, , drop = FALSE])
    per_gff <- if (is.null(gff_data)) NULL else gff_data[.breseq_seq(gff_data$seqid) == seqid, , drop = FALSE]
    gc <- calculate_gc_content(rep$reference_file,
      window_size = min(10000L, rep$genome_info$length[1L]), gc_colors = config$gc_colors)
    if (is.null(gc)) stop("Could not calculate GC statistics for ", seqid)
    filename <- file.path(plot_dir, paste0("circos_", gsub("[^[:alnum:]_-]", "_", seqid), ".pdf"))
    plot_calls <- setNames(per_calls, unname(labels[names(per_calls)]))
    create_full_circos_plot(mutation_list = plot_calls, missing_data = NULL,
      samples = unname(labels[order]), genome_info = rep$genome_info, gff_data = per_gff,
      eggnog_data = eggnog_data, gc_data = gc, output_file = filename,
      hide_intergenic_labels = config$hide_intergenic_labels, gc_colors = config$gc_colors,
      observed_tracks = TRUE, reverse_tracks = FALSE)
    if (!file.exists(filename) || file.info(filename)$size == 0) stop("Plot export failed: ", filename)
    plots[seqid] <- filename
    saveRDS(gc, file.path(work, paste0("gc_", seqid, ".rds")))
  }
  source_files <- unique(normalizePath(source_files, mustWork = TRUE))
  manifest <- data.frame(path = source_files, md5 = unname(tools::md5sum(source_files)))
  utils::write.table(manifest, file.path(work, "source_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  result <- list(output_dir = output, excel_file = excel_out, plots = plots,
    samples = order, sample_labels = labels, counts = counts,
    mutation_list = calls, merged_results = merged, display_rows = display_rows,
    reference_info = refs$info, genome_export = genome_export)
  saveRDS(result, file.path(work, "one_shot_results.rds"))
  message("Complete: ", output, "\nCreated ", length(plots), " replicon plots and the combined Excel workbook.")
  if (!is.null(genome_export)) message("Changed-sequence files: ", file.path(output, "changed_sequences"))
  invisible(result)
}

# Bundled: standalone_evidence.R
# Add explicitly audited display evidence without changing canonical call rows.
# The sidecar is an input: MC cores are not inferred again from boundary depths.
.breseq_standalone_display_deletions <- function(mutation_list, evidence_file,
                                                sample_aliases, reference_info,
                                                merged_missing = NULL) {
  if (!is.character(evidence_file) || length(evidence_file) != 1L ||
      is.na(evidence_file) || !file.exists(evidence_file))
    stop("evidence_file must name an existing audited TSV file.")
  evidence <- utils::read.delim(evidence_file, check.names = FALSE,
                              stringsAsFactors = FALSE)
  required <- c("sample", "seq_id", "start", "end", "size", "evidence_origin",
                "source_start", "source_end", "depth_ratio", "criterion",
                "raw_final_DEL_call")
  if (!all(required %in% names(evidence)))
    stop("Audited evidence TSV is missing required columns: ",
         paste(setdiff(required, names(evidence)), collapse = ", "))
  evidence$sample_id <- character(nrow(evidence))
  evidence$appended <- rep(FALSE, nrow(evidence))
  evidence$skip_reason <- rep("", nrow(evidence))
  if (!nrow(evidence)) return(list(mutation_list = mutation_list, display_rows = evidence))

  samples <- names(mutation_list)
  if (!is.list(mutation_list) || is.null(samples) || anyNA(samples) ||
      any(!nzchar(samples)) || anyDuplicated(samples))
    stop("mutation_list must have unique, non-empty sample IDs in analysis order.")
  if (!is.character(sample_aliases) || is.null(names(sample_aliases)) ||
      anyNA(sample_aliases) || anyNA(names(sample_aliases)) ||
      any(!nzchar(sample_aliases)) || any(!nzchar(names(sample_aliases))) ||
      anyDuplicated(names(sample_aliases)))
    stop("sample_aliases must map unique source aliases to actual sample IDs.")
  source_samples <- .breseq_text(evidence$sample)
  mapped <- source_samples
  alias_at <- match(source_samples, names(sample_aliases))
  mapped[!is.na(alias_at)] <- unname(sample_aliases[alias_at[!is.na(alias_at)]])
  if (anyNA(mapped) || any(!mapped %in% samples))
    stop("Audited evidence contains an unknown sample or sample alias.")
  evidence$sample_id <- mapped

  if (!all(c("seq_id", "length") %in% names(reference_info)))
    stop("reference_info must contain seq_id and length.")
  reference_ids <- .breseq_seq(reference_info$seq_id)
  reference_lengths <- suppressWarnings(as.numeric(reference_info$length))
  if (anyNA(reference_ids) || any(!nzchar(reference_ids)) || anyDuplicated(reference_ids) ||
      any(!is.finite(reference_lengths) | reference_lengths < 1 |
          reference_lengths != floor(reference_lengths)))
    stop("reference_info must have unique replicon IDs and positive integer lengths.")
  seqids <- .breseq_seq(evidence$seq_id)
  reference_at <- match(seqids, reference_ids)
  if (anyNA(reference_at)) stop("Audited evidence contains an unknown reference replicon.")
  coordinates <- lapply(evidence[c("start", "end", "size")], function(x)
    suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE))))
  start <- coordinates$start; end <- coordinates$end; size <- coordinates$size
  if (any(!is.finite(c(start, end, size))) ||
      any(start != floor(start) | end != floor(end) | size != floor(size)) ||
      any(start < 1 | end < start | size < 1 | end > reference_lengths[reference_at]))
    stop("Audited evidence coordinates must be positive, within-reference integer intervals.")
  if (any(size != end - start + 1)) stop("Audited evidence size must equal end - start + 1.")
  origin <- .breseq_text(evidence$evidence_origin)
  if (anyNA(origin) || any(!origin %in% c("MC", "COVERAGE_RATIO")))
    stop("Audited evidence_origin must be MC or COVERAGE_RATIO.")
  final_call <- toupper(.breseq_text(evidence$raw_final_DEL_call))
  if (anyNA(final_call) || any(final_call != "FALSE"))
    stop("Audited display evidence must have raw_final_DEL_call = FALSE.")
  ratio <- suppressWarnings(as.numeric(evidence$depth_ratio))
  coverage <- origin == "COVERAGE_RATIO"
  if (any(!is.finite(ratio[coverage]) | ratio[coverage] < 0 | ratio[coverage] >= 0.1) ||
      any(start[coverage] != 1 | end[coverage] != reference_lengths[reference_at[coverage]]))
    stop("COVERAGE_RATIO evidence requires depth_ratio <0.1 and the whole replicon interval.")
  if (anyDuplicated(paste(mapped, seqids, start, end, sep = "\r")))
    stop("Audited evidence contains duplicate sample/replicon intervals.")

  missing_columns <- if (!is.null(merged_missing)) .breseq_stage_columns(merged_missing) else NULL
  missing_seqids <- if (!is.null(missing_columns)) .breseq_wide_seq(merged_missing, missing_columns) else NULL
  source_keys <- function(data) {
    if (!nrow(data)) return(character())
    seq_columns <- intersect(c("seq id", "seq_id", "seqid"), names(data))
    if (!length(seq_columns) || !all(c("position", "mutation") %in% names(data)))
      stop("Mutation rows require a replicon ID, position, and mutation.")
    paste(.breseq_seq(data[[seq_columns[1]]]), .breseq_coord(data$position),
          .breseq_allele(data$mutation), sep = "\r")
  }
  seen <- character()
  for (sample in samples) {
    data <- mutation_list[[sample]]
    if (!is.data.frame(data)) stop("Every mutation_list entry must be a data frame.")
    automatic_rownames <- .row_names_info(data, 1L) < 0L
    # Prior source calls and prior display rows determine first-appearance labels.
    prior_keys <- seen
    source_keys(data)
    for (i in which(mapped == sample)) {
      seq_columns <- intersect(c("seq id", "seq_id", "seqid"), names(data))
      needed <- c("position", "mutation", "annotation", "evidence", "gene",
                  "description", "new", "ALE_version", "version_color", "label")
      if (!length(seq_columns) || !all(needed %in% names(data)))
        stop("Display rows require the standard processed mutation columns.")
      if (any(vapply(data, is.factor, logical(1))))
        stop("Processed mutation columns must not be factors when appending display evidence.")
      same_replicon <- .breseq_seq(data[[seq_columns[1]]]) == seqids[i]
      deleted_size <- suppressWarnings(as.numeric(sub("^[\u0394\u2206]([0-9,]+).*", "\\1",
        .breseq_allele(data$mutation))))
      deleted_start <- suppressWarnings(as.numeric(.breseq_coord(data$position)))
      final_del <- same_replicon & grepl("^[\u0394\u2206][0-9,]+bp$", .breseq_allele(data$mutation)) &
        is.finite(deleted_start) & is.finite(deleted_size)
      deleted_end <- deleted_start + deleted_size - 1
      if (any(final_del & deleted_start <= start[i] & deleted_end >= end[i], na.rm = TRUE)) {
        evidence$skip_reason[i] <- "already covered by same-sample DEL"
        next
      }
      if (any(final_del & deleted_start <= end[i] & deleted_end >= start[i], na.rm = TRUE))
        stop("Audited display interval partially overlaps an existing same-sample DEL.")

      gene <- if (coverage[i]) paste0("p", seqids[i]) else ""
      if (!coverage[i] && !is.null(missing_columns)) {
        missing_end <- .breseq_field(merged_missing, missing_columns, "end", sample)
        at <- which(missing_seqids == seqids[i] &
          .breseq_coord(merged_missing$start) == .breseq_coord(evidence$source_start[i]) &
          .breseq_coord(missing_end) == .breseq_coord(evidence$source_end[i]))
        if (length(at) > 1L) stop("Ambiguous source MC row for an audited display interval.")
        if (length(at)) {
          source_gene <- .breseq_field(merged_missing, missing_columns, "gene", sample)[at]
          if (.breseq_present(source_gene)) gene <- source_gene
        }
      }
      row <- data[NA_integer_, , drop = FALSE]
      for (name in names(row)) {
        row[[name]] <- if (is.logical(data[[name]])) FALSE else
          if (is.integer(data[[name]])) NA_integer_ else
            if (is.numeric(data[[name]])) NA_real_ else ""
      }
      for (name in seq_columns) row[[name]] <- seqids[i]
      row$position <- if (is.integer(data$position)) as.integer(start[i]) else
        if (is.character(data$position)) as.character(start[i]) else start[i]
      row$mutation <- paste0("\u0394", format(size[i], big.mark = ",", scientific = FALSE, trim = TRUE), " bp")
      row$annotation <- "coverage-derived DEL display"
      row$evidence <- origin[i]
      row$gene <- gene
      row$description <- if (coverage[i]) "Whole-replicon depth / chromosome depth <0.1" else
        "Conservative internal span of source missing-coverage evidence"
      row$new <- !source_keys(row) %in% c(prior_keys, source_keys(data))
      row$ALE_version <- sample
      colors <- data$version_color[.breseq_present(data$version_color)]
      if (!length(colors)) stop("A display row requires the sample's existing version_color.")
      row$version_color <- colors[1]
      row$label <- trimws(paste(row$mutation, gene))
      if ("start" %in% names(row)) row$start <- if (is.integer(data$start)) as.integer(start[i]) else start[i]
      if ("end" %in% names(row)) row$end <- if (is.integer(data$end)) as.integer(end[i]) else end[i]
      if ("type" %in% names(row)) row$type <- "DEL"
      if ("mutation_type" %in% names(row)) row$mutation_type <- if (size[i] >= 50) "large_deletion" else "deletion"
      source_row_names <- attr(data, "row.names")
      rownames(row) <- make.unique(c(rownames(data), paste0("display_", i)))[nrow(data) + 1L]
      data[nrow(data) + 1L, ] <- row
      if (automatic_rownames) rownames(data) <- NULL else if (is.integer(source_row_names))
        attr(data, "row.names") <- c(source_row_names, max(c(0L, source_row_names)) + 1L)
      evidence$appended[i] <- TRUE
    }
    mutation_list[[sample]] <- data
    seen <- union(seen, source_keys(data))
  }
  list(mutation_list = mutation_list, display_rows = evidence)
}

# Bundled: standalone_reference.R
# Prepare one reference per replicon for the standalone driver.  The GC reader
# consumes one record at a time, so the original multi-record file must not be
# reused for every replicon.  This helper deliberately does not parse features;
# the driver can supply separately parsed GFF annotations through gff_data.
#
# GenBank records retain their sequence and feature text.  Only the temporary
# LOCUS identifier is aligned to VERSION (ACCESSION/LOCUS fallback), with a final
# numeric accession version removed to match breseq/GFF sequence identifiers.
# FASTA records retain their sequence and description with the same ID rule.
# Explicit topology is retained.  FASTA without topology defaults to circular,
# as in parse_fasta_info(); GenBank without topology remains unknown.
.breseq_prepare_standalone_reference <- function(reference_file, work_dir) {
  scalar_path <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!scalar_path(reference_file) || !file.exists(reference_file) ||
      dir.exists(reference_file)) {
    stop("Reference file not found: ", paste(reference_file, collapse = ", "))
  }
  if (!scalar_path(work_dir)) stop("work_dir must be a nonempty directory path")

  extension <- tolower(tools::file_ext(reference_file))
  genbank <- extension %in% c("gb", "gbk", "gbff", "genbank")
  fasta <- extension %in% c("fa", "fna", "fasta")
  if (!genbank && !fasta) {
    stop("Standalone references must be GenBank (.gb, .gbk, .gbff, .genbank) ",
         "or nucleotide FASTA (.fa, .fna, .fasta).")
  }

  lines <- readLines(reference_file, warn = FALSE)
  starts <- grep(if (genbank) "^LOCUS[[:space:]]" else "^>", lines)
  if (!length(starts)) {
    stop(if (genbank) "No LOCUS records found in GenBank reference" else
           "No FASTA records found in reference")
  }
  ends <- c(starts[-1L] - 1L, length(lines))
  normalize_id <- function(id) sub("\\.[0-9]+$", "", id)
  first_field <- function(text) {
    fields <- strsplit(trimws(text), "[[:space:]]+")[[1L]]
    if (length(fields) && nzchar(fields[1L])) fields[1L] else ""
  }
  topology_from <- function(text, default) {
    if (grepl("\\bcircular\\b", text, ignore.case = TRUE)) "circular" else
      if (grepl("\\blinear\\b", text, ignore.case = TRUE)) "linear" else default
  }

  records <- lapply(seq_along(starts), function(i) {
    record <- lines[seq.int(starts[i], ends[i])]
    if (genbank) {
      terminator <- grep("^//[[:space:]]*$", record)
      if (length(terminator) != 1L) {
        stop("GenBank record ", i, " must have exactly one // terminator")
      }
      record <- record[seq_len(terminator)]
      origin <- grep("^ORIGIN([[:space:]]|$)", record)
      if (length(origin) != 1L || origin >= length(record) - 1L) {
        stop("GenBank record ", i, " has no usable ORIGIN sequence")
      }
      sequence <- toupper(gsub("[[:space:][:digit:]]", "",
        paste(record[seq.int(origin + 1L, length(record) - 1L)], collapse = "")))
      locus_id <- first_field(sub("^LOCUS[[:space:]]+", "", record[1L]))
      field_value <- function(field) {
        hit <- grep(paste0("^", field, "[[:space:]]+"), record)
        if (length(hit)) first_field(sub(paste0("^", field, "[[:space:]]+"),
                                        "", record[hit[1L]])) else ""
      }
      id_candidates <- c(field_value("VERSION"), field_value("ACCESSION"), locus_id)
      id <- normalize_id(id_candidates[nzchar(id_candidates)][1L])
      declared_match <- regmatches(record[1L],
        regexpr("[0-9]+[[:space:]]+bp\\b", record[1L], ignore.case = TRUE))
      declared_length <- if (length(declared_match)) {
        as.numeric(sub("[[:space:]]+bp$", "", declared_match, ignore.case = TRUE))
      } else NA_real_
      if (!is.finite(declared_length) || declared_length != nchar(sequence)) {
        stop("GenBank record ", i, " LOCUS length does not match its ORIGIN sequence")
      }
      topology <- topology_from(record[1L], "unknown")
      definition <- grep("^DEFINITION[[:space:]]", record)
      description <- ""
      if (length(definition)) {
        j <- definition[1L]
        description <- sub("^DEFINITION[[:space:]]+", "", record[j])
        while (j + 1L < origin && grepl("^ {12}[^ ]", record[j + 1L])) {
          j <- j + 1L
          description <- paste(description, trimws(record[j]))
        }
      }
      # Preserve all other record text, including feature locations and bases.
      record[1L] <- paste0("LOCUS       ", id,
        sub("^LOCUS[[:space:]]+[^[:space:]]+", "", record[1L]))
    } else {
      header <- trimws(sub("^>", "", record[1L]))
      raw_id <- first_field(header)
      id <- normalize_id(raw_id)
      description <- trimws(sub("^[^[:space:]]+", "", header))
      if (length(record) < 2L) stop("FASTA record ", i, " has no sequence")
      sequence <- toupper(gsub("[[:space:]]", "", paste(record[-1L], collapse = "")))
      topology <- topology_from(description, "circular")
      # Normalizing whitespace does not remove ambiguous nucleotide symbols.
      chunk_starts <- seq.int(1L, max(1L, nchar(sequence)), by = 80L)
      record <- c(paste0(">", id, if (nzchar(description)) paste0(" ", description)),
                  substring(sequence, chunk_starts, chunk_starts + 79L))
    }
    if (is.na(id) || !nzchar(id)) stop("Reference record ", i, " has no sequence ID")
    if (!nzchar(sequence) || grepl("[^ACGTRYSWKMBDHVN]", sequence)) {
      stop("Reference record ", i, " must contain a nonempty nucleotide sequence")
    }
    list(lines = record, info = data.frame(seq_id = id,
      length = nchar(sequence), topology = topology, description = description,
      stringsAsFactors = FALSE))
  })

  info <- do.call(rbind, lapply(records, `[[`, "info"))
  rownames(info) <- NULL
  if (anyDuplicated(info$seq_id)) {
    stop("Reference has duplicate sequence IDs after removing accession versions: ",
         paste(unique(info$seq_id[duplicated(info$seq_id)]), collapse = ", "))
  }
  if (!dir.exists(work_dir) && !dir.create(work_dir, recursive = TRUE)) {
    stop("Could not create reference work directory: ", work_dir)
  }
  work_dir <- normalizePath(work_dir, mustWork = TRUE)
  written <- character()
  complete <- FALSE
  on.exit(if (!complete) unlink(written), add = TRUE)
  replicons <- lapply(seq_along(records), function(i) {
    safe_id <- gsub("[^A-Za-z0-9_.-]", "_", info$seq_id[i])
    destination <- tempfile(paste0(sprintf("%03d", i), "_", safe_id, "_"),
                            tmpdir = work_dir, fileext = if (genbank) ".gbff" else ".fasta")
    written <<- c(written, destination)
    writeLines(records[[i]]$lines, destination, useBytes = TRUE)
    list(reference_file = destination, genome_info = info[i, , drop = FALSE],
         gff_data = NULL)
  })
  names(replicons) <- info$seq_id
  complete <- TRUE
  list(info = info, replicons = replicons)
}

  function(config) .breseq_run_standalone(config)
})
if (!identical(getOption('breseqConverter.autorun', TRUE), FALSE)) {
  .breseq_overrides <- getOption('breseqConverter.one_shot.config', list())
  ONE_SHOT_RESULTS <- breseq_one_shot(utils::modifyList(.breseq_one_shot_config, .breseq_overrides, keep.null=TRUE))
}
