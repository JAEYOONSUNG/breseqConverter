#' breseqVisualization
#'
#' @param fna_count Number of contigs to process. If NULL, uses the number of rows in contig_length.
#' @param completeness Logical indicating if the genome is complete. If TRUE, saves output as PDF after interactive display.
#' @param mutation_data A data frame with mutation data (columns: start, end, type, color, ALE_version; or position, mutation, color, ALE_version for raw data).
#' @return A genome map displayed interactively in the plot window, with optional PDF output.
#' @export


breseqVisualization <- function(fna_count = NULL, completeness = NULL, mutation_data = NULL) {
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
          dplyr::filter(!is.na(label) & label != "") %>%
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
                                     circlize::circos.barplot(
                                       value = current_sector_data[[paste0(type, "_minus_average")]] %>% as.numeric(),
                                       pos = (current_sector_data$start + current_sector_data$end) / 2,
                                       col = current_sector_data[[paste0(type, "_color")]],
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
