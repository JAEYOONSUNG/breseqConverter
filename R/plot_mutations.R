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
