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