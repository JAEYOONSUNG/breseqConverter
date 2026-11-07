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
