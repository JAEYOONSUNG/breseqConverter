one_shot_definitions <- function() {
  file <- system.file("scripts", "one_shot_legacy.R", package = "breseqConverter")
  if (!nzchar(file)) file <- test_path("..", "..", "inst", "scripts", "one_shot_legacy.R")
  expressions <- parse(file)
  env <- new.env(parent = baseenv())
  # Load only settings and function definitions, never workbook reads or drawing.
  for (expr in expressions) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        (as.character(expr[[2]]) %in% c("HIDE_INTERGENIC_LABELS", "GC_COLORS") ||
         (is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))))) {
      eval(expr, env)
    }
  }
  env
}

one_shot_find_calls <- function(expr, predicate) {
  found <- list()
  visit <- function(x) {
    if (missing(x) || !is.call(x)) return(invisible(NULL))
    if (predicate(x)) found[[length(found) + 1L]] <<- x else
      for (part in as.list(x)[-1L]) visit(part)
  }
  visit(expr)
  found
}

test_that("standalone options work without executing the legacy workflow", {
  env <- one_shot_definitions()
  expect_false(exists("breseq_result_dir", env, inherits = FALSE))
  expect_false(env$HIDE_INTERGENIC_LABELS)
  expect_identical(env$GC_COLORS, .breseq_gc_colors())
  expect_identical(eval(formals(env$breseqVisualization)$gc_colors, env), env$GC_COLORS)
  env$HIDE_INTERGENIC_LABELS <- TRUE
  expect_true(eval(formals(env$breseqVisualization)$hide_intergenic_labels, env))
  expect_identical(env$.breseq_one_shot_gc_track_colors("gc_skew", c(-1, 1, 0, NA)),
                   c("#FF0000", "#00FF00", "#00FF00", NA))
  expect_identical(env$.breseq_one_shot_gc_track_colors("gc_ratio", c(-1, 1)),
                   c("#FFFF00", "#0000FF"))
  expect_identical(env$.breseq_one_shot_gc_colors(c(skew_positive = "#123456")),
                   .breseq_gc_colors(c(skew_positive = "#123456")))
  for (invalid in list(character(), c(skew_positive = ""), "green",
                       c(unknown = "red"), c(skew_positive = NA_character_),
                       c(skew_positive = "red", skew_positive = "blue"),
                       c(ratio_negative = "not-a-color"))) {
    expect_error(env$.breseq_one_shot_gc_colors(invalid), "gc_colors")
  }
})

test_that("the standalone label path filters labels while retaining all source marks", {
  skip_if_not_installed("dplyr")
  env <- one_shot_definitions()
  env$`%>%` <- getExportedValue("dplyr", "%>%")
  data <- data.frame(
    seqid = "chr", start = 1:6, end = 1:6,
    label = paste0("label", 1:6), version_color = "black", new = TRUE,
    annotation = c("intergenic (+1/-2)", " \tINTERGENIC\u00a0(+3/-4)",
                   "coding (1/3 nt)", "intergenic-like protein", NA, "pseudogene"),
    stringsAsFactors = FALSE
  )
  env$mut_data <- data
  env$mm_h <- identity
  env$circos.genomicLabels <- function(bed, ...) { env$captured_labels <- bed }
  calls <- one_shot_find_calls(body(env$breseqVisualization), function(x) {
    (identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("label_data"))) ||
      identical(x[[1]], as.name("circos.genomicLabels"))
  })
  for (hide in c(FALSE, TRUE)) {
    env$hide_intergenic_labels <- hide
    for (call in calls) eval(call, env)
    expect_identical(env$captured_labels$start, if (hide) 3:6 else 1:6)
    expect_identical(env$mut_data, data)
  }
  without_annotation <- data[, setdiff(names(data), "annotation")]
  expect_identical(env$.breseq_one_shot_filter_intergenic_labels(without_annotation, TRUE),
                   without_annotation)
})

test_that("legacy and standalone GC rendering ignores stale cached colors", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("dplyr")
  env <- one_shot_definitions()
  env$`%>%` <- getExportedValue("dplyr", "%>%")
  fixture <- data.frame(seqid = "chr", start = c(1, 10001), end = c(10000, 20000),
    gc_skew_minus_average = c(-0.2, 0.2), gc_ratio_minus_average = c(0.1, -0.1),
    gc_skew_color = "magenta", gc_ratio_color = "cyan")
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(
    circos.track = function(panel.fun, ...) panel.fun(NULL, NULL),
    circos.barplot = function(value, pos, col, ...) { captured$plot <- list(value = value, col = col) },
    get.cell.meta.data = function(...) "chr",
    get.current.track.index = function() 0,
    .package = "circlize"
  )
  for (renderer in list(env$breseqVisualization, breseqVisualization)) {
    context <- new.env(parent = env)
    context$nucleotide <- c(chr = "unused")
    context$i <- 1L
    context$seqid <- "chr"
    context$window_chr <- fixture
    context$gc_palette <- .breseq_gc_colors(c(skew_positive = "#123456"))
    assignment <- one_shot_find_calls(body(renderer), function(x) {
      identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("plot_gc_track"))
    })
    expect_length(assignment, 1L)
    eval(assignment[[1]], context)
    for (type in c("gc_skew", "gc_ratio")) {
      context$plot_gc_track(type)
      expect_identical(captured$plot$value, fixture[[paste0(type, "_minus_average")]])
      expect_identical(captured$plot$col, if (type == "gc_skew")
        c("#FF0000", "#123456") else c("#0000FF", "#FFFF00"))
      expect_identical(context$window_chr, fixture)
    }
  }
})
