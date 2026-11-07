test_that("GC palette supports validated partial overrides", {
  defaults <- c(skew_positive = "#00FF00", skew_negative = "#FF0000",
                ratio_positive = "#0000FF", ratio_negative = "#FFFF00")
  expect_identical(.breseq_gc_colors(), defaults)
  custom <- .breseq_gc_colors(c(skew_positive = "darkgreen", ratio_negative = "#ABCDEF"))
  expect_identical(custom[c("skew_positive", "ratio_negative")],
                   c(skew_positive = "darkgreen", ratio_negative = "#ABCDEF"))
  expect_identical(custom[c("skew_negative", "ratio_positive")],
                   defaults[c("skew_negative", "ratio_positive")])
  invalid <- list("red", c(unknown = "red"), c(skew_positive = "red", skew_positive = "blue"),
                  c(skew_positive = NA_character_), setNames("red", NA_character_),
                  setNames("red", ""), c(skew_positive = "not-a-color"),
                  c(skew_positive = "#12GG34"), c(skew_positive = ""), c(skew_positive = 1))
  for (colors in invalid) expect_error(.breseq_gc_colors(colors), "gc_colors")
})

test_that("palette changes do not change GC calculation, centering or sign", {
  reference <- tempfile(fileext = ".fasta")
  on.exit(unlink(reference), add = TRUE)
  writeLines(c(">fixture", "CCCCGGGAAATT"), reference)
  actual <- suppressMessages(calculate_gc_content(reference, window_size = 4))
  expect_identical(actual$start, c(1, 5, 9))
  expect_identical(actual$end, c(4, 8, 12))
  expect_identical(actual$gc_skew, c(1, -1, 0))
  expect_identical(actual$gc_ratio, c(1, 0.75, 0))
  expect_equal(actual$gc_skew_minus_average, c(1, -1, 0))
  expect_equal(actual$gc_ratio_minus_average, c(5 / 12, 1 / 6, -7 / 12))
  expect_identical(actual$gc_skew_color, c("#00FF00", "#FF0000", "#00FF00"))
  expect_identical(actual$gc_ratio_color, c("#0000FF", "#0000FF", "#FFFF00"))
  custom <- c(skew_negative = "orange", ratio_positive = "cyan")
  recolored <- suppressMessages(calculate_gc_content(reference, 4, gc_colors = custom))
  values <- setdiff(names(actual), c("gc_skew_color", "gc_ratio_color"))
  expect_identical(recolored[values], actual[values])
  expect_identical(recolored$gc_skew_color, c("#00FF00", "orange", "#00FF00"))
  expect_identical(recolored$gc_ratio_color, c("cyan", "cyan", "#FFFF00"))
})

capture_gc_bars <- function(data, colors = NULL) {
  before <- grDevices::dev.list()
  scratch <- tempfile(fileext = ".pdf")
  grDevices::pdf(scratch)
  output <- tempfile(fileext = ".pdf")
  on.exit({
    circlize::circos.clear()
    for (device in setdiff(grDevices::dev.list(), before)) grDevices::dev.off(device)
    unlink(c(output, scratch))
  }, add = TRUE)
  bars <- list()
  local_mocked_bindings(circos.barplot = function(value, pos, col, bar_width, ...) {
    bars[[length(bars) + 1L]] <<- list(value = value, pos = pos, col = col, width = bar_width)
    invisible(NULL)
  }, .package = "circlize")
  mutations <- data.frame(position = 1000, mutation = "A\u2192G", annotation = "coding",
    gene = "gene", new = TRUE, label = "gene", ALE_version = "S1", version_color = "gray")
  suppressWarnings(suppressMessages(create_full_circos_plot(
    mutation_list = list(S1 = mutations), samples = "S1",
    genome_info = data.frame(seq_id = "chr", length = 100000),
    gc_data = data, output_file = output, gc_colors = colors
  )))
  bars
}

test_that("rendered GC bars use each metric and ignore stale cached colors", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("RColorBrewer")
  data <- data.frame(start = c(1, 10001, 20001), end = c(10000, 20000, 25886),
    gc_skew_minus_average = c(-0.2, 0, 0.3), gc_ratio_minus_average = c(0.04, -0.02, 0),
    gc_skew_color = "#123456", gc_ratio_color = "#654321")
  original <- data
  defaults <- capture_gc_bars(data)
  expect_length(defaults, 2L)
  expect_identical(defaults[[1]]$value, data$gc_skew_minus_average)
  expect_identical(defaults[[2]]$value, data$gc_ratio_minus_average)
  expect_identical(defaults[[1]]$col, c("#FF0000", "#00FF00", "#00FF00"))
  expect_identical(defaults[[2]]$col, c("#0000FF", "#FFFF00", "#0000FF"))
  expect_identical(defaults[[1]]$pos, (data$start + data$end) / 2)
  expect_identical(defaults[[2]]$pos, defaults[[1]]$pos)
  expect_identical(defaults[[1]]$width, c(10000, 10000, 5886))
  expect_identical(defaults[[2]]$width, defaults[[1]]$width)
  custom <- c(skew_positive = "darkgreen", ratio_negative = "orange")
  rendered <- capture_gc_bars(data, custom)
  expect_identical(rendered[[1]]$col, c("#FF0000", "darkgreen", "darkgreen"))
  expect_identical(rendered[[2]]$col, c("#0000FF", "orange", "#0000FF"))
  expect_identical(lapply(rendered, `[`, c("value", "pos")),
                   lapply(defaults, `[`, c("value", "pos")))
  # Legacy metric names remain supported, without requiring cached colors.
  legacy <- data[, 1:4]
  names(legacy)[3:4] <- c("gc_skew_deviation", "gc_ratio_deviation")
  expect_identical(capture_gc_bars(legacy), defaults)
  expect_identical(data, original)
})

test_that("full workflow forwards the palette to calculation and rendering", {
  for (fun in list(run_full_analysis, calculate_gc_content)) {
    expect_identical(tail(names(formals(fun)), 1), "gc_colors")
  }
  expect_identical(tail(names(formals(create_full_circos_plot)), 3),
                   c("gc_colors", "observed_tracks", "reverse_tracks"))
  input <- tempfile()
  dir.create(input)
  reference <- file.path(input, "reference.fasta")
  excel <- file.path(input, "existing.xlsx")
  writeLines(c(">chr", "CCCCGGGAAATT"), reference)
  file.create(excel)
  output <- tempfile()
  on.exit(unlink(c(input, output), recursive = TRUE), add = TRUE)
  captured <- new.env(parent = emptyenv())
  gc_fixture <- data.frame(start = 1, end = 12, gc_skew_minus_average = 0,
                          gc_ratio_minus_average = 0)
  local_mocked_bindings(
    parse_reference_info = function(...) data.frame(seq_id = "chr", length = 12),
    load_mutation_excel = function(...) list(Mutations = data.frame(position = 1)),
    process_mutation_data_for_viz = function(...) {
      list(mutation_list = list(S1 = data.frame(position = 1)), missing_data = NULL, samples = "S1")
    },
    calculate_gc_content = function(...) {
      captured$calculation <- list(...)
      gc_fixture
    },
    create_full_circos_plot = function(...) {
      captured$render <- list(...)
      invisible(NULL)
    }
  )
  colors <- c(skew_positive = "darkgreen", ratio_negative = "orange")
  suppressMessages(run_full_analysis(input_dir = input, output_dir = output,
    reference_file = reference, excel_file = excel, add_genome_map = TRUE,
    create_genbank = FALSE, gc_colors = colors))
  expect_identical(captured$calculation$gc_colors, .breseq_gc_colors(colors))
  expect_identical(captured$render$gc_colors, .breseq_gc_colors(colors))
  expect_identical(captured$render$gc_data, gc_fixture)
  expect_identical(captured$calculation$window_size, 10000)
})

test_that("flat GC tracks have a valid axis without changing zero values", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("RColorBrewer")
  flat <- data.frame(start = 1, end = 10000,
    gc_skew_minus_average = 0, gc_ratio_minus_average = 0)
  original <- flat
  rendered <- capture_gc_bars(flat)
  expect_length(rendered, 2L)
  expect_identical(rendered[[1]]$value, 0)
  expect_identical(rendered[[2]]$value, 0)
  expect_identical(rendered[[1]]$col, "#00FF00")
  expect_identical(rendered[[2]]$col, "#0000FF")
  expect_identical(flat, original)
})
