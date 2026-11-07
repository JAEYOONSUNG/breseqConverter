capture_observed_tracks <- function(mutations, samples, observed = NULL, reverse = NULL) {
  before <- grDevices::dev.list()
  scratch <- tempfile(fileext = ".pdf")
  output <- tempfile(fileext = ".pdf")
  grDevices::pdf(scratch)
  on.exit({
    circlize::circos.clear()
    for (device in setdiff(grDevices::dev.list(), before)) grDevices::dev.off(device)
    unlink(c(output, scratch))
  }, add = TRUE)
  marks <- list()
  labels <- list()
  local_mocked_bindings(
    circos.genomicLabels = function(bed, ...) {
      labels[[length(labels) + 1L]] <<- bed
      invisible(NULL)
    },
    circos.rect = function(xleft, ybottom, xright, ytop, col, ...) {
      if (identical(ybottom, 0.1) && identical(ytop, 0.9)) {
        marks[[length(marks) + 1L]] <<- data.frame(
          track = circlize::get.current.track.index(), start = xleft,
          end = xright, color = unname(col))
      }
      invisible(NULL)
    }, .package = "circlize"
  )
  args <- list(mutation_list = mutations, samples = samples,
    genome_info = data.frame(seq_id = "chr", length = 100000), output_file = output)
  if (!is.null(observed)) args$observed_tracks <- observed
  if (!is.null(reverse)) args$reverse_tracks <- reverse
  suppressWarnings(suppressMessages(do.call(create_full_circos_plot, args)))
  list(marks = if (length(marks)) do.call(rbind, marks) else data.frame(), labels = labels)
}

test_that("observed rings retain only the current sample calls and all retained marks", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("RColorBrewer")
  s1 <- data.frame(position = c(1000, 4000), mutation = c("A\u2192G", "\u03941000 bp"),
    annotation = "coding", gene = c("g1", "g2"), new = TRUE,
    label = c("g1", "g2"), version_color = "gray", stringsAsFactors = FALSE)
  s2 <- s1; s2$mutation[1] <- "+A"; s2$new <- c(TRUE, FALSE)
  s3 <- s2[2, , drop = FALSE]
  s4 <- s2; s4$new <- FALSE
  mutations <- list(WT = s1[0, ], S1 = s1, S2 = s2, S3 = s3, S4 = s4)
  original <- mutations
  observed <- capture_observed_tracks(mutations, names(mutations), TRUE)
  legacy <- capture_observed_tracks(mutations, names(mutations))
  explicit_legacy <- capture_observed_tracks(mutations, names(mutations), FALSE)
  expect_identical(legacy, explicit_legacy)
  expect_identical(mutations, original)
  # Tracks run from the last sample outside to WT inside, in both modes.
  expect_identical(as.integer(table(observed$marks$track)), c(2L, 1L, 2L, 2L))
  expect_identical(as.integer(table(legacy$marks$track)), c(7L, 5L, 4L, 2L))
  expect_identical(observed$marks$start, c(1000, 4000, 4000, 1000, 4000, 1000, 4000))
  expect_identical(observed$marks$color,
    c("#FF5C8D", "#03045E", "#03045E", "#FF5C8D", "#03045E", "#139487", "#03045E"))
  # The changed allele has only its current INS mark; the missing SNP is not
  # carried into S3, and retained/reappeared calls with new=FALSE still draw.
  expect_equal(sum(observed$marks$start == 1000 & observed$marks$color == "#139487"), 1L)
  expect_equal(sum(observed$marks$start == 1000 & observed$marks$color == "#FF5C8D"), 2L)
  expect_identical(observed$labels, legacy$labels)
  expect_equal(nrow(observed$labels[[1]]), 3L)
  forward <- capture_observed_tracks(mutations, names(mutations), TRUE, FALSE)
  expect_identical(as.integer(table(forward$marks$track)), c(2L, 2L, 1L, 2L))
  expect_identical(forward$marks$start, c(1000, 4000, 1000, 4000, 4000, 1000, 4000))
  expect_identical(forward$marks$color,
    c("#139487", "#03045E", "#FF5C8D", "#03045E", "#03045E", "#FF5C8D", "#03045E"))
  # The empty WT ring is first, so the first mutation mark is on track 2.
  expect_equal(min(forward$marks$track), 2L)
  expect_identical(forward$labels, observed$labels)
  expect_identical(mutations, original)
})

test_that("empty reference baselines draw without fabricated mutation marks", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("RColorBrewer")
  baseline <- list(WT = data.frame())
  expect_equal(nrow(capture_observed_tracks(baseline, "WT", TRUE)$marks), 0L)
  expect_equal(nrow(capture_observed_tracks(list(), "WT", TRUE)$marks), 0L)
  expect_equal(nrow(capture_observed_tracks(baseline, "WT")$marks), 0L)
  expect_identical(baseline, list(WT = data.frame()))
  expect_error(create_full_circos_plot(list(), samples = "WT", observed_tracks = NA), "observed_tracks")
  expect_error(create_full_circos_plot(list(), samples = "WT", reverse_tracks = NA), "reverse_tracks")
})

test_that("renderer closes only its PDF without creating a default graphics file", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("RColorBrewer")
  directory <- tempfile()
  dir.create(directory)
  previous_dir <- getwd()
  previous_devices <- grDevices::dev.list()
  setwd(directory)
  on.exit({
    setwd(previous_dir)
    for (device in setdiff(grDevices::dev.list(), previous_devices)) grDevices::dev.off(device)
    unlink(directory, recursive = TRUE)
  }, add = TRUE)
  output <- "reference.pdf"
  suppressWarnings(suppressMessages(create_full_circos_plot(
    mutation_list = list(WT = data.frame()), samples = "WT",
    genome_info = data.frame(seq_id = "chr", length = 100000),
    output_file = output, observed_tracks = TRUE, reverse_tracks = FALSE
  )))
  expect_true(file.exists(output))
  expect_identical(list.files(directory), output)
  expect_identical(grDevices::dev.list(), previous_devices)
  # The same owned-device cleanup also applies when rendering fails.
  local_mocked_bindings(circos.initialize = function(...) stop("forced plot error"),
                       .package = "circlize")
  expect_error(suppressMessages(create_full_circos_plot(
    mutation_list = list(WT = data.frame()), samples = "WT",
    genome_info = data.frame(seq_id = "chr", length = 100000),
    output_file = "failed.pdf", observed_tracks = TRUE
  )), "forced plot error")
  expect_false(file.exists("Rplots.pdf"))
  expect_identical(grDevices::dev.list(), previous_devices)
})
