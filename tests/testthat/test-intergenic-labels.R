test_that("intergenic filtering changes only the label subset", {
  data <- data.frame(
    annotation = c("intergenic (+1/-2)", " \tINTERGENIC\u00a0(+3/-4)",
                   "intergenic", "intergenic-like protein", "coding (1/9 nt)", NA),
    label = c("one", "two", "three", "four", "intergenic gene name", "six"),
    new = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE),
    mutation = c("A\u2192G", "+A", "\u039430 bp", "C\u2192T", "+C", "G\u2192A"),
    stringsAsFactors = FALSE
  )
  original <- data
  expect_identical(.breseq_filter_intergenic_labels(data), data)
  expect_identical(.breseq_filter_intergenic_labels(data, FALSE), data)
  expect_identical(.breseq_filter_intergenic_labels(data, TRUE), data[4:6, ])
  expect_identical(data, original)
  without_annotation <- data[, setdiff(names(data), "annotation"), drop = FALSE]
  expect_identical(.breseq_filter_intergenic_labels(without_annotation, TRUE), without_annotation)
  expect_identical(.breseq_filter_intergenic_labels(data[FALSE, ], TRUE), data[FALSE, ])
})

intergenic_plot_fixture <- function() {
  data.frame(
    position = seq(1000, 5000, by = 1000),
    mutation = c("A\u2192G", "+A", "\u039430 bp", "C\u2192T", "+C"),
    annotation = c("intergenic (+1/-2)", " INTERGENIC\u00a0(+3/-4)",
                   "coding (1/30 nt)", "pseudogene (2/90 nt)", NA),
    gene = c("flanks1", "flanks2", "intergenic_gene", "pseudo", "unknown"),
    new = TRUE, ALE_version = "S1", version_color = "#66C2A5",
    label = paste0("label", seq_len(5)), stringsAsFactors = FALSE
  )
}

# Capture the actual label and rectangle layers while allowing all plot parsing,
# sample/newness processing, coordinate derivation, and track selection to run.
capture_intergenic_plot <- function(plotter, data, hide) {
  devices_before <- grDevices::dev.list()
  scratch <- tempfile(fileext = ".pdf")
  grDevices::pdf(scratch)
  on.exit({
    devices_after <- grDevices::dev.list()
    for (device in setdiff(devices_after, devices_before)) grDevices::dev.off(device)
    unlink(scratch)
  }, add = TRUE)
  suppressPackageStartupMessages(library(circlize))
  captured <- new.env(parent = emptyenv())
  captured$labels <- list()
  captured$marks <- list()
  testthat::local_mocked_bindings(
    circos.genomicLabels = function(bed, ...) {
      captured$labels <- append(captured$labels, list(bed))
      invisible(NULL)
    },
    circos.rect = function(...) {
      captured$marks <- append(captured$marks, list(list(...)))
      invisible(NULL)
    },
    .package = "circlize"
  )
  output <- tempfile(fileext = ".pdf")
  on.exit(unlink(output), add = TRUE)
  if (identical(plotter, "full")) {
    suppressWarnings(suppressMessages(create_full_circos_plot(
      mutation_list = list(S1 = data), samples = "S1",
      genome_info = data.frame(seq_id = "chr", length = 100000),
      output_file = output, hide_intergenic_labels = hide
    )))
  } else {
    suppressWarnings(suppressMessages(plot_circos_genome(
      mutation_data = data, samples = "S1", seqid = "chr",
      genome_length = 100000, output_file = output,
      hide_intergenic_labels = hide
    )))
  }
  list(labels = captured$labels, marks = captured$marks)
}

test_that("both Circos renderers hide labels without dropping mutation marks", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  data <- intergenic_plot_fixture()
  original <- data
  for (plotter in c("full", "simple")) {
    shown <- capture_intergenic_plot(plotter, data, FALSE)
    hidden <- capture_intergenic_plot(plotter, data, TRUE)
    expect_length(shown$labels, 1L)
    expect_length(hidden$labels, 1L)
    expect_equal(nrow(shown$labels[[1]]), 5L)
    expect_equal(nrow(hidden$labels[[1]]), 3L)
    expect_identical(hidden$labels[[1]]$start, shown$labels[[1]]$start[3:5])
    expect_gt(length(shown$marks), 0L)
    expect_identical(hidden$marks, shown$marks)
    expect_identical(data, original)
  }
})

test_that("the full workflow forwards presentation options through both steps", {
  input <- tempfile()
  dir.create(input)
  writeLines("<html></html>", file.path(input, "index_sample.html"))
  output <- tempfile()
  on.exit(unlink(c(input, output), recursive = TRUE), add = TRUE)
  captured <- new.env(parent = emptyenv())
  data <- intergenic_plot_fixture()
  local_mocked_bindings(
    run_breseq_pipeline = function(...) {
      captured$pipeline <- list(...)
      list(excel_file = "already-exported.xlsx")
    },
    load_mutation_excel = function(...) list(Mutations = data),
    process_mutation_data_for_viz = function(mutation_data, sample_order = NULL) {
      captured$order <- sample_order
      list(mutation_list = list(S1 = data), missing_data = NULL, samples = "S1")
    },
    create_full_circos_plot = function(...) {
      captured$plot <- list(...)
      invisible(NULL)
    }
  )
  labels <- c(S1 = "Published strain")
  suppressMessages(run_full_analysis(
    input_dir = input, output_dir = output, add_genome_map = FALSE,
    create_genbank = FALSE, hide_intergenic_labels = TRUE,
    highlight_new = TRUE, clean_headers = TRUE,
    sample_labels = labels, sample_order = "S1"
  ))
  expect_true(captured$pipeline$hide_intergenic_labels)
  expect_true(captured$pipeline$highlight_new)
  expect_true(captured$pipeline$clean_headers)
  expect_identical(captured$pipeline$sample_labels, labels)
  expect_identical(captured$pipeline$sample_order, "S1")
  expect_identical(captured$order, "S1")
  expect_true(captured$plot$hide_intergenic_labels)
  expect_identical(captured$plot$mutation_list$S1, data)
})

stage_plot_fixture <- function(stages = c("45SJY3", "45SJY3w", "passage.3")) {
  data <- data.frame(position = 100, check.names = FALSE)
  for (i in seq_along(stages)) {
    suffix <- paste0(".", stages[i], "/output/index")
    data[[paste0("mutation", suffix)]] <- "A\u2192G"
    data[[paste0("annotation", suffix)]] <- paste("coding", stages[i])
    data[[paste0("gene", suffix)]] <- paste0("gene", i)
  }
  data
}

test_that("exact stage IDs survive wide plot and workflow parsing", {
  skip_if_not_installed("reshape2")
  skip_if_not_installed("gtools")
  stages <- c("45SJY3", "45SJY3w", "passage.3")
  data <- stage_plot_fixture(stages)
  processed <- suppressMessages(process_mutation_data_for_viz(list(Mutations = data), stages))
  expect_identical(processed$samples, stages)
  expect_identical(names(processed$mutation_list), stages)
  expect_true(processed$mutation_list[[stages[1]]]$new)
  expect_false(processed$mutation_list[[stages[2]]]$new)
  expect_identical(processed$mutation_list[[stages[2]]]$annotation, "coding 45SJY3w")
  attr(data, "breseq_sample_order") <- c("WT", stages[2:1], stages[3])
  metadata <- suppressMessages(process_mutation_data_for_viz(list(Mutations = data)))
  expect_identical(metadata$samples, attr(data, "breseq_sample_order"))
  expect_true(metadata$mutation_list[[stages[2]]]$new)
  expect_false(metadata$mutation_list[[stages[1]]]$new)
  colors <- stats::setNames(rep("black", 3), stages)
  parsed <- suppressMessages(parse_and_prepare_mutations(
    data, "chr", 10000, stages, c(SNP = "black"), colors
  ))
  expect_identical(parsed$ALE_version, stages)
  expect_identical(parsed$annotation, paste("coding", stages))
})

test_that("workflow reshaping keeps replicons separate at a shared position", {
  skip_if_not_installed("reshape2")
  skip_if_not_installed("gtools")
  data <- data.frame(position = c(1, 1), seq_id = c("chr", "plasmid"), check.names = FALSE)
  data[["mutation.S1/output/index"]] <- c("A\u2192G", "C\u2192T")
  data[["mutation.S2/output/index"]] <- c("A\u2192G", "A\u2192G")
  original <- data
  processed <- suppressMessages(process_mutation_data_for_viz(list(Mutations = data), c("S1", "S2")))
  expect_equal(nrow(processed$mutation_list$S1), 2L)
  expect_equal(nrow(processed$mutation_list$S2), 2L)
  expect_identical(processed$mutation_list$S1$mutation, c("A\u2192G", "C\u2192T"))
  expect_identical(processed$mutation_list$S2$new, c(FALSE, TRUE))
  expect_identical(processed$mutation_list$S2[["seq id"]], c("chr", "plasmid"))
  expect_false(".breseq_source_row" %in% names(processed$mutation_list$S1))
  expect_identical(data, original)
})

test_that("missing coverage stage IDs do not truncate 3w or dotted names", {
  skip_if_not_installed("reshape2")
  stages <- c("45SJY3", "45SJY3w", "passage.3")
  missing <- data.frame(start = "100", check.names = FALSE)
  for (i in seq_along(stages)) {
    missing[[paste0("end.", stages[i], "/output/index")]] <- as.character(200 + i)
    missing[[paste0("seq id.", stages[i], "/output/index")]] <- "chr"
  }
  parsed <- process_missing_regions(missing, stages)
  expect_setequal(parsed$ALE_version, stages)
  expect_equal(nrow(parsed), 3L)
})

test_that("Excel loading preserves canonical names before metadata restoration", {
  skip_if_not_installed("openxlsx")
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  data <- stage_plot_fixture()
  data[["seq id.S1/output/index"]] <- "chr"
  data[["reads (cov).1.S1/output/index"]] <- "10 (20)"
  openxlsx::write.xlsx(list(Mutations = data), file)
  captured <- new.env(parent = emptyenv())
  local_mocked_bindings(restore_canonical_export_headers = function(file, sheet, data) {
    captured$names <- names(data)
    captured$sheet <- sheet
    data
  })
  loaded <- suppressMessages(load_mutation_excel(file))
  expect_identical(captured$names, names(data))
  expect_identical(captured$sheet, "Mutations")
  expect_identical(names(loaded$Mutations), names(data))
})

test_that("the pipeline forwards plotting options and separates replicons", {
  skip_if_not_installed("openxlsx")
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("ggplot2")
  input <- tempfile()
  output <- tempfile()
  dir.create(input)
  for (sample in c("S1", "S2")) {
    writeLines("<html></html>", file.path(input, paste0("index_", sample, ".html")))
  }
  on.exit(unlink(c(input, output), recursive = TRUE), add = TRUE)
  raw <- data.frame(
    position = c("100", "100"), `seq id` = c("chr.1", "plasmid.1"),
    mutation = c("A\u2192G", "\u0394100 bp"), evidence = "RA",
    annotation = c("intergenic (+1/-2)", "coding (1/300 nt)"),
    gene = c("flanks", "plasmid_gene"), description = "product",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  captured <- new.env(parent = emptyenv())
  captured$plasmid_only <- FALSE
  captured$plots <- list()
  local_mocked_bindings(
    extract_single_html = function(html_file) {
      data <- if (captured$plasmid_only) raw[2, , drop = FALSE] else raw
      list(`Mutation predictions` = data)
    },
    plot_circos_genome = function(...) {
      captured$plots <- append(captured$plots, list(list(...)))
      invisible(NULL)
    },
    plot_mutations = function(...) NULL,
    plot_mutation_summary = function(...) NULL
  )
  local_mocked_bindings(ggsave = function(...) invisible(NULL), .package = "ggplot2")
  order <- c("S2", "S1")
  labels <- c(S1 = "Early clone", S2 = "Late clone")
  cases <- list(multiple = c(chr.1 = 10000, plasmid.1 = 2000),
                single = c(plasmid.1 = 2000))
  for (case in names(cases)) {
    lengths <- cases[[case]]
    captured$plasmid_only <- identical(case, "single")
    captured$plots <- list()
    result <- NULL
    expect_no_warning(invisible(capture.output(result <- suppressMessages(run_breseq_pipeline(
      input_dir = input, output_dir = file.path(output, case),
      genome_lengths = lengths, create_plots = TRUE, export_snapgene = FALSE,
      hide_intergenic_labels = TRUE, highlight_new = TRUE, clean_headers = TRUE,
      sample_labels = labels, sample_order = order
    )))))
    expect_identical(result$sample_order, order)
    expect_true(file.exists(result$excel_file))
    expect_length(captured$plots, length(lengths))
    expect_setequal(vapply(captured$plots, `[[`, "", "seqid"), names(lengths))
    for (call in captured$plots) {
      expect_true(call$hide_intergenic_labels)
      expect_identical(call$samples, order)
      expect_identical(call$genome_length, unname(lengths[[call$seqid]]))
      expect_equal(nrow(call$mutation_data), 1L)
      expect_identical(.breseq_wide_seq(call$mutation_data,
        .breseq_stage_columns(call$mutation_data)), .breseq_seq(call$seqid))
      expect_identical(as.character(call$mutation_data$position), "100")
    }
  }
})
