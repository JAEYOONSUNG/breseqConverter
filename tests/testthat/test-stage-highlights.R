mutation_fixture <- function(position, sequence, mutation, evidence = "RA", gene = "gene") {
  data.frame(position = position, `seq id` = sequence, mutation = mutation,
             evidence = rep(evidence, length(position)), annotation = rep("coding", length(position)),
             gene = rep(gene, length(position)), description = rep("product", length(position)),
             check.names = FALSE, stringsAsFactors = FALSE)
}

# Inspect the saved OOXML, including the style actually assigned to each cell.
# This uses base R so the regression check adds no test-only dependency.
saved_excel_styles <- function(file) {
  directory <- tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  utils::unzip(file, exdir = directory)
  read_xml <- function(path) paste(readLines(path, warn = FALSE), collapse = "")
  tags <- function(xml, tag) {
    pattern <- paste0("<", tag, "(?:\\s[^>]*?)?(?:/>|>.*?</", tag, ">)")
    regmatches(xml, gregexpr(pattern, xml, perl = TRUE))[[1]]
  }
  attribute <- function(xml, name) {
    if (!length(xml)) return(NA_character_)
    hit <- regmatches(xml, regexec(paste0("\\b", name, "=\"([^\"]*)\""), xml, perl = TRUE))[[1]]
    if (length(hit)) hit[2] else NA_character_
  }
  styles <- read_xml(file.path(directory, "xl", "styles.xml"))
  fonts <- tags(styles, "font")
  fills <- tags(styles, "fill")
  formats <- tags(tags(styles, "cellXfs"), "xf")
  sheets <- list.files(file.path(directory, "xl", "worksheets"), "^sheet[0-9]+[.]xml$", full.names = TRUE)
  lapply(sheets, function(sheet) {
    cells <- tags(read_xml(sheet), "c")
    rows <- lapply(cells, function(cell) {
      style_id <- attribute(cell, "s")
      format <- formats[if (is.na(style_id)) 1L else as.integer(style_id) + 1L]
      font <- fonts[as.integer(attribute(format, "fontId")) + 1L]
      fill <- fills[as.integer(attribute(format, "fillId")) + 1L]
      data.frame(cell = attribute(cell, "r"),
        font = attribute(tags(font, "name"), "val"),
        size = attribute(tags(font, "sz"), "val"),
        bold = grepl("<b(?:/?>|\\s)", font, perl = TRUE),
        fill = attribute(tags(fill, "fgColor"), "rgb"), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  })
}

test_that("stage parsing preserves exact IDs and distinguishes dotted field names", {
  d <- as.data.frame(setNames(rep(list(character()), 7), c(
    "position", "reads (cov).1", "mutation.45SJY3w/output/index",
    "annotation.45SJY3w/output/index", "mutation.strain.2/output/index",
    "reads (cov).1.strain.2/output/index", "reads (cov).1/output/index")), check.names = FALSE)
  m <- .breseq_stage_columns(d)
  expect_false(2L %in% m$index)
  expect_identical(m$sample, c("45SJY3w", "45SJY3w", "strain.2", "strain.2", "1"))
  expect_identical(m$field, c("mutation", "annotation", "mutation", "reads (cov).1", "reads (cov)"))
  expect_error(.breseq_stage_columns(d, c("45SJY3w", "strain.2")), "sample_order")
})

test_that("replicon-aware comparisons preserve alleles and classify transitions", {
  a <- mutation_fixture(c(1L, 1L), c("chr.1", "plasmid.1"), c("A\u2192G", "\u0394100 bp"))
  b <- mutation_fixture(c(1L, 1L), c("chr", "plasmid"), c("A\u2192T", "\u0394100 bp"), "JC", "renamed")
  c <- b[2, , drop = FALSE]
  samples <- list(WT = list(`Mutation predictions` = a[0, ]),
                  A = list(`Mutation predictions` = a), B = list(`Mutation predictions` = b),
                  C = list(`Mutation predictions` = c), D = list(`Mutation predictions` = b))
  merged <- merge_comparison_results(compare_multiple_tables(samples))[[1]]
  expect_equal(nrow(merged), 2L)
  expect_length(unique(attr(merged, "breseq_keys")), 2L)
  status <- .breseq_stage_status(merged, "Mutations", names(samples))
  seqid <- .breseq_wide_seq(merged, .breseq_stage_columns(merged))
  expect_identical(unname(status$states[seqid == "chr", ]),
                   c("not_reported", "new", "changed", "not_reported", "reappeared"))
  expect_identical(unname(status$states[seqid == "plasmid", ]),
                   c("not_reported", "new", "retained", "retained", "retained"))
  expect_identical(a$mutation, c("A\u2192G", "\u0394100 bp"))
  expect_error(compare_multiple_tables(samples, c("WT", "A")), "every input sample")
})

test_that("junction identity includes both oriented endpoints and not read metrics", {
  junction <- function(position, reads) data.frame(position = position,
    `seq id` = rep("chr", length(position)), `reads (cov)` = reads,
    gene = rep("gene", length(position)), check.names = FALSE)
  a <- junction(c("100 =", "= 200"), c("10 (1)", "20 (2)"))
  # Reverse endpoint row order for the retained event; also reuse endpoint100
  # for a different partner. The latter must remain a distinct new event.
  b <- junction(c("= 200", "100 =", "100 =", "= 300"), rep("99 (9)", 4))
  raw <- list(WT = list(`Unassigned new junction evidence` = a[0, ]),
              A = list(`Unassigned new junction evidence` = a),
              B = list(`Unassigned new junction evidence` = b))
  d <- merge_comparison_results(compare_multiple_tables(raw))[[1]]
  expect_equal(nrow(d), 4L)
  st <- .breseq_stage_status(d, "New junction", names(raw))$states
  expect_equal(sum(st[, "A"] == "new"), 2L)
  expect_equal(sum(st[, "B"] == "retained"), 2L)
  expect_equal(sum(st[, "B"] == "new"), 2L)
  expect_error(.breseq_stage_status(structure(d, breseq_keys = NULL), "New junction"), "paired comparison keys")
})

test_that("MC boundary revisions and evidence inside a prior DEL are not new losses", {
  mc <- function(start, end, seqid) data.frame(start = start, end = end,
    `seq id` = seqid, gene = rep("g1-g2", length(start)), check.names = FALSE)
  a <- mc("10\u201312", "20\u201322", "chr")
  b <- mc(c("11\u201313", "10"), c("21\u201323", "90"), c("chr", "plasmid"))
  ma <- mutation_fixture(1L, "plasmid", "\u0394100 bp")
  mb <- ma[0, ]
  raw <- list(A = list(`Mutation predictions` = ma, `Unassigned missing coverage evidence` = a),
              B = list(`Mutation predictions` = mb, `Unassigned missing coverage evidence` = b))
  merged <- merge_comparison_results(compare_multiple_tables(raw))
  st <- .breseq_stage_status(merged[[2]], "Missing", c("A", "B"), merged[[1]])$states
  expect_equal(sum(st == "new"), 1L)
  expect_equal(sum(st == "updated"), 1L)
  expect_equal(sum(st == "prior_del_region"), 1L)
})

test_that("clean workbook headers round trip and defaults retain canonical data", {
  skip_if_not_installed("openxlsx")
  a <- mutation_fixture(c(10L, 20L), rep("chr", 2), c("A\u2192G", "+T"))
  b <- a; b$mutation[1] <- "A\u2192T"
  samples <- list(`45SJY3` = list(`Mutation predictions` = a),
                  `45SJY3w` = list(`Mutation predictions` = b))
  merged <- merge_comparison_results(compare_multiple_tables(samples))
  unchanged <- merged
  normal <- tempfile(fileext = ".xlsx"); pretty <- tempfile(fileext = ".xlsx")
  plain <- tempfile(fileext = ".xlsx")
  on.exit(unlink(c(normal, pretty, plain)), add = TRUE)
  save_merged_results_to_excel(merged, normal)
  expect_identical(openxlsx::getSheetNames(normal), "Mutations")
  normal_data <- openxlsx::read.xlsx(normal, check.names = FALSE, sep.names = " ")
  expect_identical(names(normal_data), names(merged[[1]]))
  expect_equal(unname(as.matrix(normal_data)), unname(as.matrix(merged[[1]])), ignore_attr = TRUE)
  baseline <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(baseline, "Mutations")
  openxlsx::writeData(baseline, "Mutations", merged[[1]])
  openxlsx::saveWorkbook(baseline, plain)
  expect_identical(saved_excel_styles(normal), saved_excel_styles(plain))
  labels <- c(`45SJY3` = "Parent strain", `45SJY3w` = "Clone w")
  save_merged_results_to_excel(merged, pretty, highlight_new = TRUE,
    clean_headers = TRUE, sample_labels = labels, sample_order = names(samples))
  expect_true("_breseq_columns" %in% openxlsx::getSheetNames(pretty))
  wb <- openxlsx::loadWorkbook(pretty)
  expect_identical(openxlsx::sheetVisibility(wb)[2], "hidden")
  display <- openxlsx::read.xlsx(pretty, "Mutations", check.names = FALSE, skipEmptyCols = FALSE)
  restored <- restore_canonical_export_headers(pretty, "Mutations", display)
  expect_identical(names(restored), names(merged[[1]]))
  expect_identical(attr(restored, "breseq_sample_order"), names(samples))
  expect_equal(unname(as.matrix(restored)), unname(as.matrix(merged[[1]])), ignore_attr = TRUE)
  expect_error(restore_canonical_export_headers(pretty, "Mutations", display[, rev(seq_len(ncol(display)))]),
               "headers differ")
  expect_identical(merged, unchanged)
})

test_that("all exporter entry points accept opt-in presentation options", {
  skip_if_not_installed("openxlsx")
  a <- mutation_fixture(1L, "chr", "A\u2192G")
  comparisons <- compare_named_tables(list(`Mutation predictions` = a[0, ]),
    list(`Mutation predictions` = a), "WT/output/index", "sample.3w/output/index")
  single <- tempfile(fileext = ".xlsx"); directory <- tempfile()
  on.exit(unlink(c(single, directory), recursive = TRUE), add = TRUE)
  save_comparison_results_to_excel(comparisons, single, highlight_new = TRUE,
    clean_headers = TRUE, sample_labels = c(`sample.3w` = "Clone w"))
  save_results_to_excel(comparisons, directory, highlight_new = TRUE, clean_headers = TRUE)
  expect_true(file.exists(single))
  expect_true(file.exists(file.path(directory, "Mutation predictions.xlsx")))
})

test_that("stage colors show all observed calls while bold distinguishes transitions", {
  skip_if_not_installed("openxlsx")
  a <- mutation_fixture(c(10L, 20L), c("chr", "chr"), c("A\u2192G", "+T"))
  b <- a; b$mutation[1] <- "A\u2192T"
  raw <- list(WT = list(`Mutation predictions` = a[0, ]),
              cc = list(`Mutation predictions` = a),
              `1` = list(`Mutation predictions` = b),
              `3` = list(`Mutation predictions` = b[2, , drop = FALSE]),
              `6` = list(`Mutation predictions` = b))
  merged <- merge_comparison_results(compare_multiple_tables(raw))
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  save_merged_results_to_excel(merged, file, highlight_new = TRUE, clean_headers = TRUE)
  saved <- saved_excel_styles(file)
  expect_true(all(vapply(saved, function(x) all(x$font == "Arial" & x$size == "10"), logical(1))))
  cells <- saved[[1]]
  headers <- cells[grepl("^[A-Z]+1$", cells$cell), ]
  expect_true(all(headers$bold))
  expect_true(all(headers$fill == "FFFFFFFF"))
  columns <- .breseq_stage_columns(merged[[1]])
  status <- .breseq_stage_status(merged[[1]], "Mutations")$states
  colors <- c(WT = "FFF3F3F3", cc = "FFEEF3F8", `1` = "FFF0F5EF", `3` = "FFF7F2EA", `6` = "FFF3F0F6")
  for (sample in names(raw)) for (row in seq_len(nrow(merged[[1]]))) {
    addresses <- paste0(openxlsx::int2col(columns$index[columns$sample == sample]), row + 1L)
    actual <- cells[match(addresses, cells$cell), ]
    state <- status[row, sample]
    if (state == "not_reported") {
      expect_true(all(is.na(actual$fill)), info = paste(sample, row, "absent"))
    } else {
      expect_true(all(actual$fill == colors[[sample]]), info = paste(sample, row, state))
    }
    expect_true(all(actual$bold == (state %in% c("new", "changed", "reappeared"))),
      info = paste(sample, row, state))
  }
  # Clean headers alone opt into Arial 10, without stage-presence fills.
  save_merged_results_to_excel(merged, file, clean_headers = TRUE)
  clean <- saved_excel_styles(file)
  expect_true(all(vapply(clean, function(x) all(x$font == "Arial" & x$size == "10"), logical(1))))
  body <- clean[[1]][!grepl("^[A-Z]+1$", clean[[1]]$cell), ]
  expect_true(all(is.na(body$fill)))
  expect_false(any(body$bold))
})

test_that("observed WT calls are gray and MC revisions are filled but not bold", {
  skip_if_not_installed("openxlsx")
  mc <- function(start, end, seqid) data.frame(start = start, end = end,
    `seq id` = seqid, gene = rep("g1-g2", length(start)), check.names = FALSE)
  wt <- mutation_fixture(1L, "plasmid", "\u0394100 bp")
  raw <- list(WT = list(`Mutation predictions` = wt,
                       `Unassigned missing coverage evidence` = mc("10\u201312", "20\u201322", "chr")),
              cc = list(`Mutation predictions` = wt[0, ],
                        `Unassigned missing coverage evidence` = mc(c("11\u201313", "10"),
                          c("21\u201323", "90"), c("chr", "plasmid"))))
  merged <- merge_comparison_results(compare_multiple_tables(raw))
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  save_merged_results_to_excel(merged, file, highlight_new = TRUE, clean_headers = TRUE)
  saved <- saved_excel_styles(file)
  expect_true(all(vapply(saved, function(x) all(x$font == "Arial" & x$size == "10"), logical(1))))
  # The hidden metadata sheet is inserted between the two visible sheets.
  saved <- saved[match(c("Mutations", "Missing"), openxlsx::getSheetNames(file))]
  mutation_cols <- .breseq_stage_columns(merged[[1]])
  wt_cells <- paste0(openxlsx::int2col(mutation_cols$index[mutation_cols$sample == "WT"]), "2")
  expect_true(all(saved[[1]]$fill[match(wt_cells, saved[[1]]$cell)] == "FFF3F3F3"))
  mc_cols <- .breseq_stage_columns(merged[[2]])
  states <- .breseq_stage_status(merged[[2]], "Missing", mutation_data = merged[[1]])$states
  expect_setequal(states[states != "not_reported"], c("new", "updated", "prior_del_region"))
  for (row in which(states[, "cc"] != "not_reported")) {
    addresses <- paste0(openxlsx::int2col(mc_cols$index[mc_cols$sample == "cc"]), row + 1L)
    actual <- saved[[2]][match(addresses, saved[[2]]$cell), ]
    expect_true(all(actual$fill == "FFEEF3F8"))
    expect_false(any(actual$bold))
  }
})
