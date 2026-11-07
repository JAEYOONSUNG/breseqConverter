publication_fixture <- function() {
  mutation <- function(position, seqid, allele) data.frame(position = position, `seq id` = seqid,
    mutation = allele, evidence = "RA", gene = "gene", description = "product", annotation = "coding", check.names = FALSE)
  mc <- function(start, end, seqid, gene, description) data.frame(start = start, end = end,
    `seq id` = seqid, gene = gene, description = description, check.names = FALSE)
  junction <- function(position, seqid) data.frame(position = position, `seq id` = seqid,
    `reads (cov)` = "10 (1)", gene = "junction gene", check.names = FALSE)
  ma <- mutation(c(2, 30, 1), c("chr", "chr", "p1"), c("C\u2192T", "A\u2192G", "\u0394100 bp"))
  mb <- ma[1, , drop = FALSE]; mb$mutation <- "C\u2192G"; mb$annotation <- "changed coding"
  ca <- mc("10\u201312", "20\u201322", "chr", "L1-L2", "[L1],geneB")
  cb <- mc(c("11\u201313", "10"), c("21\u201323", "90"), c("chr", "p1"),
           c("L1-L2", "P1"), c("[L1],geneB", "hypothetical protein"))
  ja <- junction(c("40 =", "= 50"), c("chr", "chr"))
  jb <- junction(c("40 =", "= 50", "40 =", "= 20"), c("chr", "chr", "chr", "p1"))
  tables <- list(WT = list(`Mutation predictions` = ma[0, ], `Unassigned missing coverage evidence` = ca[0, ],
                          `Unassigned new junction evidence` = ja[0, ]),
                 A = list(`Mutation predictions` = ma, `Unassigned missing coverage evidence` = ca,
                          `Unassigned new junction evidence` = ja),
                 B = list(`Mutation predictions` = mb, `Unassigned missing coverage evidence` = cb,
                          `Unassigned new junction evidence` = jb))
  suppressMessages(merge_comparison_results(compare_multiple_tables(tables)))
}

test_that("compact tables preserve observed calls, allele changes and paired events", {
  merged <- publication_fixture(); original <- merged
  result <- .breseq_publication_tables(merged, c("WT", "A", "B"))
  expect_identical(merged, original)
  expect_identical(names(result$sheets), c("Mutations", "Missing coverage", "Junction evidence"))
  mutation <- result$sheets$Mutations
  expect_equal(dim(mutation$data), c(3, 8))
  expect_identical(mutation$data$WT, rep("Ref.", 3))
  expect_identical(mutation$data$B, c("C\u2192G", "\u2014", "\u2014"))
  expect_identical(mutation$states[, "B"], c("changed", "not_reported", "not_reported"))
  expect_match(mutation$data$`Variant annotation`[1], "A: coding")
  expect_match(mutation$data$`Variant annotation`[1], "B: changed coding")
  missing <- result$sheets[["Missing coverage"]]
  expect_equal(nrow(missing$data), 1L)
  expect_identical(unname(missing$states[1, ]), c("not_reported", "new", "updated"))
  expect_length(result$audit$omitted_mc_families, 1L)
  expect_equal(nrow(result$audit$mc_records), 3L)
  expect_match(missing$data$B, "Start: 11\u201313\\nEnd: 21\u201323")
  junction <- result$sheets[["Junction evidence"]]
  expect_equal(nrow(junction$data), 2L)
  expect_equal(sum(junction$states[, "B"] == "retained"), 1L)
  expect_equal(sum(junction$states[, "B"] == "new"), 1L)
  expect_identical(result$audit$junction_mode, "paired_events")
})

test_that("missing paired keys retain every endpoint instead of inventing junctions", {
  merged <- publication_fixture()
  attr(merged[[3]], "breseq_keys") <- NULL
  result <- .breseq_publication_tables(merged)
  expect_equal(nrow(result$sheets[["Junction evidence"]]$data), 4L)
  expect_identical(result$audit$junction_mode, "unpaired_endpoints")
  expect_match(result$sheets[["Junction evidence"]]$notes[1], "no junction pairs were inferred")
  # Cached annotation alone must not invent a reported endpoint in an empty
  # baseline or override a missing sample-specific seq-id/coverage field.
  merged[[3]][["gene.WT/output/index"]] <- "cached annotation"
  baseline <- .breseq_publication_tables(merged)
  expect_true(all(baseline$sheets[["Junction evidence"]]$states[, "WT"] == "not_reported"))
  expect_false(any(baseline$sheets[["Junction evidence"]]$data$WT == "Detected"))
  bad <- merged
  attr(bad[[3]], "breseq_keys") <- rep("broken", nrow(bad[[3]]))
  expect_error(.breseq_publication_tables(bad), "paired comparison keys")
})

test_that("FASTA plus GFF resolves partial gene tokens without assuming gene loss", {
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fasta <- file.path(directory, "reference.fasta"); gff <- file.path(directory, "reference.gff")
  writeLines(c(">chr", strrep("A", 100), ">p1", strrep("C", 100)), fasta)
  writeLines(c("##gff-version 3",
    "chr\t.\tCDS\t10\t15\t.\t+\t.\tlocus_tag=L1;product=protein%20one",
    "chr\t.\tCDS\t16\t22\t.\t+\t.\tlocus_tag=L2;gene=geneB;product=protein%20two"), gff)
  result <- .breseq_publication_tables(publication_fixture(), reference_file = fasta, gff_file = gff)
  expect_identical(result$sheets[["Missing coverage"]]$data$`Gene / region annotation`, "protein one; protein two")
  expect_identical(result$audit$mc_annotation_map$locus, c("L1", "L2"))
  expect_identical(result$audit$mc_annotation_map$partial, c(TRUE, FALSE))
  expect_match(result$sheets[["Missing coverage"]]$notes[3], "do not establish complete gene loss")
  expect_identical(result$sheets$Mutations$data$Replicon, c("chr", "chr", "p1"))
})

test_that("annotation overrides cannot alter calls and derived evidence stays separate", {
  merged <- publication_fixture()
  override <- data.frame(sheet = "Mutations", seq_id = "chr", position = 2,
                         field = "Gene(s)", value = "gene (user alias)")
  display <- data.frame(sample_id = "B", seq_id = "p1", start = 1, end = 100, size = 100,
    evidence_origin = "COVERAGE_RATIO", source_start = "1", source_end = "100", criterion = "<0.1",
    depth_ratio = 0.001, raw_final_DEL_call = FALSE, appended = TRUE)
  result <- .breseq_publication_tables(merged, plot_evidence = display, annotation_overrides = override)
  expect_identical(result$sheets$Mutations$data$`Gene(s)`[1], "gene (user alias)")
  expect_identical(result$sheets$Mutations$data$B[3], "\u2014")
  expect_equal(nrow(result$sheets[["Plot-derived evidence"]]$data), 1L)
  expect_identical(result$sheets[["Plot-derived evidence"]]$data$`Raw DEL call`, "No")
  plot_alias <- data.frame(sheet = "Plot-derived evidence", seq_id = "p1", position = 1,
    field = "Figure label", value = "Plasmid loss (coverage)")
  labelled <- .breseq_publication_tables(merged, plot_evidence = display,
    annotation_overrides = plot_alias)
  expect_identical(labelled$sheets[["Plot-derived evidence"]]$data$`Figure label`,
    "Plasmid loss (coverage)")
  plot_alias$field <- "Displayed mutation"
  expect_error(.breseq_publication_tables(merged, plot_evidence = display,
    annotation_overrides = plot_alias), "annotation fields only")
  display$appended <- FALSE
  skipped <- .breseq_publication_tables(merged, plot_evidence = display)
  expect_false("Plot-derived evidence" %in% names(skipped$sheets))
  expect_equal(result$audit$annotation_overrides$output_row, 1L)
  override$field <- "B"
  expect_error(.breseq_publication_tables(merged, annotation_overrides = override), "annotation fields only")
  override$field <- "Gene(s)"; override$position <- 99
  expect_error(.breseq_publication_tables(merged, annotation_overrides = override), "exactly one compact row")
  no_wt <- lapply(merged, function(d) d[, !grepl("[.]WT/", names(d)), drop = FALSE])
  result <- .breseq_publication_tables(no_wt)
  expect_identical(result$sample_order, c("A", "B"))
  expect_false("WT" %in% names(result$sheets$Mutations$data))
})

test_that("publication workbook has row-seven headers and usable table relationships", {
  skip_if_not_installed("openxlsx")
  file <- tempfile(fileext = ".xlsx")
  on.exit(unlink(file), add = TRUE)
  prepared <- save_publication_results_to_excel(publication_fixture(), file)
  expect_identical(openxlsx::getSheetNames(file), names(prepared$sheets))
  loaded <- openxlsx::read.xlsx(file, sheet = "Mutations", startRow = 7L,
                                check.names = FALSE, sep.names = " ")
  expect_identical(names(loaded), names(prepared$sheets$Mutations$data))
  expect_identical(loaded$B, prepared$sheets$Mutations$data$B)
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  utils::unzip(file, exdir = directory)
  table <- readLines(list.files(file.path(directory, "xl", "tables"), pattern = "[.]xml$", full.names = TRUE)[1], warn = FALSE)
  expect_match(paste(table, collapse = ""), 'showRowStripes="0"')
  relationships <- list.files(file.path(directory, "xl", "worksheets", "_rels"), full.names = TRUE)
  all_rel <- paste(unlist(lapply(relationships, readLines, warn = FALSE)), collapse = "")
  expect_false(grepl("relationships/(drawing|vmlDrawing)", all_rel))
  # Inspect effective saved cell formats, rather than only the style palette.
  tags <- function(xml, tag) regmatches(xml, gregexpr(paste0("<", tag,
    "(?:\\s[^>]*?)?(?:/>|>.*?</", tag, ">)"), xml, perl = TRUE))[[1]]
  attribute <- function(xml, name) {
    hit <- regmatches(xml, regexec(paste0("\\b", name, '=\"([^\"]*)\"'), xml, perl = TRUE))[[1]]
    if (length(hit)) hit[2] else NA_character_
  }
  xml <- paste(readLines(file.path(directory, "xl", "styles.xml"), warn = FALSE), collapse = "")
  formats <- tags(tags(xml, "cellXfs"), "xf"); fonts <- tags(xml, "font"); fills <- tags(xml, "fill")
  expect_true(all(vapply(fonts, function(font) attribute(tags(font, "name"), "val") == "Arial" &&
                          attribute(tags(font, "sz"), "val") == "10", logical(1))))
  cell_style <- function(sheet, address) {
    cells <- tags(paste(readLines(file.path(directory, "xl", "worksheets", paste0("sheet", sheet, ".xml")), warn = FALSE), collapse = ""), "c")
    cell <- cells[vapply(cells, attribute, "", name = "r") == address]
    format <- formats[as.integer(attribute(cell, "s")) + 1L]
    font <- fonts[as.integer(attribute(format, "fontId")) + 1L]
    fill <- fills[as.integer(attribute(format, "fillId")) + 1L]
    list(bold = grepl("<b(?:/?>|\\s)", font, perl = TRUE), fill = attribute(tags(fill, "fgColor"), "rgb"))
  }
  expect_identical(cell_style(1, "G8"), list(bold = TRUE, fill = "FFDEEBF7"))
  expect_identical(cell_style(1, "H8"), list(bold = TRUE, fill = "FFDAE1EF"))
  expect_identical(cell_style(1, "H9"), list(bold = FALSE, fill = "FFFFFFFF"))
  expect_identical(cell_style(1, "F8"), list(bold = FALSE, fill = "FFF2F7FD"))
  expect_identical(cell_style(2, "F8"), list(bold = FALSE, fill = "FFDAE1EF"))
  expect_identical(cell_style(3, "G8"), list(bold = FALSE, fill = "FFDAE1EF"))
  expect_identical(cell_style(3, "G9"), list(bold = TRUE, fill = "FFDAE1EF"))
})

test_that("real 42 fixture retains the publication counts and stage statuses", {
  fixture <- Sys.getenv("BRESEQCONVERTER_REAL_DATA_RDS", unset = "")
  skip_if(!nzchar(fixture) || !file.exists(fixture), "Optional real-data fixture is not configured")
  merged <- readRDS(fixture)$merged
  result <- .breseq_publication_tables(merged, c("WT", "cc", "1", "3", "6"))
  expect_equal(vapply(result$sheets, function(s) nrow(s$data), 1L),
               c(Mutations = 314L, `Missing coverage` = 6L, `Junction evidence` = 45L))
  expect_equal(sum(result$sheets$Mutations$states != "not_reported"), 883L)
  expect_equal(colSums(result$sheets$Mutations$states == "new"), c(WT = 0L, cc = 88L, `1` = 159L, `3` = 51L, `6` = 16L))
  expect_equal(sum(result$sheets$Mutations$states == "changed"), 1L)
  expect_equal(sum(result$sheets[["Missing coverage"]]$states != "not_reported"), 11L)
  expect_length(result$audit$mc_families, 12L)
  expect_length(result$audit$omitted_mc_families, 6L)
  expect_equal(sum(result$sheets[["Junction evidence"]]$states != "not_reported"), 72L)
})
