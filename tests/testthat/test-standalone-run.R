test_that("one-file driver runs without pre-existing genome objects", {
  skip_if_not_installed("circlize")
  skip_if_not_installed("ComplexHeatmap")
  input <- tempfile("standalone-input-")
  output <- tempfile("standalone-output-")
  dir.create(input)
  on.exit(unlink(c(input, output), recursive = TRUE), add = TRUE)
  writeLines(c(">chr circular", paste(rep("ACGT", 5000), collapse = ""),
               ">plasmid circular", paste(rep("GGCA", 2000), collapse = "")),
             file.path(input, "reference.fasta"))
  mutation <- data.frame(seq_id = c("chr", "plasmid"), position = c(10, 20),
                         check.names = FALSE)
  for (s in c("A", "B")) {
    mutation[[paste0("seq id.", s, "/output/index")]] <- c("chr", "plasmid")
    mutation[[paste0("mutation.", s, "/output/index")]] <- c("C\u2192G", "\u039410 bp")
    mutation[[paste0("annotation.", s, "/output/index")]] <- c("intergenic (+1/-2)", "coding")
  }
  file <- file.path(input, "merged_mutation_results.xlsx")
  openxlsx::write.xlsx(list(Mutations = mutation), file)
  cfg <- list(input_dir = input, output_dir = output, reference_file = NULL,
    excel_file = file, gff_file = NULL, eggnog_file = NULL,
    baseline_id = "WT", sample_order = c("WT", "A", "B"), sample_labels = NULL,
    sample_aliases = NULL, display_deletions_file = NULL,
    hide_intergenic_labels = TRUE, gc_colors = NULL, create_genbank = TRUE, create_dna = FALSE)
  before <- tools::md5sum(c(file, file.path(input, "reference.fasta")))
  devices <- grDevices::dev.list()
  on.exit({
    for (device in setdiff(grDevices::dev.list(), devices)) grDevices::dev.off(device)
  }, add = TRUE)
  result <- suppressMessages(.breseq_run_standalone(cfg))
  expect_setequal(names(result$plots), c("chr", "plasmid"))
  expect_true(all(file.exists(result$plots)))
  expect_true(all(file.info(result$plots)$size > 0))
  expect_true(file.exists(result$excel_file))
  expect_identical(result$counts$raw_calls, c(0L, 2L, 2L))
  expect_identical(result$counts$plot_calls, result$counts$raw_calls)
  expect_equal(nrow(result$genome_export$manifest), 6L)
  expect_true(all(file.exists(result$genome_export$manifest$genbank_file)))
  expect_equal(result$genome_export$manifest$new_length,
               c(20000, 8000, 20000, 7990, 20000, 7990))
  expect_true(all(result$genome_export$manifest$dna_status == "not_requested"))
  edits <- result$genome_export$edits
  expect_equal(sum(edits$ref_length - edits$alt_length), 20)
  expect_identical(before, tools::md5sum(names(before)))
  expect_false(file.exists(file.path(input, "plots")))
  cfg$output_dir <- file.path(input, "forbidden-output")
  expect_error(.breseq_run_standalone(cfg), "outside input_dir")
  expect_false(dir.exists(cfg$output_dir))
})

test_that("unresolved coverage remains evidence even without a plot sidecar", {
  data <- data.frame(start = c("100\u2013110", "500"), check.names = FALSE)
  data[["seq id.A/output/index"]] <- c("chr", "plasmid")
  data[["end.A/output/index"]] <- c("190\u2013200", "600")
  data[["gene.A/output/index"]] <- c("gene1", "gene2")
  original <- data
  unresolved <- .breseq_standalone_unresolved_coverage(data, c("WT", "A"))
  expect_equal(nrow(unresolved), 2L)
  expect_identical(unresolved$sample, c("A", "A"))
  expect_identical(unresolved$source_start, data$start)
  expect_true(all(grepl("not applied", unresolved$interpretation)))
  expect_identical(data, original)
})
