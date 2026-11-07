standalone_evidence_fixture <- function() {
  row <- function(sample, seqid, position, mutation) {
    data <- data.frame(
      position = as.integer(position), `seq id` = seqid, mutation = mutation,
      annotation = "coding", evidence = "RA", gene = "original_gene",
      description = "original product", new = FALSE, ALE_version = sample,
      version_color = "#123456", label = paste("original", sample),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    # Processed source frames can retain nonconsecutive integer source-row names.
    attr(data, "row.names") <- as.integer(position)
    data
  }
  list(
    mutations = list(S1 = row("S1", "chr", 10L, "A\u2192G"),
      S3 = row("S3", "plasmid", 1L, "\u03941,000 bp"),
      S6 = row("S6", "chr", 20L, "C\u2192T")),
    reference = data.frame(seq_id = c("chr", "plasmid"), length = c(10000, 1000)),
    aliases = c(`1` = "S1", `3` = "S3", `6` = "S6"),
    evidence = data.frame(
      sample = c("6", "1"), seq_id = c("plasmid", "chr"),
      start = c(1L, 110L), end = c(1000L, 120L), size = c(1000L, 11L),
      evidence_origin = c("COVERAGE_RATIO", "MC"),
      source_start = c("1", "100\u2013110"), source_end = c("1000", "130\u2013120"),
      depth_ratio = c(0.0019, NA), criterion = c("<0.1", "audited MC conservative core"),
      raw_final_DEL_call = FALSE, stringsAsFactors = FALSE),
    missing = data.frame(start = "100\u2013110", `seq id.S1/output/index` = "chr",
      `end.S1/output/index` = "130\u2013120", `gene.S1/output/index` = "MC_gene",
      check.names = FALSE)
  )
}

with_standalone_evidence <- function(evidence, action) {
  file <- tempfile(fileext = ".tsv")
  on.exit(unlink(file), add = TRUE)
  utils::write.table(evidence, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  action(file)
}

test_that("audited evidence preserves p1 loss in the observed sample and all source rows", {
  fixture <- standalone_evidence_fixture()
  original <- fixture$mutations
  result <- with_standalone_evidence(fixture$evidence, function(file) {
    .breseq_standalone_display_deletions(fixture$mutations, file, fixture$aliases,
      fixture$reference, fixture$missing)
  })
  expect_identical(fixture$mutations, original)
  for (sample in names(original)) {
    expect_identical(result$mutation_list[[sample]][seq_len(nrow(original[[sample]])), ],
                     original[[sample]])
  }
  expect_identical(vapply(result$mutation_list, nrow, integer(1)), c(S1 = 2L, S3 = 1L, S6 = 2L))
  added <- result$mutation_list$S6[2, ]
  expect_identical(added[["seq id"]], "plasmid")
  expect_identical(added$position, 1L)
  expect_identical(added$mutation, "\u03941,000 bp")
  expect_identical(added$annotation, "coverage-derived DEL display")
  expect_identical(added$evidence, "COVERAGE_RATIO")
  expect_false(added$new) # The same full deletion was already present in S3.
  expect_identical(added$ALE_version, "S6")
  expect_identical(result$mutation_list$S1$gene[2], "MC_gene")
  expect_true(result$mutation_list$S1$new[2])
  expect_identical(result$display_rows$sample_id, c("S6", "S1"))
  expect_true(all(result$display_rows$appended))
  expect_identical(result$display_rows$sample, c(6L, 1L))
})

test_that("covered source deletions and empty audited sidecars add no duplicate rows", {
  fixture <- standalone_evidence_fixture()
  evidence <- fixture$evidence[1, ]
  evidence$sample <- "3"
  result <- with_standalone_evidence(evidence, function(file) {
    .breseq_standalone_display_deletions(fixture$mutations, file, fixture$aliases, fixture$reference)
  })
  expect_identical(result$mutation_list, fixture$mutations)
  expect_false(result$display_rows$appended)
  expect_match(result$display_rows$skip_reason, "already covered")
  result <- with_standalone_evidence(evidence[FALSE, ], function(file) {
    .breseq_standalone_display_deletions(fixture$mutations, file, fixture$aliases, fixture$reference)
  })
  expect_identical(result$mutation_list, fixture$mutations)
  expect_equal(nrow(result$display_rows), 0L)
})

test_that("invalid or unattributed display evidence fails before appending rows", {
  fixture <- standalone_evidence_fixture()
  check <- function(data, pattern) with_standalone_evidence(data, function(file) {
    expect_error(.breseq_standalone_display_deletions(fixture$mutations, file,
      fixture$aliases, fixture$reference), pattern)
  })
  data <- fixture$evidence; data$depth_ratio[1] <- 0.1
  check(data, "depth_ratio <0.1")
  data <- fixture$evidence; data$size[1] <- 999
  check(data, "size must equal")
  data <- fixture$evidence; data$sample[1] <- "unknown"
  check(data, "unknown sample")
  data <- fixture$evidence; data$seq_id[1] <- "unknown"
  check(data, "unknown reference")
  data <- fixture$evidence; data$start[2] <- 0
  check(data, "positive")
  data <- fixture$evidence; data$raw_final_DEL_call[1] <- TRUE
  check(data, "raw_final_DEL_call = FALSE")
  data <- fixture$evidence; data$evidence_origin[1] <- "RA"
  check(data, "evidence_origin")
  data <- fixture$evidence; data$start[1] <- 2; data$size[1] <- 999
  check(data, "whole replicon")
})
