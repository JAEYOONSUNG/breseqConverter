snapgene_export_fixture <- function(directory) {
  dir.create(directory, recursive = TRUE)
  path <- file.path(directory, "edited reference.gbk")
  writeLines(c("LOCUS       TEST 4 bp DNA circular SYN 01-JAN-2026",
               "ORIGIN", "        1 acgt", "//"), path)
  path
}

snapgene_export_test_bytes <- function() {
  c(as.raw(c(9, 0, 0, 0, 14)), charToRaw("SnapGene"),
    as.raw(c(0, 1, 0, 15, 0, 20, 0, 0, 0, 0, 5, 1)), charToRaw("ACGT"))
}

test_that("SnapGene conversion returns an explicit unavailable result", {
  directory <- tempfile("snapgene-unavailable-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- snapgene_export_fixture(directory)
  output <- file.path(directory, "edited.dna")
  result <- .breseq_convert_genbank_to_dna(input, output,
                                          file.path(directory, "not-installed"))
  expect_false(result$success)
  expect_identical(result$status, "unavailable")
  expect_null(result$output_file)
  expect_false(file.exists(output))
  expect_error(.breseq_convert_genbank_to_dna(input, file.path(directory, "bad.gbk")), ".dna")
  writeLines(rep(readLines(input), 2L), input)
  expect_error(.breseq_convert_genbank_to_dna(input, output), "one complete")
})

test_that("SnapGene conversion leaves an open application and existing outputs alone", {
  directory <- tempfile("snapgene-session-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- snapgene_export_fixture(directory)
  output <- file.path(directory, "edited.dna")
  local_mocked_bindings(find_snapgene = function(...) file.path(R.home("bin"), "R"),
    .breseq_snapgene_is_running = function() TRUE,
    .breseq_snapgene_run = function(...) stop("Must not invoke the application"))
  before <- tools::md5sum(input)
  result <- .breseq_convert_genbank_to_dna(input, output)
  expect_identical(result$status, "application_running")
  expect_false(result$success)
  expect_false(file.exists(output))
  expect_identical(tools::md5sum(input), before)
  writeBin(charToRaw("existing output"), output)
  output_before <- tools::md5sum(output)
  expect_identical(.breseq_convert_genbank_to_dna(input, output)$status, "output_exists")
  expect_identical(tools::md5sum(output), output_before)
})

test_that("SnapGene conversion uses separate arguments and verifies genuine DNA headers", {
  directory <- tempfile("snapgene path with spaces-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- snapgene_export_fixture(directory)
  output <- file.path(directory, "edited sample's output.dna")
  seen <- NULL
  local_mocked_bindings(find_snapgene = function(...) file.path(R.home("bin"), "R"),
    .breseq_snapgene_is_running = function() FALSE,
    .breseq_snapgene_run = function(executable, args, timeout) {
      seen <<- args
      writeBin(snapgene_export_test_bytes(), args[match("--output", args) + 1L])
      list(unavailable = FALSE, pending = FALSE, status = 0L, message = "")
    })
  before <- tools::md5sum(input)
  result <- .breseq_convert_genbank_to_dna(input, output)
  expect_true(result$success)
  expect_identical(result$status, "converted")
  expect_identical(result$output_file, normalizePath(output))
  expect_identical(result$exit_status, 0L)
  expect_identical(seen[match("--input", seen) + 1L], normalizePath(input))
  expect_identical(seen[match("--convert", seen) + 1L], "SnapGene")
  expect_false(any(grepl("quit|kill|register", seen)))
  expect_true(.breseq_snapgene_dna_header(output))
  expect_identical(tools::md5sum(input), before)
  expect_false(any(grepl("snapgene-conversion", list.files(directory, all.files = TRUE))))

  wrong_type <- snapgene_export_test_bytes()
  wrong_type[15L] <- as.raw(2L)
  writeBin(wrong_type, output)
  expect_false(.breseq_snapgene_dna_header(output))
})

test_that("invalid and pending CLI results are never reported as converted", {
  directory <- tempfile("snapgene-failure-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  input <- snapgene_export_fixture(directory)
  output <- file.path(directory, "edited.dna")
  mode <- "invalid"
  local_mocked_bindings(find_snapgene = function(...) file.path(R.home("bin"), "R"),
    .breseq_snapgene_is_running = function() FALSE,
    .breseq_snapgene_run = function(executable, args, timeout) {
      writeLines("not a binary DNA file", args[match("--output", args) + 1L])
      list(unavailable = FALSE, pending = mode == "pending",
           status = if (mode == "failed") 2L else 0L, message = "test outcome")
    })
  invalid <- .breseq_convert_genbank_to_dna(input, output)
  expect_identical(invalid$status, "invalid_output")
  expect_false(invalid$success)
  expect_false(file.exists(output))
  mode <- "failed"
  failed <- .breseq_convert_genbank_to_dna(input, output)
  expect_identical(failed$status, "conversion_failed")
  expect_identical(failed$exit_status, 2L)
  expect_false(file.exists(output))
  mode <- "pending"
  pending <- .breseq_convert_genbank_to_dna(input, output)
  expect_identical(pending$status, "pending")
  expect_false(pending$success)
  expect_true(file.exists(pending$pending_file))
  expect_false(file.exists(output))
})
