sequence_calls <- function(position, mutation, ...) {
  data.frame(position = position, mutation = mutation, ..., check.names = FALSE,
             stringsAsFactors = FALSE)
}

test_that("sequence edits use reference coordinates and report lifted bounds", {
  reference <- "ACGTACGT"
  calls <- sequence_calls(c(7, 2, 4), c("\u03942\u00a0bp", "C\u2192T", "\u200b+AA"))
  original <- calls
  result <- .breseq_apply_sequence_edits(reference, calls)
  expect_identical(result$sequence, "ATGTAAAC")
  expect_identical(calls, original)
  expect_true(result$complete)
  expect_identical(result$edits$type, c("SNP", "INS", "DEL"))
  expect_identical(result$edits$ref_start, c(2, 5, 7))
  expect_identical(result$edits$ref_end, c(2, 4, 8))
  expect_identical(result$edits$new_start, c(2, 5, 9))
  expect_identical(result$edits$new_end, c(2, 6, 8))
  expect_identical(result$edits$ref_sequence, c("C", "", "GT"))
  expect_identical(result$report$edit_index, c(3L, 1L, 2L))
  expect_true(all(result$report$status == "applied"))
  expect_identical(.breseq_apply_sequence_edits("ACGT", sequence_calls(4, "+AA"))$sequence, "ACGTAA")
  expect_identical(.breseq_apply_sequence_edits("ACGT", sequence_calls(1, "1 bp\u2192TG"))$sequence, "TGCGT")
  expect_identical(.breseq_apply_sequence_edits("ACGT", sequence_calls(2, "-CG"))$sequence, "AT")
})

test_that("repeat display positions are the underlying INS or DEL coordinates", {
  reference <- paste0("ACG", strrep("T", 8), "G", strrep("AC", 3), "GCCATGG")
  calls <- sequence_calls(c(20, 11, 2, 17, 25),
    c("2 bp\u2192TA", "(T)8\u21929", "C\u2192A", "(AC)3\u21922", "+GG"))
  result <- .breseq_apply_sequence_edits(reference, calls)
  expect_identical(result$sequence, paste0("AAG", strrep("T", 9), "G", strrep("AC", 2), "GTAATGGGG"))
  expect_identical(result$edits$ref_length, c(1, 0, 2, 2, 0))
  expect_equal(result$edits$alt_length, c(1, 1, 0, 2, 2))
  expect_identical(result$edits$ref_start, c(2, 12, 17, 20, 26))
  expect_identical(result$edits$ref_end, c(2, 11, 18, 21, 25))
  expect_error(.breseq_apply_sequence_edits(reference, sequence_calls(11, "(T)7\u21928")),
               "more preceding repeat copies")
  expect_error(.breseq_apply_sequence_edits(reference, sequence_calls(11, "(A)8\u21929")),
               "repeat sequence/count")
  expect_error(.breseq_apply_sequence_edits(reference, sequence_calls(11, "(T)8\u21928")),
               "copy numbers")
})

test_that("duplicates remain auditable and whole-replicon DEL returns no sequence", {
  calls <- sequence_calls(c("2", "2", "4"), c("C\u2192A", "C->A", "+T"))
  result <- .breseq_apply_sequence_edits("ACGT", calls)
  expect_identical(result$sequence, "AAGTT")
  expect_equal(nrow(result$edits), 2L)
  expect_identical(result$edits$source_rows, c("1,2", "3"))
  expect_identical(result$report$status, c("applied", "duplicate", "applied"))
  expect_identical(result$report$edit_index, c(1L, 1L, 2L))
  lost <- .breseq_apply_sequence_edits(strrep("A", 1000), sequence_calls("1", "\u03941,000 bp"))
  expect_identical(lost$sequence, "")
  expect_identical(lost$edits$new_start, 1)
  expect_identical(lost$edits$new_end, 0)
  unchanged <- .breseq_apply_sequence_edits("ACGTN", data.frame())
  expect_identical(unchanged$sequence, "ACGTN")
  expect_equal(nrow(unchanged$edits), 0L)
})

test_that("unresolved edits stop strictly and never return a partial genome", {
  calls <- sequence_calls(c("2", "3", "4", "8", "1=", "5"),
    c("C\u2192A", "A\u2192T", "+300 bp", "\u03943 bp", "G\u2192A", "IS1 (+)"))
  inspected <- .breseq_apply_sequence_edits("ACGTACGT", calls, strict = FALSE)
  expect_false(inspected$complete)
  expect_null(inspected$sequence)
  expect_identical(inspected$report$status,
    c("validated", "reference_mismatch", "unsupported", "invalid_position", "invalid_position", "unsupported"))
  error <- tryCatch(.breseq_apply_sequence_edits("ACGTACGT", calls), error = identity)
  expect_s3_class(error, "breseq_sequence_edit_error")
  expect_identical(error$report, inspected$report)
  expect_null(error$sequence)
  expect_error(.breseq_apply_sequence_edits("ACGT", sequence_calls(1, "A\u2192A")), "does not change")
  expect_error(.breseq_apply_sequence_edits("ACGT", sequence_calls(2, "-AA")), "Reference allele mismatch")
  expect_error(.breseq_apply_sequence_edits("ACGT", sequence_calls(0, "+A")), "outside the reference")
  expect_error(.breseq_apply_sequence_edits("ACGT", sequence_calls(1, "\u03940 bp")), "positive integer")
  expect_error(.breseq_apply_sequence_edits("ACGT", sequence_calls(c(1, 2), c("A\u2192G", "C\u2192A"),
    seq_id = c("chr", "plasmid"))), "multiple replicons")
})

test_that("conflicting edits are rejected while distinct boundary edits are unambiguous", {
  conflicts <- list(
    sequence_calls(c(2, 3), c("\u03943 bp", "G\u2192A")),
    sequence_calls(c(2, 3), c("\u03943 bp", "+A")),
    sequence_calls(c(2, 2), c("+A", "+G")),
    sequence_calls(c(2, 2), c("C\u2192A", "C\u2192T"))
  )
  for (calls in conflicts) {
    result <- .breseq_apply_sequence_edits("ACGTACGT", calls, strict = FALSE)
    expect_identical(result$report$status, c("overlap", "overlap"))
    expect_null(result$sequence)
    expect_error(.breseq_apply_sequence_edits("ACGTACGT", calls), "Overlapping edits")
  }
  boundary <- .breseq_apply_sequence_edits("ACGT", sequence_calls(c(2, 2, 3), c("C\u2192A", "+T", "\u03941 bp")))
  expect_identical(boundary$sequence, "AATT")
  expect_identical(boundary$edits$new_start, c(2, 3, 4))
})

test_that("synthetic reconstruction agrees with independent gdtools APPLY", {
  executable <- unname(Sys.which("gdtools"))
  skip_if(!nzchar(executable), "gdtools is optional for cross-validation")
  directory <- tempfile(); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  reference <- paste0("ACG", strrep("T", 8), "G", strrep("AC", 3), "GCCATGG")
  fasta <- file.path(directory, "reference.fasta")
  gd <- file.path(directory, "calls.gd")
  output <- file.path(directory, "mutated.fasta")
  log <- file.path(directory, "gdtools.log")
  writeLines(c(">fixture", reference), fasta)
  # Literal independent GenomeDiff calls, not generated from the engine's edits.
  writeLines(c("#=GENOME_DIFF 1.0", "SNP\t1\t.\tfixture\t2\tA", "INS\t2\t.\tfixture\t11\tT",
    "DEL\t3\t.\tfixture\t17\t2", "SUB\t4\t.\tfixture\t20\t2\tTA",
    "INS\t5\t.\tfixture\t25\tGG"), gd)
  status <- system2(executable, c("APPLY", "-f", "FASTA", "-r", shQuote(fasta),
                                 "-o", shQuote(output), shQuote(gd)), stdout = log, stderr = log)
  expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  lines <- readLines(output, warn = FALSE)
  independently_applied <- toupper(paste0(lines[!grepl("^>", lines)], collapse = ""))
  calls <- sequence_calls(c(20, 11, 2, 17, 25),
    c("2 bp\u2192TA", "(T)8\u21929", "C\u2192A", "(AC)3\u21922", "+GG"))
  expect_identical(.breseq_apply_sequence_edits(reference, calls)$sequence, independently_applied)
})
