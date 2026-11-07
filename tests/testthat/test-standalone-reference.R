test_that("standalone GenBank references retain every replicon and its own GC", {
  directory <- tempfile("standalone-genbank-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  reference <- file.path(directory, "multi.gbff")
  original <- c(
    "LOCUS       chromosome_alias 12 bp DNA circular BCT 01-JAN-2026",
    "DEFINITION  First chromosome",
    "            with a continued definition.",
    "ACCESSION   NZ_TEST001", "VERSION     NZ_TEST001.1",
    "FEATURES             Location/Qualifiers",
    "     CDS             complement(join(1..3,7..9))",
    "ORIGIN", "        1 ccccggggaaaa", "//",
    "LOCUS       plasmid_alias 12 bp DNA linear BCT 01-JAN-2026",
    "DEFINITION  Second replicon.",
    "ACCESSION   NZ_TEST002", "VERSION     NZ_TEST002.4",
    "ORIGIN", "        1 aaaattccaaaa", "//"
  )
  writeLines(original, reference)
  before <- tools::md5sum(reference)
  output_dir <- file.path(directory, "derived", "references")
  prepared <- .breseq_prepare_standalone_reference(reference, output_dir)
  expect_identical(prepared$info$seq_id, c("NZ_TEST001", "NZ_TEST002"))
  expect_identical(prepared$info$length, c(12L, 12L))
  expect_identical(prepared$info$topology, c("circular", "linear"))
  expect_identical(prepared$info$description[1L],
                   "First chromosome with a continued definition.")
  expect_identical(names(prepared$replicons), prepared$info$seq_id)
  expect_identical(tools::md5sum(reference), before)
  expect_identical(readLines(reference), original)

  gc <- lapply(prepared$replicons, function(replicon) {
    expect_true(file.exists(replicon$reference_file))
    expect_identical(dirname(replicon$reference_file), normalizePath(output_dir))
    expect_equal(nrow(replicon$genome_info), 1L)
    expect_null(replicon$gff_data)
    expect_length(grep("^LOCUS", readLines(replicon$reference_file)), 1L)
    suppressMessages(calculate_gc_content(replicon$reference_file, window_size = 4))
  })
  expect_identical(unique(gc[[1L]]$seqid), "NZ_TEST001")
  expect_identical(unique(gc[[2L]]$seqid), "NZ_TEST002")
  expect_identical(gc[[1L]]$gc_ratio, c(1, 1, 0))
  expect_identical(gc[[2L]]$gc_ratio, c(0, 0.5, 0))
  expect_true(any(grepl("complement(join(1..3,7..9))",
    readLines(prepared$replicons[[1L]]$reference_file), fixed = TRUE)))
})

test_that("standalone FASTA references preserve topology, bases, and record order", {
  directory <- tempfile("standalone-fasta-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  reference <- file.path(directory, "multi.fna")
  writeLines(c(">chr.2 chromosome [topology=circular]", " ccCC GGGG ", "aaaa",
               ">plasmid.1 plasmid [topology=linear]", "aaaa ttcc aaaa"), reference)
  before <- tools::md5sum(reference)
  prepared <- .breseq_prepare_standalone_reference(reference, file.path(directory, "split"))
  expect_identical(prepared$info$seq_id, c("chr", "plasmid"))
  expect_identical(prepared$info$topology, c("circular", "linear"))
  expect_identical(prepared$info$length, c(12L, 12L))
  expect_identical(prepared$info$description[2L], "plasmid [topology=linear]")
  gc <- lapply(prepared$replicons, function(replicon) {
    suppressMessages(calculate_gc_content(replicon$reference_file, 4))
  })
  expect_identical(gc[[1L]]$gc_ratio, c(1, 1, 0))
  expect_identical(gc[[2L]]$gc_ratio, c(0, 0.5, 0))
  expect_identical(unique(gc[[2L]]$seqid), "plasmid")
  expect_identical(tools::md5sum(reference), before)

  writeLines(c(">ambiguous", "ACGTRYSWKMBDHVN"), reference)
  ambiguous <- .breseq_prepare_standalone_reference(reference, file.path(directory, "split"))
  expect_identical(ambiguous$info$length, 15L)
  expect_identical(ambiguous$info$topology, "circular")
  expect_identical(readLines(ambiguous$replicons[[1L]]$reference_file)[2L],
                   "ACGTRYSWKMBDHVN")
})

test_that("standalone reference fallback IDs are supported and invalid records rejected", {
  directory <- tempfile("standalone-reference-errors-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  reference <- file.path(directory, "fallback.gbk")
  writeLines(c("LOCUS       fallback_locus 4 bp DNA BCT 01-JAN-2026",
               "ACCESSION   ACCESSION_ONLY.2", "ORIGIN", "        1 acgt", "//",
               "LOCUS       LOCAL_ID.3 4 bp DNA BCT 01-JAN-2026",
               "ORIGIN", "        1 aacc", "//"), reference)
  prepared <- .breseq_prepare_standalone_reference(reference, file.path(directory, "split"))
  expect_identical(prepared$info$seq_id, c("ACCESSION_ONLY", "LOCAL_ID"))
  expect_identical(prepared$info$topology, c("unknown", "unknown"))

  writeLines(c("LOCUS       bad 5 bp DNA circular BCT 01-JAN-2026",
               "ORIGIN", "        1 acgt", "//"), reference)
  expect_error(.breseq_prepare_standalone_reference(reference, directory), "length does not match")
  writeLines(c("LOCUS       bad 4 bp DNA circular BCT 01-JAN-2026",
               "ORIGIN", "        1 acgt"), reference)
  expect_error(.breseq_prepare_standalone_reference(reference, directory), "terminator")

  fasta <- file.path(directory, "bad.fasta")
  writeLines(c(">same.1", "ACGT", ">same.2", "ACGT"), fasta)
  expect_error(.breseq_prepare_standalone_reference(fasta, directory), "duplicate sequence IDs")
  writeLines(c(">empty", ">next", "ACGT"), fasta)
  expect_error(.breseq_prepare_standalone_reference(fasta, directory), "no sequence")
  writeLines(c(">protein", "ACGTE"), fasta)
  expect_error(.breseq_prepare_standalone_reference(fasta, directory), "nucleotide sequence")
})
