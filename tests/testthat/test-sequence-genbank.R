changed_genbank_fixture <- function() {
  sequence <- paste0("AAAA", "ATGAAATTTTAA", "CCCC", "TTAGGGCAT", "GGGG")
  features <- c(
    "     source          1..33", "                     /organism=\"Example bacterium\"",
    "     gene            5..16", "                     /locus_tag=\"g1\"",
    "     CDS             5..16", "                     /locus_tag=\"g1\"",
    "                     /product=\"reference product\"", "                     /transl_table=11",
    "                     /protein_id=\"REF1\"", "                     /translation=\"MKF\"",
    "     gene            complement(21..29)", "                     /locus_tag=\"g2\"",
    "     CDS             complement(21..29)", "                     /locus_tag=\"g2\"",
    "                     /transl_table=11", "                     /protein_id=\"REF2\"",
    "                     /translation=\"MP\"",
    "     tRNA            17..20", "                     /anticodon=(pos:18..20,aa:Leu)",
    "     misc_feature    join(2..3,30..32)", "                     /label=\"joined\"",
    "     misc_feature    3^4", "     misc_feature    <2..>3",
    "     misc_feature    21..23", "                     /label=\"fully deleted\"")
  record <- function(id, dna, f = character()) c(
    sprintf("LOCUS       %-16s %11d bp    DNA     circular BCT 01-JAN-2020", id, nchar(dna)),
    paste0("DEFINITION  Reference ", id, "."), paste0("ACCESSION   ", id),
    paste0("VERSION     ", id, ".1"), "SOURCE      Example bacterium",
    "  ORGANISM  Example bacterium", "            Bacteria.",
    "FEATURES             Location/Qualifiers", f,
    paste0("CONTIG      join(", id, ".1:1..", nchar(dna), ")"),
    "BASE COUNT  1 a 1 c 1 g 1 t", "ORIGIN", paste("        1", tolower(dna)), "//")
  file <- tempfile(fileext = ".gbff")
  writeLines(c(record("chr", sequence, features), record("pRaw", "ACGT"), record("pCov", "TGCA")), file)
  list(file = file, sequence = sequence)
}

test_that("changed GenBank files contain edited bases and lifted reference features", {
  fixture <- changed_genbank_fixture(); output <- tempfile()
  on.exit(unlink(c(fixture$file, output), recursive = TRUE), add = TRUE)
  calls <- data.frame(position = c(3, 8, 21, 1),
    mutation = c("+CC", "A\u2192G", "\u03943 bp", "\u03944 bp"),
    `seq id` = c("chr", "chr", "chr", "pRaw"), evidence = c("RA", "RA", "JC", "MC"),
    annotation = "coding", check.names = FALSE)
  # This plot-only DEL must never alter the exported sequence.
  display <- calls[1, ]; display$position <- 5; display$mutation <- "\u03943 bp"
  display$annotation <- "coverage-derived DEL display"; display$evidence <- "MC"
  calls <- rbind(calls, display)
  original <- calls
  loss <- data.frame(sample = "S1", seq_id = "pCov", evidence_origin = "COVERAGE_RATIO", depth_ratio = 0.001)
  result <- export_changed_sequences_to_genbank(fixture$file, list(S1 = calls), output,
    samples = c("WT", "S1"), replicon_losses = loss)
  expect_identical(calls, original)
  manifest <- result$manifest
  expect_equal(nrow(manifest), 6L)
  expect_true(all(manifest$status[manifest$sample == "WT"] == "reference_baseline"))
  chr <- manifest[manifest$sample == "S1" & manifest$seq_id == "chr", ]
  expect_equal(chr$new_length, 32)
  expect_true(file.exists(chr$file))
  expect_true(endsWith(chr$file, "_called_variants.gbk"))
  expect_true(file.exists(chr$fasta_file))
  record <- .breseq_gb_read_records(chr$file)[[1]]
  expected <- paste0(substr(fixture$sequence, 1, 3), "CC", substr(fixture$sequence, 4, 7),
                     "G", substr(fixture$sequence, 9, 20), substr(fixture$sequence, 24, 33))
  expect_identical(record$sequence, expected)
  expect_match(record$header[1], "32 bp", fixed = TRUE)
  locations <- vapply(record$features, `[[`, "", "location")
  expect_true(all(c("1..32", "7..18", "complement(<23..28)", "join(2..3,29..31)", "5^6", "<2..>3") %in% locations))
  features <- record$features
  g1 <- Filter(function(f) f$key == "CDS" && f$location == "7..18", features)[[1]]
  q1 <- paste(unlist(g1$qualifiers), collapse = "")
  expect_match(q1, '/translation="MEF"', fixed = TRUE)
  expect_false(grepl("/pseudo|/protein_id", q1))
  g2 <- Filter(function(f) f$key == "CDS" && grepl("23..28", f$location, fixed = TRUE), features)[[1]]
  q2 <- paste(unlist(g2$qualifiers), collapse = "")
  expect_match(q2, "/pseudo", fixed = TRUE)
  expect_false(grepl("/translation|/protein_id", q2))
  trna <- Filter(function(f) f$key == "tRNA", features)[[1]]
  expect_match(paste(unlist(trna$qualifiers), collapse = ""), "pos:20..22", fixed = TRUE)
  markers <- Filter(function(f) f$key == "misc_feature" &&
    any(grepl('/label="\\[', unlist(f$qualifiers))), features)
  expect_length(markers, 3L)
  expect_setequal(vapply(markers, `[[`, "", "location"), c("4..5", "10", "22"))
  marker_text <- vapply(markers, function(f) paste(trimws(unlist(f$qualifiers)), collapse = ""), "")
  expect_true(any(grepl('[SNP:S1] A\u2192G', marker_text, fixed = TRUE)))
  expect_true(any(grepl('/ApEinfo_fwdcolor="#139487"', marker_text, fixed = TRUE)))
  expect_true(any(grepl('/ApEinfo_fwdcolor="#FF5C8D"', marker_text, fixed = TRUE)))
  expect_true(any(grepl('/ApEinfo_fwdcolor="#03045E"', marker_text, fixed = TRUE)))
  expect_true(any(grepl("deleted reference interval 21..23", marker_text, fixed = TRUE)))
  expect_true(any(result$feature_audit$status == "dropped_fully_deleted"))
  expect_identical(manifest$status[manifest$sample == "S1" & manifest$seq_id == "pRaw"], "lost_confirmed")
  expect_identical(manifest$status[manifest$sample == "S1" & manifest$seq_id == "pCov"], "lost_by_coverage")
  lost <- manifest$sample == "S1" & manifest$seq_id != "chr"
  expect_true(all(manifest$new_length[lost] == 0 & manifest$file[lost] == ""))
  expect_equal(sum(result$edits$sample == "S1"), 4L)
})

test_that("FASTA reconstruction and unsupported calls are explicit and safe", {
  reference <- tempfile(fileext = ".fasta"); output <- tempfile()
  on.exit(unlink(c(reference, output), recursive = TRUE), add = TRUE)
  writeLines(c(">chr", "ACGTACGT"), reference)
  good <- data.frame(position = 2, mutation = "C\u2192A", seq_id = "chr")
  bad <- data.frame(position = 3, mutation = "A\u2192C", seq_id = "chr")
  result <- export_changed_sequences_to_genbank(reference, list(good = good, bad = bad), output)
  expect_identical(result$manifest$status, c("exported_called_variants", "unresolved_calls"))
  expect_identical(.breseq_gb_read_records(result$manifest$file[1])[[1]]$sequence, "AAGTACGT")
  expect_identical(result$manifest$file[2], "")
  expect_match(result$manifest$message[2], "Reference allele mismatch", fixed = TRUE)
  expect_true(any(result$call_report$status == "reference_mismatch"))
  expect_error(export_changed_sequences_to_genbank(reference, list(good = good), output,
    replicon_losses = data.frame(sample = "good", seq_id = "chr", evidence_origin = "COVERAGE_RATIO", depth_ratio = 0.2)),
    "below 0.1")
})

test_that("translation preserves strand and flags frameshift or unsupported genetic codes", {
  expect_identical(.breseq_gb_qualifier("note", "test"), paste0(strrep(" ", 21), '/note="test"'))
  expect_false(grepl(intToUtf8(92), .breseq_gb_qualifier("note", "test"), fixed = TRUE))
  node <- .breseq_gb_location("complement(join(1..3,4..9))")
  expect_identical(.breseq_gb_translation(node, "TTAGGGCAT", 11L), "MP")
  expect_error(.breseq_gb_translation(.breseq_gb_location("1..8"), "ATGAAATA", 11L), "multiple of three")
  expect_error(.breseq_gb_translation(.breseq_gb_location("1..9"), "ATGAAATAA", 4L), "unsupported")
  expect_error(.breseq_gb_translation(.breseq_gb_location("<1..9"), "ATGAAATAA", 11L), "partial")
})
