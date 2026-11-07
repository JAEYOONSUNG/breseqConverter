test_that("comma-formatted deletion lengths retain their full size", {
  reported <- c("\u03945,886 bp", "\u039412,411 bp", "\u039429,170 bp", "\u039449 bp", "\u039450 bp")
  original <- reported
  expect_identical(classify_mutations(reported),
                   c("large_deletion", "large_deletion", "large_deletion", "deletion", "large_deletion"))
  expect_identical(reported, original)
})

test_that("oriented junction coordinates are parsed without changing their identity", {
  reported <- c("186731\u00a0=", "=\u00a0186733", "1668721 =", "= 1668733",
                "> 1,234", "1,234 <", "= 2,000 =")
  original <- reported
  parsed <- lapply(reported, parse_position)
  expected <- c(186731L, 186733L, 1668721L, 1668733L, 1234L, 1234L, 2000L)
  expect_identical(vapply(parsed, `[[`, integer(1), "start"), expected)
  expect_identical(vapply(parsed, `[[`, integer(1), "end"), expected)
  expect_identical(reported, original)
  expect_false(identical(reported[5], reported[6]))
})

test_that("ordinary position ranges remain supported and malformed text stays missing", {
  expect_identical(parse_position("1,234\u20132,345"), list(start = 1234L, end = 2345L))
  expect_identical(parse_position("1234 - 2345"), list(start = 1234L, end = 2345L))
  expect_identical(parse_position("1234-"), list(start = 1234L, end = 1234L))
  expect_identical(parse_position(1234), list(start = 1234L, end = 1234L))
  for (missing in list(NA_character_, NULL, "", "unknown", "chr1:1234", "1234-2345-3456")) {
    expect_identical(parse_position(missing), list(start = NA_integer_, end = NA_integer_))
  }
})

test_that("wide classification uses the first observed mutation in each row", {
  data <- data.frame(position = c(1, 2, 3), check.names = FALSE)
  data[["mutation.S1/output/index"]] <- c(NA, "A\u2192G", "")
  data[["annotation.S1/output/index"]] <- c(NA, "coding (1/9 nt)", "")
  data[["mutation.S2/output/index"]] <- c("\u039455,393 bp", "C\u2192T", "+A")
  data[["annotation.S2/output/index"]] <- c("", "intergenic (+1/-2)", "pseudogene (1/9 nt)")
  original <- data
  result <- add_mutation_types(data)
  expect_identical(result$mutation_type, c("large_deletion", "SNP", "insertion"))
  expect_identical(result[names(original)], original)
  expect_identical(data, original)
})

test_that("public classification and GenBank palettes remain separate", {
  expect_identical(classify_mutation("\u03945,886 bp", "coding"), "large_deletion")
  expect_identical(unname(get_mutation_color("large_deletion")), "#c0392b")
  expect_identical(.classify_genbank_mutation("\u03945,886 bp"), "DEL")
  expect_identical(unname(.get_genbank_mutation_color("DEL")), "#03045E")
  expect_identical(unname(.get_genbank_mutation_color("large_deletion")), "#03045E")
})
