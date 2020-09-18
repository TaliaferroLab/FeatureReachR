context("make length df")
library(RNAreachR)

hs_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.v33.annotation.gff3.gz", package = "RNAreachR"))
longest_hs <- make_longest_df(hs_filtered_TxDb)
median_hs <- make_median_df(hs_filtered_TxDb)

test_that("filter Tx Output is correct dimension and types", {
  expect_is(longest_hs, "data.frame")
  expect_is(median_hs, "data.frame")
  expect_true(all(sapply(longest_hs, class) == "character"))
  expect_true(all(sapply(median_hs, class) == "character"))
  expect_equal(nrow(longest_hs), 28321)
  expect_equal(ncol(longest_hs), 5)
  expect_equal(nrow(median_hs), 17987)
  expect_equal(ncol(median_hs), 5)
})

test_that("filter Tx Outputs correct values", {
  expect_equal(length(longest_hs$gene_id), unique(length(longest_hs$gene_id)))
  expect_equal(length(median_hs$gene_id), unique(length(median_hs$gene_id)))
  expect_equal(as.character(longest_hs[1,]), c("ENSG00000000003", "ENST00000373020.9", "ENST00000373020.9", "ENST00000373020.9", "ENST00000373020.9"))
  expect_equal(as.character(median_hs[1,]), c("ENSG00000000003", "ENST00000373020.9", "ENST00000373020.9", "ENST00000373020.9", "ENST00000373020.9"))

})
