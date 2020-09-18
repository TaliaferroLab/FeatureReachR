context("get length or GC Outputs")
library(RNAreachR)

case <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Enhanced.fasta", package = "RNAreachR"))

L <- get_length(case)
GC <- get_GC(case)

test_that("Output is correct dimension and types", {
  expect_is(L, "data.frame")
  expect_is(L$length, "integer")
  expect_is(L$gene, "factor")
  expect_equal(nrow(L), length(case))
  expect_equal(ncol(L), 2)

  expect_is(GC, "data.frame")
  expect_is(GC$GC, "numeric")
  expect_is(GC$gene, "factor")
  expect_equal(nrow(GC), length(case))
  expect_equal(ncol(GC), 2)

})

test_that("Output produces expected values", {
  expect_equal(any(is.na(L)), FALSE)
  expect_equal(all(names(case) == L$gene), TRUE)
  expect_equal(L$length[1:10], c(200, 200, 200, 200, 200, 200, 200, 200, 200, 200))

  expect_equal(any(is.na(GC)), FALSE)
  expect_equal(all(names(case) == GC$gene), TRUE)
  expect_equal(GC$GC[1:10], c(0.340, 0.310, 0.290, 0.620, 0.280, 0.630, 0.595, 0.400, 0.600, 0.535))

})
