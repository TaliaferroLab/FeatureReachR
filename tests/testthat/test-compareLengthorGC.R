context("Compare Outputs")
library(RNAreachR)

case <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Enhanced.fasta", package = "RNAreachR"))
ctrl <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Control.fasta", package = "RNAreachR"))
L <- length_compare(case, ctrl)
GC <- GC_compare(case, ctrl)

test_that("Output is correct dimension and types", {
  expect_is(L, "data.frame")
  expect_is(L$wilcox.p, "numeric")
  expect_is(L$mean_case, "numeric")
  expect_is(L$mean_ctrl, "numeric")
  expect_is(L$mean_FC, "numeric")
  expect_is(L$CliffDelta, "numeric")
  expect_is(L$lowerCD, "numeric")
  expect_is(L$upperCD, "numeric")
  expect_equal(nrow(L), 1)
  expect_equal(ncol(L), 7)

  expect_is(GC, "data.frame")
  expect_is(GC$wilcox.p, "numeric")
  expect_is(GC$mean_case, "numeric")
  expect_is(GC$mean_ctrl, "numeric")
  expect_is(GC$mean_FC, "numeric")
  expect_is(GC$CliffDelta, "numeric")
  expect_is(GC$lowerCD, "numeric")
  expect_is(GC$upperCD, "numeric")
  expect_equal(nrow(GC), 1)
  expect_equal(ncol(GC), 7)

})

test_that("Output produces expected values", {
  expect_equal(any(is.na(L)), FALSE)
  expect_equal(L$wilcox.p, 0.3577381, tolerance = .0001)
  expect_equal(L$mean_case,  198.3822, tolerance = .0001)
  expect_equal(L$mean_ctrl, 198.2121, tolerance = .0001)
  expect_equal(L$mean_FC, 1.000858, tolerance = .0001)
  expect_equal(L$CliffDelta, 0.007511896)
  expect_equal(L$lowerCD, -0.007366476)
  expect_equal(L$upperCD, 0.02238694)

  expect_equal(any(is.na(GC)), FALSE)
  expect_equal(GC$wilcox.p, 0.04549015, tolerance = .0001)
  expect_equal(GC$mean_case, 0.4575343, tolerance = .0001)
  expect_equal(GC$mean_ctrl, 0.4492366, tolerance = .0001)
  expect_equal(GC$mean_FC, 1.018471, tolerance = .0001)
  expect_equal(GC$CliffDelta, 0.04716512)
  expect_equal(GC$lowerCD, 0.001848139)
  expect_equal(GC$upperCD,  0.09228878)

})

test_that("GC or length compare detects incorrect inputs", {
  expect_warning(GC_compare(case,case), "some sequences in case set are also in the control set. This is not recommended.")
  expect_warning(length_compare(case,case), "some sequences in case set are also in the control set. This is not recommended.")

})
