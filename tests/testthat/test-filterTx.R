context("filter tx")
library(RNAreachR)

NF <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "RNAreachR"), filter = FALSE, protein.coding = FALSE)
Filt <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "RNAreachR"), filter = TRUE, protein.coding = FALSE)
PC <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "RNAreachR"), filter = TRUE, protein.coding = TRUE)

test_that("filter Tx Output is correct dimension and types", {
  expect_is(NF, "TxDb")
  expect_is(Filt, "TxDb")
  expect_is(PC, "TxDb")
  expect_equal(length(GenomicFeatures::transcripts(NF)), 138835)
  expect_equal(length(GenomicFeatures::transcripts(Filt)), 68086)
  expect_equal(length(GenomicFeatures::transcripts(PC)), 35483)
})

