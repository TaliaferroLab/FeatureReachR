context("get control, tx2gene and RNA breakdown")
library(FeatureReachR)

mm_case_genes <- c("ENSMUSG00000097392", "ENSMUSG00000025607", "ENSMUSG00000030671", "ENSMUSG00000034764", "ENSMUSG00000116215", "ENSMUSG00000039556", "ENSMUSG00000066510", "ENSMUSG00000018160", "ENSMUSG00000114306", "ENSMUSG00000028277", "ENSMUSG00000037216", "ENSMUSG00000032299")
mm_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"), filter = TRUE, protein.coding = FALSE)
longest_mm <- make_longest_df(mm_filtered_TxDb)
control_tx <- get_ctrl_tx(longest_mm, mm_case_genes, "UTR3")

test_that("get control Output is correct dimension and types", {
  expect_is(control_tx, "character")
  expect_equal(length(control_tx), 26719)
  expect_false(all(control_tx %in% mm_case_genes))

})

mm_case_tx <- c("ENSMUST00000159265.1", "ENSMUST00000027032.5", "ENSMUST00000130201.7", "ENSMUST00000157375.1")
hs_case_tx <- c("ENST00000456328.2", "ENST00000338338.9", "ENST00000478641.5", "ENST00000356026.10", "ENST00000607222.1", "ENST00000342066.8")
hs_case_genes <- c("ENSG00000000419.12", "ENSG00000001167.14", "ENSG00000000938.13")

test_that("get control detects incorrect inputs", {
  expect_error(get_ctrl_tx(longest_mm, mm_case_tx, "UTR3"), "Please ensure the gene list contains gene IDs and not transcript IDs.")
  expect_error(get_ctrl_tx(longest_mm, hs_case_genes, "UTR3"), "Please ensure the gene list and gff are of the same species.")

})

mm_RB <- RNA_breakdown(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"), mm_case_tx)
hs_RB <- RNA_breakdown(system.file("extdata", "gencode.v33.annotation.gff3.gz", package = "FeatureReachR"), hs_case_tx)
mm_madeup_Tx <- c("ENST01.1, ENST02.2", "ENST03.3")

test_that("RNA Breakdown Output is correct dimension and types", {
  expect_is(mm_RB, "data.frame")
  expect_is(mm_RB$transcript_type, "character")
  expect_is(mm_RB$count, "integer")
  expect_is(mm_RB$percent, "numeric")
  expect_equal(ncol(mm_RB), 3)
  expect_equal(sum(mm_RB$percent), 100, tolerance = 1)
  expect_equal(sum(mm_RB$count), length(mm_case_tx))

  expect_is(hs_RB, "data.frame")
  expect_is(hs_RB$transcript_type, "character")
  expect_is(hs_RB$count, "integer")
  expect_is(hs_RB$percent, "numeric")
  expect_equal(ncol(hs_RB), 3)
  expect_equal(sum(hs_RB$percent), 100, tolerance = 1)
  expect_equal(sum(hs_RB$count), length(hs_case_tx))
})

test_that("RNA Breakdown detects incorrect inputs", {
  expect_error(RNA_breakdown(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"), mm_case_genes), "Please ensure the transcript list contains transcript IDs and not gene IDs.")
  expect_error(RNA_breakdown(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"), hs_case_tx), "Please ensure the transcript list and gff are of the same species.")
  expect_error(RNA_breakdown(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"), mm_madeup_Tx), "Please ensure the transcript list and gff are of the same species.")
  expect_error(RNA_breakdown("testdata/partial_vM20.GFF3", mm_case_tx), "No transcript IDs in tx_list were not found in gff.")


})

mm_tx <- gene2Tx(longest_mm, mm_case_genes, "UTR3")

test_that("Tx2gene Output is correct dimension and types", {
  expect_is(mm_tx, "character")
  expect_equal(length(mm_tx), 12)

})

test_that("Tx2gene detects incorrect inputs", {
  expect_error(gene2Tx(longest_mm, mm_case_tx, "UTR3"), "Please ensure the gene list contains gene IDs and not transcript IDs.")
  expect_error(gene2Tx(longest_mm, hs_case_genes, "UTR3"), "Please ensure the gene list and gff are of the same species.")
})
