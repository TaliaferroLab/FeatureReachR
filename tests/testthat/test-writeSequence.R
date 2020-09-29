context("writeSequence")
library(FeatureReachR)

hs_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.v33.annotation.gff3.gz", package = "FeatureReachR"))
hs_case_tx <- c("ENST00000456328.2", "ENST00000338338.9", "ENST00000356026.10", "ENST00000607222.1", "ENST00000342066.8")
write_Sequence(hs_filtered_TxDb, hs_case_tx, "CDS", "testdata/hs_CDS_test", "both")
write_Sequence(hs_filtered_TxDb, hs_case_tx, "whole", "testdata/hs_whole_test", "gff3")
write_Sequence(hs_filtered_TxDb, hs_case_tx, "promoter", promoter_size = 7500, "testdata/hs_promoter_test", "gff3")

mm_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"))
mm_case_tx <- c("ENSMUST00000159265.1", "ENSMUST00000027032.5", "ENSMUST00000130201.7", "ENSMUST00000157375.1")
write_Sequence(mm_filtered_TxDb, mm_case_tx, "UTR5", "testdata/mm_UTR5_test", "both")
write_Sequence(mm_filtered_TxDb, mm_case_tx, "UTR3", "testdata/mm_UTR3_test", "fa")
write_Sequence(mm_filtered_TxDb, mm_case_tx, "promoter", promoter_size = 1500, "testdata/mm_promoter_test", "fa")

hs_CDS_test_fasta <- Biostrings::readDNAStringSet("testdata/hs_CDS_test.fa")
hs_CDS_test_gff <- rtracklayer::import("testdata/hs_CDS_test.gff3")
hs_whole_test_gff <- rtracklayer::import("testdata/hs_whole_test.gff3")
hs_promoter_test_gff <- rtracklayer::import("testdata/hs_promoter_test.gff3")
mm_UTR5_test_fasta <- Biostrings::readDNAStringSet("testdata/mm_UTR5_test.fa")
mm_UTR5_test_gff <- rtracklayer::import("testdata/mm_UTR5_test.gff3")
mm_UTR3_test_fasta <- Biostrings::readDNAStringSet("testdata/mm_UTR3_test.fa")
mm_promoter_test_fasta <- Biostrings::readDNAStringSet("testdata/mm_promoter_test.fa")


hs_CDS_ctrl_fasta <- Biostrings::readDNAStringSet("testdata/200916_hs_CDS_output.fa")
hs_CDS_ctrl_gff <- rtracklayer::import("testdata/200916_hs_CDS_output.gff3")
hs_whole_ctrl_gff <- rtracklayer::import("testdata/200916_hs_whole_output.gff3")
hs_promoter_ctrl_gff <- rtracklayer::import("testdata/200916_hs_promoter_output.gff3")
mm_UTR5_ctrl_fasta <- Biostrings::readDNAStringSet("testdata/200916_mm_UTR5_output.fa")
mm_UTR5_ctrl_gff <- rtracklayer::import("testdata/200916_mm_UTR5_output.gff3")
mm_UTR3_ctrl_fasta <- Biostrings::readDNAStringSet("testdata/200916_mm_UTR3_output.fa")
mm_promoter_ctrl_fasta <- Biostrings::readDNAStringSet("testdata/200916_mm_promoter_output.fa")


test_that("getTxout Output is expected", {
  expect_equal(hs_CDS_test_fasta, hs_CDS_ctrl_fasta)
  expect_equal(hs_CDS_test_gff, hs_CDS_ctrl_gff)
  expect_equal(hs_whole_test_gff, hs_whole_ctrl_gff)
  expect_equal(hs_promoter_test_gff, hs_promoter_ctrl_gff)
  expect_equal(mm_UTR5_test_fasta, mm_UTR5_ctrl_fasta)
  expect_equal(mm_UTR5_test_gff, mm_UTR5_ctrl_gff)
  expect_equal(mm_UTR3_test_fasta, mm_UTR3_ctrl_fasta)
  expect_equal(mm_promoter_test_fasta, mm_promoter_ctrl_fasta)

})


mm_case_genes <- c("ENSMUSG00000097392", "ENSMUSG00000025607", "ENSMUSG00000030671", "ENSMUSG00000034764", "ENSMUSG00000116215", "ENSMUSG00000039556", "ENSMUSG00000066510", "ENSMUSG00000018160", "ENSMUSG00000114306", "ENSMUSG00000028277", "ENSMUSG00000037216", "ENSMUSG00000032299")
hs_case_genes <- c("ENSG00000000419.12", "ENSG00000001167.14", "ENSG00000000938.13")
mm_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"), filter = TRUE, protein.coding = FALSE)
mm_noCDS_tx <- c("ENSMUST00000223470.2", "ENSMUST00000172461.1", "ENSMUST00000183867.6")
unknown_tx <- c("FBtr0078434", "FBtr0073382", "FBtr0344279", "FBtr0076263", "FBtr0333391", "FBtr0078287")
unknown_TxDb <- AnnotationDbi::loadDb("testdata/dmTxDb")


test_that("getTxout detects incorrect inputs", {
  expect_error(write_Sequence(mm_filtered_TxDb, mm_case_genes, "CDS", "mydata/mm_test", "both"), "Please ensure the transcript list contains transcript IDs and not gene IDs.")
  expect_error(write_Sequence(mm_filtered_TxDb, hs_case_tx, "CDS", "mydata/mm_test", "both"), "Please ensure the transcript list and gff are of the same species.")
  expect_error(write_Sequence(mm_filtered_TxDb, mm_noCDS_tx, "CDS", "mydata/mm_test", "both"), "no sequences of type \"CDS\" found in transcript list")
  expect_error(write_Sequence(mm_filtered_TxDb, hs_case_tx, "CDS", "mydata/mm_test", "both"), "Please ensure the transcript list and gff are of the same species.")
  expect_error(write_Sequence(mm_filtered_TxDb, mm_case_tx, "CDS", "mydata/mm_test", "Both"), "not appropriate format, please use \"fa\", \"gff3\" or \"both\"")
  expect_error(write_Sequence(unknown_TxDb, unknown_tx, "CDS", "mydata/mm_test", "both"), "Species type is unknown. Please ensure gff and transcript list are either mouse or human")
  expect_error(write_Sequence(mm_filtered_TxDb, mm_case_tx, "CSD", "mydata/mm_test", "both"), "not appropriate seq type. please use one of \"whole\", \"CDS\", \"UTR5\", \"UTR3\" or, \"promoter\"")

})

