context("kmer and motif by gene Outputs")
library(FeatureReachR)

case <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Enhanced.fasta", package = "FeatureReachR"))

k5 <- kmer_by_gene(case[1:10], 5)
k6 <- kmer_by_gene(case[1:10], 6)

test_that("kmer by gene Output is correct dimension and types", {
  expect_is(k5, "data.frame")
  expect_is(k5$gene, "character")
  expect_equal(all(sapply(k5[,-1], class) == "integer"), TRUE)
  expect_equal(ncol(k5), 4^5+1)
  expect_equal(nrow(k5), length(case[1:10]))

  expect_is(k6, "data.frame")
  expect_is(k6$gene, "character")
  expect_equal(all(sapply(k6[,-1], class) == "integer"), TRUE)
  expect_equal(ncol(k6), 4^6+1)
  expect_equal(nrow(k6), length(case[1:10]))

})

test_that("kmer by gene Output produces expected values", {
  expect_equal(all(nchar(colnames(k5[,-1])) == 5), TRUE)
  expect_equal(all(names(case[1:10]) %in% k5$gene), TRUE)
  expect_equal(k5$AAAAA[1:10], c(4, 0, 2, 0, 0, 0, 1, 0, 0, 0))
  expect_equal(k5$AAAAC[1:10], c(1, 1, 0, 0, 1, 0, 1, 0, 0, 0))

  expect_equal(all(nchar(colnames(k6[,-1])) == 6), TRUE)
  expect_equal(all(names(case[1:10]) %in% k6$gene), TRUE)
  expect_equal(k6$AAAAAA[1:10], c(2, 0, 1, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(k6$AAAAAC[1:10], c(1, 0, 0, 0, 0, 0, 1, 0, 0, 0))
})

C <- motif_by_gene(CISBPRNA_hs_PWM[1:10], case[1:10])
R <- motif_by_gene(RBNS_PWM[1:10], case[1:10])

test_that("motif by gene Output is correct dimension and types", {
  expect_is(C, "data.frame")
  expect_is(C$motif, "character")
  expect_equal(all(sapply(C[,-1], class) == "integer"), TRUE)
  expect_equal(ncol(C), length(case[1:10])+1)
  expect_equal(nrow(C), length(CISBPRNA_hs_PWM[1:10]))

  expect_is(R, "data.frame")
  expect_is(R$motif, "character")
  expect_equal(all(sapply(R[,-1], class) == "integer"), TRUE)
  expect_equal(ncol(R), length(case[1:10])+1)
  expect_equal(nrow(R), length(RBNS_PWM[1:10]))

})

test_that("Motif by gene Output produces expected values", {
  expect_equal(all(colnames(C[,-1]) %in% names(case)), TRUE)
  expect_equal(all(names(CISBPRNA_hs_PWM[1:10]) %in% C$motif), TRUE)
  expect_equal(dplyr::pull(C[c(1:10),2]), c(0, 0, 1, 1, 0, 2, 0, 0, 3, 1))
  expect_equal(dplyr::pull(C[c(1:10),3]), c(2, 0, 2, 3, 0, 1, 1, 0, 0, 1))

  expect_equal(all(colnames(R[,-1]) %in% names(case)), TRUE)
  expect_equal(all(names(RBNS_PWM[1:10]) %in% R$motif), TRUE)
  expect_equal(dplyr::pull(R[c(1:10),2]), c(1, 1, 0, 1, 0, 2, 2, 2, 0, 1))
  expect_equal(dplyr::pull(R[c(1:10),3]), c(2, 8, 3, 4, 6, 2, 2, 4, 7, 1))
})

PWMnonList <- RBNS_PWM[[1]]
long_PWM <- list(one = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("ACACA", "CACAC", "AAACA"))),
                 two = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("TTTAA", "TTTTAA", "ATTAA"))),
                 three = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("GCATA", "CCATAA", "GGGAT"))))

test_that("motif by gene detects incorrect inputs", {
  expect_error(motif_by_gene(PWMnonList, case), "PWM_list must be a list")
  expect_error(motif_by_gene(long_PWM, case), "Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")", fixed = TRUE)

})
