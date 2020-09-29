context("motif estimates from kmer")
library(FeatureReachR)

R4 <- relate_kmer_PWM(4, RBNS_PWM)
RBNS_5mer <- RBNS_PWM[lapply(RBNS_PWM, ncol)==5]
R5 <- relate_kmer_PWM(5, RBNS_5mer[1:5])

test_that("kmer relation Output is correct dimension and types", {
  expect_is(R4, "data.frame")
  expect_is(R4$motif, "character")
  expect_equal(all(sapply(R4[,-1], class) == "numeric"), TRUE)
  expect_equal(ncol(R4), 4^4+1)
  expect_equal(nrow(R4), 348)

  expect_is(R5, "data.frame")
  expect_is(R5$motif, "character")
  expect_equal(all(sapply(R5[,-1], class) == "numeric"), TRUE)
  expect_equal(ncol(R5), 4^5+1)
  expect_equal(nrow(R5), 5)

})

test_that("kmer relation Output produces expected values", {
  expect_equal(all(nchar(colnames(R4[,-1])) == 4), TRUE)
  expect_equal(all(names(RBNS_PWM) %in% substr(R4$motif, 1, nchar(R4$motif)-2)), TRUE)
  expect_equal(dplyr::pull(R4[1:50,5]), c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(dplyr::pull(R4[1:50,14]), c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

  expect_equal(all(nchar(colnames(R5[,-1])) == 5), TRUE)
  expect_equal(all(names(RBNS_5mer[1:5]) %in% R5$motif), TRUE)
  expect_equal(as.numeric(R5[1,2:50]), c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0))
  expect_equal(as.numeric(R5[3,2:50]), c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

})

PWMnonList <- RBNS_PWM[[1]]
RNA_PWM <- RBNS_PWM[1:3]
rownames(RNA_PWM[[1]]) <- c("A", "C", "G", "U")
rownames(RNA_PWM[[2]]) <- c("A", "C", "G", "U")
rownames(RNA_PWM[[3]]) <- c("A", "C", "G", "U")
long_PWM <- list(one = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("ACACA", "CACAC", "AAACA"))),
                 two = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("TTTAA", "TTTTAA", "ATTAA"))),
                 three = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("GCATA", "CCATAA", "GGGAT"))))

test_that("kmer relation detects incorrect inputs", {
  expect_error(relate_kmer_PWM(4, PWMnonList), "PWM_list must be a list")
  expect_error(relate_kmer_PWM(4, RNA_PWM), "'rownames(pwm)' must be the 4 DNA bases ('DNA_BASES')", fixed = TRUE)
  expect_error(relate_kmer_PWM(4, long_PWM), "Error: Ensure all PWM matricies have exactly 4 rows (A, C, G and T)", fixed = TRUE)

})

enriched_sixmers <- c("AAGGAA", "ACACAC", "AGAAGG", "AGAGAG", "AGAGGG", "AGGAAG", "AGGAGG", "AGGGAG", "CACACA", "GAAGGA", "GAGAAG", "GAGAGA", "GAGGAG", "GAGGGA", "GAGGGG", "GGAAGG", "GGAGGA", "GGAGGG", "GGGAGG")
E <- estimate_motif_from_kmer(enriched_sixmers, "RBNS")

test_that("motif by gene Output is correct dimension and types", {
  expect_is(E, "data.frame")
  expect_is(E$motif, "character")
  expect_is(E$input_kmer, "numeric")
  expect_is(E$all_kmer, "numeric")
  expect_is(E$input_freq, "numeric")
  expect_is(E$all_freq, "numeric")
  expect_is(E$log2FC, "numeric")
  expect_is(E$p_val, "numeric")
  expect_is(E$p_adj, "numeric")
  expect_equal(ncol(E), 8)
  expect_equal(nrow(E), length(RBNS_PWM))

})

test_that("kmer estimate Output produces expected values", {
  expect_equal(all(E$motif %in% names(RBNS_PWM)), TRUE)
  expect_equal(E$input_kmer[1:10], c(13, 14, 12, 8, 6, 6, 6, 4, 6, 5))
  expect_equal(E$all_kmer[1:10], c(147, 228, 120, 100, 72, 80, 92, 31, 128, 118))
  expect_equal(E$input_freq[1:10], c(0.6842105, 0.7368421, 0.6315789, 0.4210526, 0.3157895, 0.3157895, 0.3157895, 0.2105263, 0.3157895, 0.2631579), tolerance = 0.0001)
  expect_equal(E$all_freq[1:10], c(0.0328673, 0.05248958, 0.02649007, 0.02256561, 0.01618837, 0.0181506, 0.02109394, 0.006622517, 0.02992396, 0.02771646), tolerance = 0.0001)
  expect_equal(E$log2FC[1:10], c(4.447405, 3.910529, 4.634716, 4.297121, 4.358057, 4.2015, 3.997349, 5.035161, 3.530203, 3.391581), tolerance = 0.0001)
  expect_equal(E$p_val[1:10], c(6.738581e-15, 3.916085e-14, 3.16627e-14, 1.070571e-08, 8.752772e-07, 1.567161e-06, 3.394179e-06, 1.581993e-05, 2.094496e-05, 0.0001866864), tolerance = 0.0001)
  expect_equal(E$p_adj[1:10], c(8.827541e-13, 1.710024e-12, 1.710024e-12, 3.50612e-07, 2.293226e-05, 3.421634e-05, 6.351963e-05, 0.0002590513, 0.0003048656, 0.002445592), tolerance = 0.0001)

})

mixed_kmers <- c("AACCAA", "AACC", "ACCCA", "ACCCCAA")

test_that("kmer estimate detects incorrect inputs", {
  expect_error(estimate_motif_from_kmer(mixed_kmers, "RBNS"), "kmers in kmer list are not the same length")
  expect_error(estimate_motif_from_kmer(enriched_sixmers, "Custom"), "motif_set must be either \"CISBPRNA_mm\", \"CISBPRNA_hs\", \"JASPAR_mm\", \"JASPAR_hs\", \"RBNS\" or \"custom\"", fixed = TRUE)
  expect_error(estimate_motif_from_kmer(enriched_sixmers, "custom"), "A custom motif_by_kmer_matrix is required for estimating custom motif occurance from kmers \n\n         use ouput from relate_kmer_PWM()", fixed = TRUE)

})

enriched_fourmers <- c("AGAC", "ATAT", "CAGT", "GGGG", "GCCG", "TGCT", "GGAG", "CTGA", "TTTT", "TACT", "ACTT", "AACC", "AAGG", "CCAT", "GAGA", "GTGT", "TGGT")
enriched_fivemers <- c("AGTAC", "CATAT", "CAAGT", "GGTGG", "GCGCG", "TGCCT", "GGAGA", "CTGAT", "TTTAT", "CTACT", "ATCTT", "GAACC", "AACGG", "CCCAT", "GAGAA", "GATGT", "TGGTT")
enriched_sevenmers <- c("AGTACTA", "CAGCTAT", "CCCAAGT", "GGTATGG", "GCGCGCG", "TGTACCT", "GGCCAGA", "CTTAGAT", "TTTTTTT", "TTAGACT", "ACCATTT", "AACCACC", "AATAAGG", "CCATCAT", "GGATAGA", "GTGTTAT", "TGCCCGT")

test_that("tables are accessible", {
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fourmers, "CISBPRNA_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fourmers, "CISBPRNA_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fourmers, "JASPAR_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fourmers, "JASPAR_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fourmers, "RBNS"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fivemers, "CISBPRNA_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fivemers, "CISBPRNA_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fivemers, "JASPAR_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fivemers, "JASPAR_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fivemers, "RBNS"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sixmers, "CISBPRNA_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sixmers, "CISBPRNA_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sixmers, "JASPAR_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sixmers, "JASPAR_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sixmers, "RBNS"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sevenmers, "CISBPRNA_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sevenmers, "CISBPRNA_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sevenmers, "JASPAR_mm"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sevenmers, "JASPAR_hs"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_sevenmers, "RBNS"))),FALSE)
  expect_equal(all(is.na(estimate_motif_from_kmer(enriched_fivemers, "custom", R5))),FALSE)
})

