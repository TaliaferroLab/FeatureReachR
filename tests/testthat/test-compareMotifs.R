context("kmer compare Outputs")
library(RNAreachR)

case <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Enhanced.fasta", package = "RNAreachR"))
ctrl <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Control.fasta", package = "RNAreachR"))

cisbp <- motif_compare(CISBPRNA_hs_PWM[1:10], case, ctrl[1:1000])
RBNS <- motif_compare(RBNS_PWM[1:10], case, ctrl[1:1000])

test_that("Output is correct dimension and types", {
  expect_is(cisbp, "data.frame")
  expect_is(cisbp$motif, "character")
  expect_is(cisbp$ctrl, "integer")
  expect_is(cisbp$case, "integer")
  expect_is(cisbp$ctrl_freq, "numeric")
  expect_is(cisbp$case_freq, "numeric")
  expect_is(cisbp$log2FC, "numeric")
  expect_is(cisbp$ctrl_tot, "integer")
  expect_is(cisbp$case_tot, "integer")
  expect_is(cisbp$pval, "numeric")
  expect_is(cisbp$p_adj, "numeric")
  expect_equal(nrow(cisbp), length(CISBPRNA_hs_PWM[1:10]))
  expect_equal(ncol(cisbp), 10)

  expect_is(RBNS, "data.frame")
  expect_is(RBNS$motif, "character")
  expect_is(RBNS$ctrl, "integer")
  expect_is(RBNS$case, "integer")
  expect_is(RBNS$ctrl_freq, "numeric")
  expect_is(RBNS$case_freq, "numeric")
  expect_is(RBNS$log2FC, "numeric")
  expect_is(RBNS$ctrl_tot, "integer")
  expect_is(RBNS$case_tot, "integer")
  expect_is(RBNS$pval, "numeric")
  expect_is(RBNS$p_adj, "numeric")
  expect_equal(nrow(RBNS), length(RBNS_PWM[1:10]))
  expect_equal(ncol(RBNS), 10)

})

test_that("Output produces expected values", {
  expect_equal(any(is.na(cisbp[,c(1:5,7:10)])), FALSE) #log2FC can be NAN if no counts made
  expect_equal(all(names(CISBPRNA_hs_PWM[1:10]) %in% cisbp$motif), TRUE)
  expect_equal(cisbp$ctrl[1:10], c(558, 1602, 310, 1140, 680, 714, 744, 2114, 818, 828))
  expect_equal(cisbp$case[1:10], c(466, 897, 228, 789, 476, 435, 458, 1360, 520, 550))
  expect_equal(cisbp$ctrl_freq[1:10], c(0.002819948, 0.008095979, 0.001566638, 0.005761184, 0.003436496, 0.00360832, 0.00375993, 0.01068346, 0.004133902, 0.004184439))
  expect_equal(cisbp$case_freq[1:10], c(0.00354834, 0.006830175, 0.001736098, 0.006007812, 0.003624485, 0.003312292, 0.003487425, 0.01035567, 0.003959522, 0.004187955))
  expect_equal(cisbp$log2FC[1:10], c(0.3314759, -0.2452831, 0.1481767, 0.0604745, 0.07683794, -0.1234976, -0.1085439, -0.04495761, -0.06217811, 0.001211966), tolerance = .0001)
  expect_equal(cisbp$ctrl_tot[1:10], c(8950, 7906, 9198, 8368, 8828, 8794, 8764, 7394, 8690, 8680))
  expect_equal(cisbp$case_tot[1:10], c(5713, 5282, 5951, 5390, 5703, 5744, 5721, 4819, 5659, 5629))
  expect_equal(cisbp$pval[1:10], c(4.031473e-05, 9.310465e-05, 0.1509896, 0.1490519, 0.1998503, 0.2725896, 0.3568087, 0.7529026, 0.7037642, 0.6861791), tolerance = .0001)
  expect_equal(cisbp$p_adj[1:10], c(0.0004031473, 0.0004655232, 0.3774741, 0.3774741, 0.3997005, 0.4543161, 0.5097267, 0.7529026, 0.7529026, 0.7529026), tolerance = .0001)

  expect_equal(any(is.na(RBNS[,c(1:5,7:10)])), FALSE) #log2FC can be NAN if no counts made
  expect_equal(all(names(RBNS_PWM[1:10]) %in% RBNS$motif), TRUE)
  expect_equal(RBNS$ctrl[1:10], c(4156, 831, 2077, 2621, 5344, 1749, 3319, 3188, 3865, 1180))
  expect_equal(RBNS$case[1:10], c(2534, 578, 1264, 1720, 3341, 1139, 2133, 1986, 2465, 754))
  expect_equal(RBNS$ctrl_freq[1:10], c(0.02100305, 0.0041996, 0.01049647, 0.01324567, 0.02700681, 0.008838869, 0.01677313, 0.0161111, 0.01953243, 0.005963331))
  expect_equal(RBNS$case_freq[1:10], c(0.01929505, 0.00440116, 0.009624683, 0.01309688, 0.02543993, 0.008672875, 0.01624165, 0.01512233, 0.01876965, 0.005741306))
  expect_equal(RBNS$log2FC[1:10], c(-0.122368, 0.06763213, -0.1250936, -0.01629767, -0.08622891, -0.02735143, -0.04645355, -0.09137489, -0.05746965, -0.05473932), tolerance = .0001)
  expect_equal(RBNS$ctrl_tot[1:10], c(24174, 27499, 26253, 25709, 22986, 26581, 25011, 25142, 24465, 27150))
  expect_equal(RBNS$case_tot[1:10], c(15380, 17336, 16650, 16194, 14573, 16775, 15781, 15928, 15449, 17160))
  expect_equal(RBNS$pval[1:10], c(0.1186971, 0.07553338, 0.2687028, 0.2135907, 0.5740209, 0.4301321, 0.5342328, 0.5857138, 0.7284855, 0.8300984), tolerance = .0001)
  expect_equal(RBNS$p_adj[1:10], c(0.5934853, 0.5934853, 0.6717569, 0.6717569, 0.7321423, 0.7321423, 0.7321423, 0.7321423, 0.8094284, 0.8300984), tolerance = .0001)

  })

PWMnonList <- RBNS_PWM[[1]]
long_PWM <- list(one = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("ACACA", "CACAC", "AAACA"))),
                 two = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("TTTAA", "TTTTAA", "ATTAA"))),
                 three = Biostrings::consensusMatrix(Biostrings::DNAStringSet(c("GCATA", "CCATAA", "GGGAT"))))

test_that("motif_compare detects incorrect inputs", {
  expect_warning(motif_compare(RBNS_PWM, case, case), "some sequences in case set are also in the control set. This is not recommended.")
  expect_error(motif_compare(PWMnonList, case, ctrl), "PWM_list must be a list")
  expect_error(motif_compare(long_PWM, case, ctrl), "Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")", fixed = TRUE)

})
