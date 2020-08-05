context("Compare Outputs")
library(RNAreachR)

case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
x <- length_compare(case_fasta, ctrl_fasta)

expect_is(x, data.frame)
expect_is(x$p.value, numeric)
expect_is(x$)
