context("ggplot wrappers")
library(FeatureReachR)

case <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Enhanced.fasta", package = "FeatureReachR"))
ctrl <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Control.fasta", package = "FeatureReachR"))


kmer_stats <- kmer_compare(case,ctrl,4)
kmer_stats_ns <- kmer_compare(case[1:10], ctrl,4)


data(RBFOX2_RBNS_compare, package="FeatureReachR")
motif_stats <- RBFOX2_RBNS_compare
motif_stats_ns <- motif_compare(RBNS_PWM[1:5], case, ctrl)

test_that("Plot returns ggplot object",{
  expect_is(kmer_plot(kmer_stats), "ggplot")
  expect_is(motif_plot(motif_stats), "ggplot")
  expect_is(kmer_plot(kmer_stats_ns), "ggplot")
  expect_is(motif_plot(motif_stats_ns), "ggplot")
  expect_is(length_plot(case, ctrl), "ggplot")
  expect_is(GC_plot(case, ctrl), "ggplot")
})

