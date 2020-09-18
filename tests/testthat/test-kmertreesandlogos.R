context("kmer trees and logos")
library(RNAreachR)

case <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Enhanced.fasta", package = "RNAreachR"))
ctrl <- Biostrings::readDNAStringSet(system.file("extdata", "DownstreamIntron.Control.fasta", package = "RNAreachR"))

enriched_sixmers <- c("AAGGAA", "ACACAC", "AGAAGG", "AGAGAG", "AGAGGG", "AGGAAG", "AGGAGG", "AGGGAG", "CACACA", "GAAGGA", "GAGAAG", "GAGAGA", "GAGGAG", "GAGGGA", "GAGGGG", "GGAAGG", "GGAGGA", "GGAGGG", "GGGAGG")
clustered_sixmers <- c("ACACAC", "ACACCC", "CCACAC")
alldiff_sixmers <- c("GGGGGG", "TTTTTT", "AAAAAA")


test_that("kmer grouping fx returns correct typs", {
  expect_is(kmer2logo(enriched_sixmers, "My Plot Title"), "ggplot")
  expect_is(kmer2tree(enriched_sixmers, "My Tree Title"), "NULL")
  expect_is(kmer2PWM(enriched_sixmers), "list")
  expect_equal(length(kmer2PWM(enriched_sixmers)), 6)
})

test_that("kmer grouping fx detects incorrect inputs", {
  expect_output(kmer2logo(clustered_sixmers), "All kmers are within a Levenshtein distance of 2")
  expect_output(kmer2logo(alldiff_sixmers), "No kmers are within a Levenshtein distance of 2, each kmer is its own cluster")
  expect_warning(kmer2tree(clustered_sixmers), "all kmers are within the same custer and cannot be plotted as a tree")
  expect_warning(kmer2tree(alldiff_sixmers), "each kmer is in its own cluster and cannot be plotted as a tree")
  expect_output(kmer2PWM(clustered_sixmers), "All kmers are within a Levenshtein distance of 2")
  expect_output(kmer2PWM(alldiff_sixmers), "No kmers are within a Levenshtein distance of 2, each kmer is its own cluster")

})
