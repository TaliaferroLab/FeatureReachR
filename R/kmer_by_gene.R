#' Count Kmers within each sequence of a DNAStringSet
#'
#' This function simply returns the counts of each kmer for each gene in a
#' DNAStringSet.
#'
#' @param DNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param k An Integer. The length of kmers to count. Values 3 < k < 8 are
#'   recommended.
#' @return a dataframe containing kmer counts across each sequences within the DNAStringSet
#' @seealso \code{\link{write_Sequence}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{oligonucleotideFrequency}}, \code{\link{kmer_compare}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' case_fasta_kmer_counts <- kmer_by_gene(case_fasta, 6)
#' @export

kmer_by_gene <- function(DNAStringSet, k){

  print("counting kmers...", quote = FALSE)
  kmer_counts <- Biostrings::oligonucleotideFrequency(DNAStringSet, width = k) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(gene = names(DNAStringSet)) %>%
    dplyr::select(gene, dplyr::everything())

  print("counting complete.", quote= FALSE)

  return(kmer_counts)
}
