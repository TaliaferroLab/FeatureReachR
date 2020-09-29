#' Compare kmer content between two DNAStringSets
#'
#' This function compares kmer content between between case and control
#' DNAStringSets
#'
#' @param caseDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param  ctrlDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return a summary table including the pvalue as calculated by
#'   \code{wilcox.test}, then corrected by \code{p.adjust} with the \code{method = "BH"}.
#'   The mean values for each DNAStringSet, and the log2 Fold
#'   change between those means
#' @seealso \code{\link{write_Sequence}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{oligonucleotideFrequency}},\code{\link{kmer_by_gene}},
#'    \code{\link[effsize]{cliff.delta}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
#' kmer_compare(case_fasta, ctrl_fasta, 6)
#' @export
#'
kmer_compare <- function(caseDNAStringSet, ctrlDNAStringSet, k){

  if (any(names(caseDNAStringSet) %in% names(ctrlDNAStringSet))){
    warning("some sequences in case set are also in the control set. This is not recommended.")
  }

  print("counting kmers...", quote = FALSE)
  case_kmer <- Biostrings::oligonucleotideFrequency(caseDNAStringSet, width = k) %>%
    colSums() %>%
    data.frame(kmer = names(.), case = .) %>%
    dplyr::as_tibble()
  ctrl_kmer <- Biostrings::oligonucleotideFrequency(ctrlDNAStringSet, width = k) %>%
    colSums() %>%
    data.frame(kmer = names(.), ctrl = .) %>%
    dplyr::as_tibble()
  print("counting complete.", quote= FALSE)

  #compare kmers between case and ctrl takes

  fisher <- function(a, b, c, d){
    mat <- matrix(c(a, b, c, d), nr = 2)
    fisher.test(mat, alternative = "two.sided")$p.value
  }

  print("calculating kmer statistics...", print = FALSE)

  kmer_stats <- dplyr::left_join(ctrl_kmer, case_kmer) %>%
    na.omit() %>%
    dplyr::mutate(ctrl_freq = ctrl / sum(ctrl),
                  case_freq = case / sum(case),
                  log2FC = log2(case_freq/ctrl_freq),
                  ctrl_tot = sum(ctrl)-ctrl,
                  case_tot = sum(case)-case) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pval = fisher(case, ctrl, case_tot, ctrl_tot)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p_adj = p.adjust(pval, method = "BH")) %>%
    dplyr::arrange(p_adj)

  print("calculations complete.", quote = FALSE)

  return(kmer_stats)
}
