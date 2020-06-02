#' Compare Motif occurances between two DNAStringSets
#'
#' This function compares motif occurances between between case and control
#' DNAStringSets.
#'
#' @param PWM_list a named list of position weight matrices. A list of named
#'   matrices where the name is unique, and the matrix columns represent
#'   positions, and the rows are the probability of each base "A", "C", "G", and
#'   "T". Each column should sum to 1. RNAreachR has three PWM_lists build in:
#'   "CISBPRNA_mm_PWM", "CISBPRNA_hs_PWM", and "RBNS_PWM".
#' @param caseDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param ctrlDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return a summary table including the pvalue as calculated by
#'   \code{wilcox.test}, then corrected by \code{p.adjust} with the \code{method = "BH"}.
#'   The mean values for each DNAStringSet, and the log2 fold
#'   change between those means.
#' @seealso \code{\link{getTxOut}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{countPWM}}, \code{\link{motif_by_gene}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
#'
#' motif_compare(CISBPRNA_mm_PWM, case_fasta, ctrl_fasta)
#' motif_compare(RBNS_PWM, case_fasta, ctrl_fasta)
#'
#' @export
motif_compare <- function(PWM_list, caseDNAStringset, ctrlDNAStringSet){

  print("Checking PWM data", quote = FALSE)

  if (typeof(PWM_list) != "list"){
    stop("PWM_list must be a list")
  }

  if (any(as.numeric(lapply(PWM_list, nrow)) != 4)){
    stop("Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")")
  }
  motifs <- names(PWM_list) %>% dplyr::as_tibble() %>% dplyr::rename("motif" = value) %>% dplyr::mutate(PWM = unname(PWM_list))

  fisher <- function(a, b, c, d){
    mat <- matrix(c(a, b, c, d), nr = 2)
    fisher.test(mat, alternative = "two.sided")$p.value
  }

  print("Counting motif occurances and calculating statistics...", quote = FALSE)


  motifs <- motifs %>%
    mutate(case = suppressWarnings(unlist(lapply(PWM, function(x) lapply(caseDNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% sum()))),
           ctrl = suppressWarnings(unlist(lapply(PWM, function(x) lapply(ctrlDNAStringSet, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% sum()))),
           case_freq = case / sum(Biostrings::width(caseDNAStringset)),
           ctrl_freq = ctrl / sum(Biostrings::width(ctrlDNAStringSet)),
           log2FC = log2(case_freq/ctrl_freq),
           case_tot = sum(case)-case,
           ctrl_tot = sum(ctrl)-ctrl) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pval = fisher(case, ctrl, case_tot, ctrl_tot),
                  p_adj = p.adjust(pval, method = "BH", nrow(motifs))) %>%
    dplyr::select(motif, case, ctrl, case_freq, ctrl_freq, log2FC, case_tot, ctrl_tot, pval, p_adj) %>%
    dplyr::arrange(p_adj)

  return(motifs)

  print("All Finished.", quote = FALSE)
}
