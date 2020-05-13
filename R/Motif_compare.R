#' Compare Motif occurances between two DNAStringSets
#'
#' This function compares motif occurances between between case and control
#' DNAStringSets. Begin with a folder containing PWM for motifs of intrest
#' (downloaded from CisBPRNA called "pwms_all_motifs") also need RBP info for
#' the motifs (dowloaded from CisBPRNA called "RBP_Information_all_motifs")
#' cisbpRNA: \url{http://cisbp-rna.ccbr.utoronto.ca/bulk.php} ensure RBP info
#' and PFMs are selected for downloading an entire species archive
#'
#' @param motif_path A path. To the directory containing Motif PWM txt files.
#' @param RBPinfo A path. To the RBPinfo table downloaded from cisbpRNA
#' @param caseDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param ctrlDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param PWM_list a named list of matrices. A list of named matrices where the
#'   name is unique, and the matrix columns represent positions, and the rows
#'   are the probability of each base "A", "C", "G", and "T". Each column should
#'   sum to 1.
#' @return a summary table including the pvalue as calculated by
#'   \code{wilcox.test}, then corrected by \code{p.adjust} with the \code{method = "BH"}.
#'   The mean values for each DNAStringSet, and the log2 Fold
#'   change between those means
#' @seealso \code{\link{getTxOut}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{countPWM}}, \code{\link{cisBPRNA_by_gene}},
#'   \code{\link{RBNS_by_gene}}, \code{\link{custom_PWM_by_gene}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
#'
#' cisBPRNA_compare("mydata/CisBPRNAdat/mm/pwms_all_motifs",
#' "mydata/CisBPRNAdat/mm/RBP_Information_all_motifs.txt", case_fasta, ctrl_fasta)
#'
#' RBNS_commpare("mydata/RBNSdat/RBNS_PWMs", case_fasta, ctrl_fasta)
#'
#' custom_PWM_list <- readRDS("mydata/example_custom_PWM_list")
#' custom_PWM_compare(custom_PWM_list, case_fasta, ctrl_fasta)
#'
#' @export
#' @describeIn Counts motif occurances in case and ctrl DNAStringSets then
#' produces a summary table.
#' cisbpRNA is an online library of RNA binding proteins and their motifs.
#' Motifs are derived from both direct and inferred binding.
cisBPRNA_compare <- function(motif_path, RBPinfo, caseDNAStringset, ctrlDNAStringSet){

  #get paths to each motif pwm

  print("getting Motif and RBP data...", quote = FALSE)

  motif_paths <- list.files(path = motif_path, full.names = TRUE)
  motif_info <- file.info(motif_paths)
  motif_info <- motif_info[motif_info$size != 0, ]
  motifs <- motif_info %>%
    dplyr::as_tibble(rownames = "PATH") %>%
    dplyr::mutate(motif = stringr::str_match(PATH, "motifs/(.*?).txt")[,2]) %>%
    dplyr::select(PATH, motif)

  #merge motif paths with RBP info

  RBP_info <- read.table(RBPinfo, header = TRUE, sep = "\t")
  RBP_info <- RBP_info %>%
    dplyr::as_tibble() %>%
    dplyr::select(Motif_ID, RBP_Name) %>%
    dplyr::filter(Motif_ID != ".") %>%
    dplyr::group_by(Motif_ID) %>%
    dplyr::summarise(RBP_name = dplyr::first(RBP_Name))

  motifs <- dplyr::left_join(motifs, RBP_info, by = c("motif" = "Motif_ID"))

  fisher <- function(a, b, c, d){
    mat <- matrix(c(a, b, c, d), nr = 2)
    fisher.test(mat, alternative = "two.sided")$p.value
  }

  print("Counting motif occurances and calculating statistics...", quote = FALSE)

  motifs <- motifs %>%
    mutate(PWM = lapply(PATH, function(x) t(read.table(x, header = TRUE, row.names = 1, col.names = c("pos", "A", "C", "G", "T")))),
           case = suppressWarnings(unlist(lapply(PWM, function(x) lapply(caseDNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% sum()))),
           ctrl = suppressWarnings(unlist(lapply(PWM, function(x) lapply(ctrlDNAStringSet, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% sum()))),
           case_freq = case / sum(Biostrings::width(caseDNAStringset)),
           ctrl_freq = ctrl / sum(Biostrings::width(ctrlDNAStringSet)),
           log2FC = log2(case_freq/ctrl_freq),
           case_tot = sum(case)-case,
           ctrl_tot = sum(ctrl)-ctrl) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pval = fisher(case, ctrl, case_tot, ctrl_tot),
                  p_adj = p.adjust(pval, method = "BH", nrow(motifs))) %>%
    dplyr::select(RBP_name, motif, case, ctrl, case_freq, ctrl_freq, log2FC, case_tot, ctrl_tot, pval, p_adj)

  return(motifs)
}
#' @describeIn Counts motif occurances in case and ctrl DNAStringSets then
#' produces a summary table.
#' RNA bind and seq is an In vitro technique to define RBP binding sites.
#' These assays are performed wiht purified RBP and random synthetic RNAs.
RBNS_compare <- function(motif_path, caseDNAStringset, ctrlDNAStringSet){

  #get paths to each motif pwm

  print("getting Motif and RBP data...", quote = FALSE)

  motif_paths <- list.files(path = motif_path, full.names = TRUE)
  motif_info <- file.info(path = motif_paths)
  motifs <- motif_info %>%
    dplyr::as_tibble(rownames = "PATH") %>%
    dplyr::mutate(motif = stringr::str_match(PATH, "PWMs/(.*?).PWM")[,2]) %>%
    dplyr::select(PATH, motif)

  fisher <- function(a, b, c, d){
    mat <- matrix(c(a, b, c, d), nr = 2)
    fisher.test(mat, alternative = "two.sided")$p.value
  }

  print("Counting motif occurances and calculating statistics...", quote = FALSE)


  motifs <- motifs %>%
    mutate(PWM = lapply(PATH, function(x) t(read.table(x, skip = 1, row.names = 1, header = TRUE, col.names = c("pos","A", "C", "G", "T")))),
           case = suppressWarnings(unlist(lapply(PWM, function(x) lapply(caseDNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% sum()))),
           ctrl = suppressWarnings(unlist(lapply(PWM, function(x) lapply(ctrlDNAStringSet, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% sum()))),
           case_freq = case / sum(Biostrings::width(caseDNAStringset)),
           ctrl_freq = ctrl / sum(Biostrings::width(ctrlDNAStringSet)),
           log2FC = log2(case_freq/ctrl_freq),
           case_tot = sum(case)-case,
           ctrl_tot = sum(ctrl)-ctrl) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pval = fisher(case, ctrl, case_tot, ctrl_tot),
                  p_adj = p.adjust(pval, method = "BH", nrow(motifs))) %>%
    dplyr::select(motif, case, ctrl, case_freq, ctrl_freq, log2FC, case_tot, ctrl_tot, pval, p_adj)

  return(motifs)

}
#' @describeIn Counts custom PWM motif occurances in case and ctrl DNAStringSets then
#' produces a summary table.
#' This is useful for creating shorter PWM lists of special interest which will
#' count motifs much faster.
#'
custom_PWM_compare <- function(PWM_list, caseDNAStringset, ctrlDNAStringSet){

  print("Checking PWM data", quote = FALSE)

  if (any(as.numeric(lapply(PWM_list, nrow)) != 4)){
    stop("Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")")
  }
  motifs <- names(PWM_list) %>% as_tibble() %>% rename(value = "motif") %>% mutate(PWM = unname(PWM_list))

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

}

