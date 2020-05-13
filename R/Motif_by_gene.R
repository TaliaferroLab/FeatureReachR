#' Count Motif occurances within each sequence of a DNAStringSet
#'
#' This function simply returns every motif counts for each gene. Begin with
#' folder containing PWM for motifs of intrest (downloaded from CisBPRNA called
#' "pwms_all_motifs") also need RBP info for the motifs (dowloaded from CisBPRNA
#' called "RBP_Information_all_motifs")
#' cisbpRNA:
#' \url{http://cisbp-rna.ccbr.utoronto.ca/bulk.php} ensure RBP info and PFMs are
#' selected for downloading an entire species archive
#'
#' @param motif_path A path. To the directory containing Motif PWM txt files.
#' @param RBPinfo A path. To the RBPinfo table downloaded from cisbpRNA
#' @param DNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param PWM_list a named list of matrices. A list of named matrices where the
#'   name is unique, and the matrix columns represent positions, and the rows
#'   are the probability of each base "A", "C", "G", and "T". Each column should
#'   sum to 1.
#' @return a dataframe containing motif counts across each sequences within the
#'   DNAStringSet
#' @seealso \code{\link{getTxOut}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{countPWM}}, \code{\link{cisBPRNA_compare}},
#'    \code{\link{RBNS_compare}}, \code{\link{custom_PWM_compare}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#'
#' cisBPRNA_by_gene("mydata/CisBPRNAdat/mm/pwms_all_motifs",
#' "mydata/CisBPRNAdat/mm/RBP_Information_all_motifs.txt", case_fasta)
#'
#' RBNS_by_gene("mydata/RBNSdat/RBNS_PWMs", case_fasta)
#'
#' custom_PWM_list <- readRDS("mydata/example_custom_PWM_list")
#' custom_PWM_by_gene(custom_PWM_list, case_fasta)
#'
#' @export
#' @describeIn Counts motif occurances for each sequence in a DNAStringSet.
#' cisbpRNA is an online library of RNA binding proteins and their motifs.
#' Motifs are derived from both direct and inferred binding.
cisBPRNA_by_gene <- function(motif_path, RBPinfo, DNAStringset){

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

  print("Getting PWM data", quote = FALSE)

  motifs <- motifs %>%
    dplyr::mutate(PWM = lapply(PATH, function(x) t(read.table(x, header = TRUE, row.names = 1, col.names = c("pos", "A", "C", "G", "T")))))

  print("Counting motif occurances...", quote = FALSE)

  counts_list <- lapply(motifs$PWM, function(x) lapply(DNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% dplyr::as_tibble(rownames = "gene"))
  names(counts_list) <- motifs$motif

  print("Counting motif occurances complete.", quote = FALSE)

  motif_by_gene <- dplyr::bind_rows(counts_list, .id = "motif") %>%
    tidyr::spread(gene, value = value) %>%
    dplyr::left_join(., RBP_info, by = c("motif" = "Motif_ID")) %>%
    dplyr::select(motif, RBP_name, dplyr::everything())

  return(motif_by_gene)
}
#' @describeIn Counts motif occurances for each sequence in a DNAStringSet.
#' RNA bind and seq is an In vitro technique to define RBP binding sites.
#' These assays are performed wiht purified RBP and random synthetic RNAs.
RBNS_by_gene <- function(motif_path, DNAStringset){

  #get paths to each motif pwm

  print("getting Motif and RBP data...", quote = FALSE)

  motif_paths <- list.files(path = motif_path, full.names = TRUE)
  motif_info <- file.info(path = motif_paths)
  motifs <- motif_info %>%
    dplyr::as_tibble(rownames = "PATH") %>%
    dplyr::mutate(motif = stringr::str_match(PATH, "PWMs/(.*?).PWM")[,2]) %>%
    dplyr::select(PATH, motif)

  print("Getting PWM data", quote = FALSE)

  motifs <- motifs %>%
    dplyr::mutate(PWM = lapply(PATH, function(x) t(read.table(x, skip = 1, row.names = 1, header = TRUE, col.names = c("pos","A", "C", "G", "T")))))

  print("Counting motif occurances...", quote = FALSE)

  counts_list <- lapply(motifs$PWM, function(x) lapply(DNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% dplyr::as_tibble(rownames = "gene"))
  names(counts_list) <- motifs$motif

  print("Counting motif occurance complete.", quote = FALSE)

  motif_by_gene <- dplyr::bind_rows(counts_list, .id = "motif") %>% tidyr::spread(gene, value = value)

  return(motif_by_gene)
}
#' @describeIn Counts custom PWM motif occurances for each sequence in a DNAStringSet.
#' This is useful for creating shorter PWM lists of special interest which will
#' count motifs much faster.
#'
custom_PWM_by_gene <- function(PWM_list, DNAStringset){

  print("Checking PWM data", quote = FALSE)

  if (any(as.numeric(lapply(PWM_list, nrow)) != 4)){
    stop("Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")")
  }
  motifs <- names(PWM_list) %>% as_tibble() %>% rename(value = "motif") %>% mutate(PWM = unname(PWM_list))

  print("Counting motif occurances...", quote = FALSE)

  counts_list <- lapply(motifs$PWM, function(x) lapply(DNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% dplyr::as_tibble(rownames = "gene"))
  names(counts_list) <- motifs$motif

  print("Counting motif occurance complete.", quote = FALSE)

  motif_by_gene <- dplyr::bind_rows(counts_list, .id = "motif") %>% tidyr::spread(gene, value = value)

  return(motif_by_gene)
}

