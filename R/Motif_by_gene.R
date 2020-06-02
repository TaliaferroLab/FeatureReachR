#' Count Motif occurances within each sequence of a DNAStringSet
#'
#' This function simply returns the number of times a PWM matches with each
#' sequence in a DNAStringSet object.
#'
#' @param PWM_list a named list of position weight matrices. A list of named
#'   matrices where the name is unique, and the matrix columns represent
#'   positions, and the rows are the probability of each base "A", "C", "G", and
#'   "T". Each column should sum to 1. RNAreachR has three PWM_lists build in:
#'   "CISBPRNA_mm_PWM", "CISBPRNA_hs_PWM", and "RBNS_PWM".
#' @param DNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return a dataframe containing motif counts across each sequences within the
#'   DNAStringSet
#' @seealso \code{\link{getTxOut}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{countPWM}}, \code{\link{motif_compare}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#'
#' custom_PWM_by_gene(CISBPRNA_mm_PWM, case_fasta)
#' custom_PWM_by_gene(RBNS_PWM, case_fasta)
#'
#' @export
motif_by_gene <- function(PWM_list, DNAStringset){

  print("Checking PWM data", quote = FALSE)

  if (typeof(PWM_list) != "list"){
    stop("PWM_list must be a list")
  }

  if (any(as.numeric(lapply(PWM_list, nrow)) != 4)){
    stop("Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")")
  }
  motifs <- names(PWM_list) %>% dplyr::as_tibble() %>% dplyr::rename("motif" = value) %>% dplyr::mutate(PWM = unname(PWM_list))

  print("Counting motif occurances...", quote = FALSE)

  counts_list <- lapply(motifs$PWM, function(x) lapply(DNAStringset, function(y) Biostrings::countPWM(x, y)) %>% unlist() %>% dplyr::as_tibble(rownames = "gene"))
  names(counts_list) <- motifs$motif

  print("Counting motif occurance complete.", quote = FALSE)

  motif_by_gene <- dplyr::bind_rows(counts_list, .id = "motif") %>% tidyr::spread(gene, value = value)

  return(motif_by_gene)
}

