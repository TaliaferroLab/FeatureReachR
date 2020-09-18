#' Relate Kmers and PWM list
#'
#' Create a logical matrix relating position weight matrices to kmers.
#' This creates the input for \code{estimate_motif_from_kmer}
#' These matrices take time to make. Logical matrices relating RBNS
#' and CISBPRNA motifs to 4-, 5-, 6- and 7-mers are included within
#' RNAreachR.
#'
#' @param k An integer. Recommended values of K are 3 < k < 8.
#' @param PWM_list a named list of position weight matrices. A list of named
#'  matrices where the name is unique, and the matrix columns represent
#'  positions, and the rows are the probability of each base "A", "C", "G", and
#'  "T". Each column should sum to 1. RNAreachR has three PWM_lists build in:
#'  "CISBPRNA_mm_PWM", "CISBPRNA_hs_PWM", and "RBNS_PWM".
#' @return A logial matrix indicating if a motif can match to each kmer.
#' @seealso \code{\link[Biostrings]{countPWM}}, \code{\link{motif_compare}}
#' @examples
#' #this table is already saved in RNAreachR
#' relate_Kmer_PWM(4, RBNS_PWM)
#' relate_kmer_PWM(6, custom_PWM_list)
#' @export
relate_kmer_PWM <- function(k, PWM_list){

  ##make DNAStringSets for every kmer
  x <- c("T", "A", "C", "G")

  kmers <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer) %>%
    Biostrings::DNAStringSet()

  names(kmers) <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer)

  #Check PWM_list
  print("Checking PWM data", quote = FALSE)

  if (typeof(PWM_list) != "list"){
    stop("PWM_list must be a list")
  }

  if (any(as.numeric(lapply(PWM_list, nrow)) != 4)){
    stop("Error: Ensure all PWM matricies have exactly 4 rows (A, C, G and T)")
  }

  #this function makes every PWM the same length as k or shorter by tiling across longer PWM matrices
  mattiles <- function(Mat, k) {

    if(any(lapply(Mat, ncol) > k)){
      x <- c(1:(ncol(Mat[[1]])-(k-1)))-1
      list <- lapply(x, function(x) Mat[[1]][,(1+x):(k+x)])
      names(list) <- lapply(x, function(x) paste(names(Mat), "_", x + 1, sep = ""))
      list
    }
    else
      Mat
  }

  x <- 1:length(PWM_list)
  tiled_motifs <- unlist(lapply(x, function(x) mattiles(PWM_list[x], k)), recursive = FALSE)


  print("Counting motif occurances...", quote = FALSE)

  counts_list <- lapply(tiled_motifs, function(x)
    lapply(kmers, function(y)
      Biostrings::countPWM(x, y)) %>% unlist() %>% dplyr::as_tibble(rownames = "kmer"))

  print("Counting motif occurances complete.", quote = FALSE)

  motif_by_kmer <- dplyr::bind_rows(counts_list, .id = "motif") %>%
    tidyr::spread(kmer, value = value) %>%
    dplyr::select(motif, dplyr::everything())

  motif_by_kmer <- motif_by_kmer %>% dplyr::mutate_if(is.numeric, ~1 * (. > 0))

  return(motif_by_kmer)
}
