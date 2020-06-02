#' Get length or GC of each sequence in a DNAStringSet object
#'
#' These functions return data frames containing the length or GC content for
#' each sequence in a DNAStringSet object.
#'
#' @param DNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return A data frame with the length or GC content for each sequence in a
#'   DNAStringSet object
#' @seealso \code{\link{getTxOut}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{letterFrequency}}
#' @export
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' get_length(case_fasta)
#' get_GC(case_fasta)
#' @describeIn retrieves the length of each sequence in a
#'   DNAStringSet object and returns a dataframe.
get_length <- function(DNAStringSet){
  dplyr::as_tibble(data.frame(gene = names(DNAStringSet),
                              length = Biostrings::width(DNAStringSet)))
}
#' @export
#' @describeIn retrieves both the length and GC content of each
#'   sequence in a DNAStringSet object and returns a dataframe.
get_GC <- function(DNAStringSet){
  GC <- data.frame(gene = names(DNAStringSet),
                   GC = Biostrings::letterFrequency(DNAStringSet, "GC") / Biostrings::width(DNAStringSet)) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(GC = G.C)
  return(GC)
}
