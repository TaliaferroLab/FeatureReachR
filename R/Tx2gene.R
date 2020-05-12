#' Get transcript list from gene list
#'
#' \code{gene2Tx} takes a gene list and returns a list of transcripts with the
#' longest or median length feature of interest.
#'
#' @param length_df A data frame. Created by \code{\link{make_longest_df}} or
#'   \code{\link{make_median_df}}.
#' @param gene_list A character list. Must contain Ensembl gene IDs and be
#'   compatible with the gff (the same species).
#' @param seq_type One of "whole", "CDS", "5pUTR", or "3pUTR".
#' @seealso \code{\link{filter_Tx}}, \code{\link{make_longest_df}},
#'   \code{\link{make_median_df}}
#' @return \code{gene2Tx} returns a character list of Ensembl transcript IDs
#'   associated with the length and sequence type of interest.
#' @importFrom magrittr "%>%"
#' @examples
#' hs_filtered_TxDb <- filter_Tx("mydata/Gencodedat/gencode.v33.annotation.gff3.gz")
#' longest_hs <- make_longest_df(hs_filtered_TxDb)
#' hs_gene <- c("ENSG00000000419.12", "ENSG00000001167.14", "ENSG00000000938.13")
#' hs_tx <- gene2Tx(longest_hs, hs_gene, "UTR3")
#' @export
gene2Tx <- function(length_df, gene_list, seq_type) {
  #Check that gene list contains genes

  if (all(grepl("ENST", gene_list)) | all(grepl("ENSMUST", gene_list))) {
    stop("Please ensure the gene list contains gene IDs and not transcript IDs.")
  }

  #Make sure gff and gene list are the same species

  ifelse (all(grepl("ENSG", gene_list)),
          (list_species = "human"),
          ifelse (all(grepl("ENSMUSG", gene_list)),
                  (list_species = "mouse"),
                  (list_species = "unknown")))

  ifelse (grepl("ENSG", length_df$gene_id[1]),
          (gff_species = "human"),
          ifelse (grepl("ENSMUSG", length_df$gene_id[1]),
                  (gff_species = "mouse"),
                  (gff_species = "unknown")))

  if (list_species != gff_species) {
    stop("Please ensure the gene list and gff are of the same species.")
  }

  #create longest transcript list from gene list

  tx_list <- length_df %>%
    dplyr::filter(gene_id %in% gene_list) %>%
    dplyr::pull(., seq_type)
  return(tx_list)

}

