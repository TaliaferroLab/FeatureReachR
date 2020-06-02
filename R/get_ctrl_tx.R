#' Get Control Transcript list from Case gene list
#'
#' This will return transcripts representing all other genes not in
#' \code{gene_list} This requires creating a \code{length_df} see
#' \code{make_longest} or \code{make_median}. Generating a control transcript
#' list with this function isn't recommended as not all genes will be exprsesed
#' in a given system.It is best to use expression data from your system to
#' create a list of expressed control transcripts.
#'
#' @param length_df A dataframe relating genes and transcritps by their length.
#'   Should be the output from \code{make_longest} or \code{make_median}.
#' @param gene_list A character list of ensembl gene IDs.
#' @param seq_type One of "whole", "CDS", "5UTR" or "3UTR".
#' @return a character list of ensembl  transcript IDs.
#' @examples
#' mm_case_gene <- c("ENSMUSG00000097392", "ENSMUSG00000025607", "ENSMUSG00000030671", "ENSMUSG00000034764", "ENSMUSG00000116215", "ENSMUSG00000039556", "ENSMUSG00000066510", "ENSMUSG00000018160", "ENSMUSG00000114306", "ENSMUSG00000028277", "ENSMUSG00000037216", "ENSMUSG00000032299") #12genes
#' get_ctrl_tx(longest_mm, mm_case_gene, "UTR3")
#' @export

get_ctrl_tx <- function(length_df, gene_list, seq_type){
  #Check that gene list contains genes

  if (all(grepl("ENST", gene_list)) | all(grepl("ENSMUST", gene_list))){
    stop("Please ensure the gene list contains gene IDs and not transcript IDs.")
  }
  #Make sure df and gene list are the same species

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

  if (list_species != gff_species){
    stop("Please ensure the gene list and gff are of the same species.")
  }

  #get longest transcript list of alll other genes not in gene list
  `%out%` = Negate(`%in%`)
  ctrl_tx_list <- length_df %>% dplyr::filter(gene_id %out% gene_list) %>% dplyr::pull(., seq_type)
  return(ctrl_tx_list)

}
