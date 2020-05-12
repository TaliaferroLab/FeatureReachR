#' Filter out undesired transcripts from gencode gff files
#'
#' \code{filter_Tx} returns a filtered TxDb object
#'
#' @param gff A path. To the human or mouse gencode annotation gff file. The
#'   most current annotations can be downloaded from:
#'   \url{https://www.gencodegenes.org/human/} or
#'   \url{https://www.gencodegenes.org/mouse/} There are two filter options, a
#'   general filter and a protein coding filter.
#' @param filter A logical scalar. Should transcripts with undefined ends and
#'   support levels higher than 2 (low confidence transcripts) be removed? The
#'   default is to remove low confidence transcripts (\code{filter == TRUE})
#' @param protein.coding A logical scalar. Should all transcripts that lack a
#'   protein coding tag be removed? The default is to keep all transcripts
#'   regardless of protein coding status (\code{protein.coding == FALSE})
#' @return The output from \code{filter_Tx()} is a TxDb object containing a gene
#'   model of the filtered transcripts This Txdb object is required for
#'   downstream functions. See
#'   \url{https://www.rdocumentation.org/packages/GenomicFeatures/versions/1.24.4/topics/TxDb-class}
#'   for more information
#' @examples
#' # No filtering
#' hs_TxDb <- filter_Tx("mydata/Gencodedat/gencode.v33.annotation.gff3.gz", filter = FALSE, protein.coding = FALSE)
#' mm_TxDb <- filter_Tx("mydata/Gencodedat/gencode.vM20.annotation.gff3.gz", filter = FALSE, protein.coding = FALSE)
#'
#' # Filter out low confidence transcripts
#' hs_filtered_TxDb <- filter_Tx("mydata/Gencodedat/gencode.v33.annotation.gff3.gz")
#' mm_filtered_TxDb <- filter_Tx("mydata/Gencodedat/gencode.vM20.annotation.gff3.gz")
#'
#' # Filter out low confidence transcripts and non-coding transcripts
#' hs_coding_TxDb <- filter_Tx("mydata/Gencodedat/gencode.v33.annotation.gff3.gz", protein.coding = TRUE)
#' mm_coding_TxDb <- filter_Tx("mydata/Gencodedat/gencode.vM20.annotation.gff3.gz", protein.coding = TRUE)
#' @export
filter_Tx <- function(gff, filter = TRUE, protein.coding = FALSE) {
  # import gencode gff as GRanges object
  gff <- rtracklayer::import(gff)

  # filter out low confidence transcripts
  if (filter == TRUE) {
    print("filtering out low confidence transcripts...", quote = FALSE)
    # remove low confidence transcripts
    gff <- gff[gff$transcript_support_level == "1" | gff$transcript_support_level == "2" | is.na(gff$transcript_support_level)]
    # remove transcripts with undefined starts or ends
    NF <- c("cds_start_NF", "cds_end_NF", "mRNA_start_NF", "mRNA_end_NF")
    `%out%` = Negate(`%in%`)
    gff <- gff[as.list(gff$tag) %out% NF]
    print("filtering complete.", quote = FALSE)
  }

  # filter out non-coding transcripts
  if (protein.coding == TRUE) {
    print("filtering out non-coding transcripts...", quote = FALSE)
    gff <- gff[gff$transcript_type == "protein_coding" | is.na(gff$transcript_type)]
    print("filtering complete.", quote = FALSE)
  }

  #create txdb gene model from filtered gff
  if (all(grepl("ENSG", gff$gene_id))) {
    txdbGFF <- GenomicFeatures::makeTxDbFromGRanges(gff)
    GenomeInfoDb::seqlevelsStyle(txdbGFF) <- "NCBI"
  } else if (all(grepl("ENSMUSG", gff$gene_id))) {
    txdbGFF <- GenomicFeatures::makeTxDbFromGRanges(gff)
  }

  return(txdbGFF)

}


