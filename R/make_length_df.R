#' Create transcript length reference dataframe
#'
#' Because creating these dataframes takes a minute or two and its likely a user
#' will be interested in multiple seq types (whole transcript, CDS, 5'UTR,
#' 3'UTR) it is better to create and save this dataframe once.
#'
#' @param TxDb_gff A TxDb object. Using \code{\link{filter_Tx}} to create this
#'   TxDb object is recommended.
#' @return Both \code{make_longest_df} and \code{make_median_df} create a
#'   dataframe relating Ensembl gene IDs and ensembl Transcript IDs. Because the
#'   longest CDS and longest 3'UTR  for a gene may not belong to the same
#'   transcript, there is a column for each sequence type.
#' @seealso \code{\link{filter_Tx}}, \code{\link{gene2tx}}
#' @examples
#' hs_filtered_TxDb <- filter_Tx("mydata/Gencodedat/gencode.v33.annotation.gff3.gz")
#' longest_hs <- make_longest_df(hs_filtered_TxDb)
#' median_hs <- make_median_df(hs_filtered_TxDb)
#'
#' @describeIn make_length_df creates a dataframe of the
#'   transcripts with the longest features for each gene in a TxDb object
make_longest_df <- function(TxDb_gff){
  #identify transcript with longest feature for each gene

  len_df <- GenomicFeatures::transcriptLengths(TxDb_gff,
                                               with.cds_len = TRUE,
                                               with.utr5_len = TRUE,
                                               with.utr3_len = TRUE) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

  longest_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::top_n(1, tx_len) %>%
    dplyr::rename("whole" = "tx_name") %>%
    dplyr::select(gene_id, whole)

  longest_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::top_n(1, cds_len) %>%
    dplyr::rename("CDS" = "tx_name") %>%
    dplyr::select(gene_id, CDS) %>%
    dplyr::left_join(., longest_tx)

  longest_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::top_n(1, utr5_len) %>%
    dplyr::rename("UTR5" = "tx_name") %>%
    dplyr::select(gene_id, UTR5) %>%
    dplyr::left_join(., longest_tx)

  longest_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::top_n(1, utr3_len) %>%
    dplyr::rename("UTR3" = "tx_name") %>%
    dplyr::select(gene_id, UTR3) %>%
    dplyr::left_join(., longest_tx)

  longest_tx <- longest_tx %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(whole = dplyr::first(whole),
                     CDS = dplyr::first(CDS),
                     UTR5 = dplyr::first(UTR5),
                     UTR3 = dplyr::first(UTR3))

}

#' @describeIn make_length_df creates a dataframe of the
#'   transcripts with the median length features for each gene in a TxDb object
make_median_df <- function(TxDb_gff){
  #identify transcript with longest feature for each gene

  len_df <- GenomicFeatures::transcriptLengths(TxDb_gff,
                                               with.cds_len = TRUE,
                                               with.utr5_len = TRUE,
                                               with.utr3_len = TRUE) %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(gene_id = sub("\\..*", "", gene_id))

  #ignore isoforms with no feature when calculating median (many transcripts lack CDS, UTR3 and UTR5 features)
  len_df[len_df == 0] <- NA

  median_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(abs(tx_len - median(tx_len, na.rm = TRUE)) == min(abs(tx_len - median(tx_len, na.rm = TRUE)), na.rm = TRUE)) %>%
    dplyr::rename("whole" = tx_name) %>%
    dplyr::select(gene_id, whole)

  median_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(abs(cds_len - median(cds_len, na.rm = TRUE)) == min(abs(cds_len - median(cds_len, na.rm = TRUE)), na.rm = TRUE)) %>%
    dplyr::rename("CDS" = tx_name) %>%
    dplyr::select(gene_id, CDS) %>%
    dplyr::left_join(., median_tx)

  median_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(abs(utr5_len - median(utr5_len, na.rm = TRUE)) == min(abs(utr5_len - median(utr5_len, na.rm = TRUE)), na.rm = TRUE)) %>%
    dplyr::rename("UTR5" = tx_name) %>%
    dplyr::select(gene_id, UTR5) %>%
    dplyr::left_join(., median_tx)

  median_tx <- len_df %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(abs(utr3_len - median(utr3_len, na.rm = TRUE)) == min(abs(utr3_len - median(utr3_len, na.rm = TRUE)), na.rm = TRUE)) %>%
    dplyr::rename("UTR3" = tx_name) %>%
    dplyr::select(gene_id, UTR3) %>%
    dplyr::left_join(., median_tx)


  median_tx <- median_tx %>%
    dplyr::group_by(gene_id) %>%
    dplyr:: summarize(whole = dplyr::first(whole),
                      CDS = dplyr::first(CDS),
                      UTR5 = dplyr::first(UTR5),
                      UTR3 = dplyr::first(UTR3))

}

