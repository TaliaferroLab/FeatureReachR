#' Summarize RNA types in transcript list
#'
#' \code{RNA_breakdown} counts RNA types present in a transcript list and can be
#' used to better understand your case and control sets of transcripts. If you
#' instead have sets of genes, we recommend using expression data to find
#' expresed transcripts within your experiment. Alternatively, use
#' \code{\link{gene2tx()}} to create transcript sets from either
#' \code{\link{make_longest_df}} or \code{\link{make_median_df}}
#'
#' @param gff A path. To the human or mouse gencode annotation gff file. The
#'   most current annotations can be downloaded from:
#'   \url{https://www.gencodegenes.org/human/} or
#'   \url{https://www.gencodegenes.org/mouse/} There are two filter options, a
#'   general filter and a protein coding filter.
#' @param tx_list A character list. Must contain Ensembl transcript IDs and be
#'   compatible with the gff (the same species)
#' @return \code{RNA_breakdown} creates a summary table containing both the
#' number of transcripts categorized within an RNA type (count) and the
#' percent of the total set an RNA type represents (percent).
#' @seealso \code{\link{gene2tx()}},
#' \code{\link{make_longest_df}}, \code{\link{make_median_df}}
#' @importFrom magrittr "%>%"
#' @examples
#' mm_tx <- c("ENSMUST00000159265.1", "ENSMUST00000027032.5", "ENSMUST00000130201.7", "ENSMUST00000157375.1")
#' RNA_breakdown("mydata/Gencodedat/gencode.vM20.annotation.gff3.gz", mm_tx)
#' @export
RNA_breakdown <- function(gff, tx_list) {
  # Check that tx_list doesn't contain Ensembl gene IDs
  print("Ensuring compatibility of GFF and transcript list...", quote = FALSE)
  if (any(grepl("ENSG", tx_list)) | any(grepl("ENSMUSG", tx_list))) {
    stop("Please ensure the transcript list contains transcript IDs and not gene IDs.")
  }
  # Determine the species of tx_list
  ifelse(all(grepl("ENST", tx_list)),
         (list_species = "human"),
         ifelse(all(grepl("ENSMUST", tx_list)),
                (list_species = "mouse"),
                (list_species = "unknown")))
  # Import gencode gff as GRanges object
  gff <- rtracklayer::import(gff)
  # Determine species of gff annotations
  gff_species = "unknown"
  if (all(grepl("ENSG", gff$gene_id))) {
    gff_species = "human"
  } else if (all(grepl("ENSMUSG", gff$gene_id))) {
    gff_species = "mouse"
  }
  # Ensure tx_list and gff are of the same species
  if (list_species != gff_species) {
    stop("Please ensure the transcript list and gff are of the same species.")
  }
  print("Gff and transcript list are compatible.", quote = FALSE)

  print("Getting RNA type information...")
  gff <- gff[gff$type == "transcript"]
  gff <- gff[gff$transcript_id %in% tx_list]

  if (length(gff) == 0) {
    stop("No transcript IDs in tx_list were not found in gff.")
  }

  RNA_breakdown <- gff@elementMetadata %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(transcript_type) %>%
    dplyr::summarize(count = n(),
                     percent = n() / nrow(.) * 100)

  return(RNA_breakdown)
}
