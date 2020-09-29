#' Export sequence of interest from a TxDb object
#'
#' \code{write_Sequence} writes gff3 files and/or fasta files from genome TxDb objects
#' subsetted with a transcript list. Sequences extracted can be the entire
#' transcripts (whole) or just the CDS, 3`UTR or 5`UTR sequences. If you instead
#' have sets of genes, we recommend using expression data to find expresed
#' transcripts within your experiment. Alternatively, use
#' \code{\link{gene2tx()}} to create transcript sets from either
#' \code{\link{make_longest_df}} or \code{\link{make_median_df}}
#' @param TxDb_gff A TxDb object. Using \code{\link{filter_Tx}} to create this
#'   TxDb object is recommended.
#' @param tx_list A character list. Must contain Ensembl transcript IDs and be
#'   compatible with the gff (the same species).
#' @param seq_type One of "whole", "CDS", "5pUTR", "3pUTR", or "promoter".
#' @param promoter_size numeric. default is 2000 bases upstream of transcription start.
#' @param file_name A single String. The file type (".fa", ".gff3") will be
#'   appended behind the file name.
#' @param output_type One of "fa", "gff3", or "both".
#' @return \code{write_Sequence} creates a named fasta and or gff3 file in the working
#'   directory.
#' @seealso \code{\link{filter_Tx}}, \code{\link{gene2tx}},
#' \code{\link{make_longest_df}}, \code{\link{make_median_df}}
#' @examples
#' #generate mydata/hs_test.fa and mydata/hs_test.gff3 containing CDS information for each transcript in hs_tx
#' hs_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.v33.annotation.gff3.gz", package = "FeatureReachR"))
#' hs_tx <- c("ENST00000456328.2", "ENST00000338338.9", "ENST00000356026.10", "ENST00000607222.1", "ENST00000342066.8")
#' write_Sequence(hs_filtered_TxDb, hs_tx, "CDS", "mydata/hs_test", "both")
#'
#' #generate mydata/mm_test.fa and mydata/mm_test.gff3 containing 5'UTR information for each transcript in mm_tx
#' mm_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"))
#' mm_tx <- c("ENSMUST00000159265.1", "ENSMUST00000027032.5", "ENSMUST00000130201.7", "ENSMUST00000157375.1")
#' write_Sequence(mm_filtered_TxDb, mm_tx, "UTR5", "mydata/mm_test", "both")
#' @export
write_Sequence <- function(TxDb_gff, tx_list, seq_type, file_name, output_type, promoter_size = 2000) {
  #Check that transcript list contains transcripts

  print("Ensuring compatibility of GFF and transcript list...", quote = FALSE)

  if (any(grepl("ENSG", tx_list)) | any(grepl("ENSMUSG", tx_list))) {
    stop("Please ensure the transcript list contains transcript IDs and not gene IDs.")
  }

  #Make sure gff and transcript list are the same species

  ifelse(all(grepl("ENST", tx_list)),
         (list_species = "human"),
         ifelse(all(grepl("ENSMUST", tx_list)),
                (list_species = "mouse"),
                (list_species = "unknown")))

  ifelse(grepl("ENSG", GenomicFeatures::genes(TxDb_gff)[1]$gene_id),
         (gff_species = "human"),
         ifelse(grepl("ENSMUSG",GenomicFeatures::genes(TxDb_gff)[1]$gene_id),
                (gff_species = "mouse"),
                (gff_species = "unknown")))

  if (list_species != gff_species) {
    stop("Please ensure the transcript list and gff are of the same species.")
  }

  print("GFF and transcript list are compatible.", quote = FALSE)

  #Subset the gffs by transcript and seq type

  tx_gff <- NULL

  if (seq_type == "whole") {
    print("filtering whole transcripts...", quote = FALSE)
    tx_gff <- GenomicFeatures::exonsBy(TxDb_gff, "tx", use.names = TRUE)
    tx_gff <- tx_gff[names(tx_gff) %in% tx_list]
    print("filtering complete.", quote = FALSE)

  } else if (seq_type == "CDS") {
    print("filtering CDS...", quote = FALSE)
    tx_gff <- GenomicFeatures::cdsBy(TxDb_gff, "tx", use.names = TRUE)
    tx_gff <- tx_gff[names(tx_gff) %in% tx_list]
    print("filtering complete.", quote = FALSE)

  } else if (seq_type == "UTR5") {
    print("filtering 5' UTRs...", quote = FALSE)
    tx_gff <- GenomicFeatures::fiveUTRsByTranscript(TxDb_gff, use.names = TRUE)
    tx_gff <- tx_gff[names(tx_gff) %in% tx_list]
    print("filtering complete.", quote = FALSE)

  } else if (seq_type == "UTR3") {
    print("filtering 3' UTRs...", quote = FALSE)
    tx_gff <- GenomicFeatures::threeUTRsByTranscript(TxDb_gff, use.names = TRUE)
    tx_gff <- tx_gff[names(tx_gff) %in% tx_list]
    print("filtering complete.", quote = FALSE)

  } else if (seq_type == "promoter") {
    print("filtering 3' UTRs...", quote = FALSE)
    tx_gff <- GenomicFeatures::promoters(TxDb_gff, upstream = promoter_size, downstream = 0, use.names = TRUE)
    tx_gff <- tx_gff[names(tx_gff) %in% tx_list]
    print("filtering complete.", quote = FALSE)

  } else
    stop("not appropriate seq type. please use one of \"whole\", \"CDS\", \"UTR5\", \"UTR3\" or, \"promoter\"")

  if(is.null(tx_gff) == TRUE | length(tx_gff) == 0) {
    stop("no sequences of type \"", seq_type, "\" found in transcript list", sep = "")
  }

  #Get the sequences from subset GFF

  if (list_species == "human" & gff_species == "human") {
    print("extracting sequences...", quote = FALSE)
    if (seq_type == "promoter") {
      seq <- GenomicFeatures::extractUpstreamSeqs(BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens, tx_gff)
    } else
      seq <- GenomicFeatures::extractTranscriptSeqs(BSgenome.Hsapiens.NCBI.GRCh38::Hsapiens, tx_gff)

  } else if (list_species == "mouse" & gff_species == "mouse") {
    print("extracting sequences...", quote = FALSE)
    if (seq_type == "promoter") {
      seq <- GenomicFeatures::extractUpstreamSeqs(BSgenome.Mmusculus.UCSC.mm10::Mmusculus, tx_gff)
    } else
      seq <- GenomicFeatures::extractTranscriptSeqs(BSgenome.Mmusculus.UCSC.mm10::Mmusculus, tx_gff)

  } else
    stop("Species type is unknown. Please ensure gff and transcript list are either mouse or human")

  names(seq) <- paste(names(seq), seq_type, sep = "_")

  print("extracting sequences complete.", quote = FALSE)

  #output sequences in desired format

  if (output_type == "fa") {
    print("writing fasta...", quote = FALSE)
    Biostrings::writeXStringSet(seq, paste(file_name, ".fa", sep = ""), format = "fasta")
    print("fasta writen.", quote = FALSE)

  } else if (output_type == "gff3") {
    print("writing GFF3...", quote = FALSE)
    rtracklayer::export(tx_gff, paste(file_name, ".gff3", sep = ""), format = "gff3")
    print("GFF3 written.", quote = FALSE)

  } else if (output_type == "both") {
    print("writing fasta...", quote = FALSE)
    Biostrings::writeXStringSet(seq, paste(file_name, ".fa", sep = ""), format = "fasta")
    print("fasta writen.", quote = FALSE)

    print("writing GFF3...", quote = FALSE)
    rtracklayer::export(tx_gff, paste(file_name, ".gff3", sep = ""), format = "gff3")
    print("GFF3 written.", quote = FALSE)

  } else
    stop("not appropriate format, please use \"fa\", \"gff3\" or \"both\"", quote = FALSE)


  print(paste("Querried ", length(tx_gff), " of ", length(tx_list), " requested genes.", sep = ""), quote = FALSE)
  print(paste("Querried ", sum(Biostrings::width(seq))/1000, " KB of sequence.", sep = ""), quote = FALSE)
}
