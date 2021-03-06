% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTxout.R
\name{write_Sequence}
\alias{write_Sequence}
\title{Export sequence of interest from a TxDb object}
\usage{
write_Sequence(TxDb_gff, tx_list, seq_type, file_name, output_type)
}
\arguments{
\item{TxDb_gff}{A TxDb object. Using \code{\link{filter_Tx}} to create this
TxDb object is recommended.}

\item{tx_list}{A character list. Must contain Ensembl transcript IDs and be
compatible with the gff (the same species).}

\item{seq_type}{One of "whole", "CDS", "5pUTR", or "3pUTR".}

\item{file_name}{A single String. The file type (".fa", ".gff3") will be
appended behind the file name.}

\item{output_type}{One of "fa", "gff3", or "both".}
}
\value{
\code{write_Sequence} creates a named fasta and or gff3 file in the working
  directory.
}
\description{
\code{write_Sequence} writes gff3 files and/or fasta files from genome TxDb objects
subsetted with a transcript list. Sequences extracted can be the entire
transcripts (whole) or just the CDS, 3`UTR or 5`UTR sequences. If you instead
have sets of genes, we recommend using expression data to find expresed
transcripts within your experiment. Alternatively, use
\code{\link{gene2tx()}} to create transcript sets from either
\code{\link{make_longest_df}} or \code{\link{make_median_df}}
}
\examples{
#generate mydata/hs_test.fa and mydata/hs_test.gff3 containing CDS information for each transcript in hs_tx
hs_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.v33.annotation.gff3.gz", package = "FeatureReachR"))
hs_tx <- c("ENST00000456328.2", "ENST00000338338.9", "ENST00000356026.10", "ENST00000607222.1", "ENST00000342066.8")
write_Sequence(hs_filtered_TxDb, hs_tx, "CDS", "mydata/hs_test", "both")

#generate mydata/mm_test.fa and mydata/mm_test.gff3 containing 5'UTR information for each transcript in mm_tx
mm_filtered_TxDb <- filter_Tx(system.file("extdata", "gencode.vM20.annotation.gff3.gz", package = "FeatureReachR"))
mm_tx <- c("ENSMUST00000159265.1", "ENSMUST00000027032.5", "ENSMUST00000130201.7", "ENSMUST00000157375.1")
write_Sequence(mm_filtered_TxDb, mm_tx, "UTR5", "mydata/mm_test", "both")
}
\seealso{
\code{\link{filter_Tx}}, \code{\link{gene2tx}},
\code{\link{make_longest_df}}, \code{\link{make_median_df}}
}
