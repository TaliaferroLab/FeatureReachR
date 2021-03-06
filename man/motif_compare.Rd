% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Motif_compare.R
\name{motif_compare}
\alias{motif_compare}
\title{Compare Motif occurances between two DNAStringSets}
\usage{
motif_compare(PWM_list, caseDNAStringset, ctrlDNAStringSet)
}
\arguments{
\item{PWM_list}{a named list of position weight matrices. A list of named
matrices where the name is unique, and the matrix columns represent
positions, and the rows are the probability of each base "A", "C", "G", and
"T". Each column should sum to 1. FeatureReachR has three PWM_lists build in:
"CISBPRNA_mm_PWM", "CISBPRNA_hs_PWM", and "RBNS_PWM".}

\item{ctrlDNAStringSet}{A DNAStringSet Object. Use
\code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.}

\item{caseDNAStringSet}{A DNAStringSet Object. Use
\code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.}
}
\value{
a summary table including the pvalue as calculated by
  \code{wilcox.test}, then corrected by \code{p.adjust} with the \code{method = "BH"}.
  The mean values for each DNAStringSet, and the log2 fold
  change between those means.
}
\description{
This function compares motif occurances between between case and control
DNAStringSets.
}
\examples{
case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")

motif_compare(CISBPRNA_mm_PWM, case_fasta, ctrl_fasta)
motif_compare(RBNS_PWM, case_fasta, ctrl_fasta)

}
\seealso{
\code{\link{write_Sequence}}, \code{\link[Biostrings]{readDNAStringSet}},
  \code{\link[Biostrings]{countPWM}}, \code{\link{motif_by_gene}}
}
