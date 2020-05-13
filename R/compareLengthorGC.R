#' Compare the average Length or GC content of two DNAStringSet objects
#'
#' \code{length_compare} and \code{GC_compare} compare the length and GC content
#' between case and control DNAStringSets respectively.
#'
#' @param caseDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param  ctrlDNAStringSet A DNAStringSet Object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return a summary table including the pvalue as calculated by
#'   \code{wilcox.test}, the mean values for each DNAStringSet, the log2 Fold
#'   change between those means, the effect size as calculated by
#'   \code{\link[effsize]{cliff.delta}}
#' @seealso \code{\link{getTxOut}}, \code{\link[Biostrings]{readDNAStringSet}},
#'   \code{\link[Biostrings]{letterFrequency}},
#'   \code{\link[effsize]{cliff.delta}}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
#' length_compare(case_fasta, ctrl_fasta)
#' GC_compare(case_fasta, ctrl_fasta)
#' @export
#' @describeIn compares the average length of two different DNAStringSet objects
#'   and returns a summary statistic table of the comparison.


length_compare <- function(caseDNAStringSet, ctrlDNAStringSet){

  if (any(names(caseDNAStringSet) %in% names(ctrlDNAStringSet))){
    warning("some sequences in case set are also in the control set. This is not recommended.")
  }

  wilcox.p <- wilcox.test(Biostrings::width(caseDNAStringSet), Biostrings::width(ctrlDNAStringSet))$p.value
  mean_case <- mean(Biostrings::width(caseDNAStringSet))
  mean_ctrl <- mean(Biostrings::width(ctrlDNAStringSet))
  mean_FC <- mean_case/mean_ctrl
  CliffDelta <- effsize::cliff.delta(Biostrings::width(caseDNAStringSet), Biostrings::width(ctrlDNAStringSet))$estimate
  lowerCD <- effsize::cliff.delta(Biostrings::width(caseDNAStringSet), Biostrings::width(ctrlDNAStringSet))$conf.int[1]
  upperCD <- effsize::cliff.delta(Biostrings::width(caseDNAStringSet), Biostrings::width(ctrlDNAStringSet))$conf.int[2]

  data.frame(wilcox.p, mean_case, mean_ctrl, mean_FC, CliffDelta, lowerCD, upperCD)

}

#' @describeIn compares the average GC content of two different DNAStringSet
#'   objects and returns a summary statistic table of the comparison.

GC_compare <- function(caseDNAStringSet, ctrlDNAStringSet){

  if (any(names(caseDNAStringSet) %in% names(ctrlDNAStringSet))){
    warning("some sequences in case set are also in the control set. This is not recommended.")
  }

  GC_case <- Biostrings::letterFrequency(caseDNAStringSet, "GC") / Biostrings::width(caseDNAStringSet)
  GC_ctrl <- Biostrings::letterFrequency(ctrlDNAStringSet, "GC") / Biostrings::width(ctrlDNAStringSet)
  wilcox.p <- wilcox.test(GC_case, GC_ctrl)$p.value
  mean_case <- mean(GC_case)
  mean_ctrl <- mean(GC_ctrl)
  mean_FC <- mean_case/mean_ctrl
  CliffDelta <- effsize::cliff.delta(GC_case, GC_ctrl)$estimate
  lowerCD <- effsize::cliff.delta(GC_case, GC_ctrl)$conf.int[1]
  upperCD <- effsize::cliff.delta(GC_case, GC_ctrl)$conf.int[2]

  data.frame(wilcox.p, mean_case, mean_ctrl, mean_FC, CliffDelta, lowerCD, upperCD)
}
