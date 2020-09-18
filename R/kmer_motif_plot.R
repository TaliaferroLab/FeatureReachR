#' ggplot Wrappers for Kmer/Motif Volcano Plots
#'
#' This creates Volcano plots comparing the enrichment of kmers or motifs in
#' case versus control groups. It plots the magnitude and significance of the
#' enrichment/depletion.
#'
#' @param case A DNAStringSet object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param ctrl A DNAStringSet object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return a ggplot2 visual
#' @seealso \code{kmer_compare}, \code{motif_compare},
#'   \code{estimate_kmer_from_motif}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
#' kmer_stats <- kmer_compare(case_fasta, ctrl_fasta)
#' kmer_plot(kmer_stats)
#'
#' RBNS_stats <- motif_compare(RBNS_PWM, case_fasta, ctrl_fasta)
#' motif_plot(RBNS_stats)
#' @export
#' @describeIn Plots all kmers with their magnitude (log2FoldChange) and
#'   significance
kmer_plot <- function(kmer_stats, sig_cutoff = 0.05){
  p <- kmer_stats %>% dplyr::mutate(sig = ifelse(p_adj < sig_cutoff, as.character(sig_cutoff), "ns"))
  if (length(unique(p$sig)) == 2) {
    p %>% ggplot2::ggplot(ggplot2::aes(x = log2FC, y = -log(p_adj), alpha = sig, col = sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = c("Red", "Black")) +
      ggplot2::scale_alpha_manual(values = c(1, 0.1)) +
      ggplot2::geom_text(data = subset(p, sig == as.character(sig_cutoff)), ggplot2::aes(label = kmer), nudge_y = 1) +
      cowplot::theme_cowplot()
  } else {
    p %>% ggplot2::ggplot(ggplot2::aes(x = log2FC, y = -log(p_adj), alpha = sig, col = sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = "Black") +
      ggplot2::scale_alpha_manual(values = 0.1) +
      cowplot::theme_cowplot()
    }
}
#' @export
#' @describeIn Plots all motifs with their magnitude (log2FoldChange) and
#'   significance.
motif_plot <- function(RBP_stats, sig_cutoff = 0.05){
  p <- RBP_stats %>% dplyr::mutate(sig = ifelse(p_adj < sig_cutoff, as.character(sig_cutoff), "ns"))
  if (length(unique(p$sig)) == 2) {
    p %>% ggplot2::ggplot(ggplot2::aes(x = log2FC, y = -log(p_adj), alpha = sig, col = sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = c("Red", "Black")) +
      ggplot2::scale_alpha_manual(values = c(1, 0.1)) +
      ggplot2::geom_text(data = subset(p, sig == as.character(sig_cutoff)), ggplot2::aes(label = motif), nudge_y = 1) +
      cowplot::theme_cowplot()
  } else {
    p %>% ggplot2::ggplot(ggplot2::aes(x = log2FC, y = -log(p_adj), alpha = sig, col = sig)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_manual(values = "Black") +
      ggplot2::scale_alpha_manual(values = 0.1) +
      cowplot::theme_cowplot()
  }
}
