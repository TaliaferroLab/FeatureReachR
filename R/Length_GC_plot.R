#' ggplot Wrappers for length/GC Violin Plots
#'
#' This creates Violin plots comparing the mean length or GC of  case and
#' control groups. It displays sample size and the p value of a wilcox test.
#'
#' @param case A DNAStringSet object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @param ctrl A DNAStringSet object. Use
#'   \code{\link[Biostrings]{readDNAStringSet}} on a fasta file to create.
#' @return a ggplot2 visual
#' @seealso \code{get_length}, \code{get_GC}
#' @examples
#' case_fasta <- Biostrings::readDNAStringSet("example_case.fa")
#' ctrl_fasta <- Biostrings::readDNAStringSet("example_ctrl.fa")
#' length_plot(case_fasta, ctrl_fasta)
#' GC_plot(case_fasta, ctrl_fasta)
#' @export
#' @describeIn Plots length of all case and ctrl sequences as violin plots.
length_plot <- function(case, ctrl){
  case_length <- RNAreachR::get_length(case) %>%
    dplyr::mutate(group = "case",
                  label = paste("case \nn = ", nrow(.), sep = ""))

  ctrl_length <- RNAreachR::get_length(ctrl) %>%
    dplyr::mutate(group = "ctrl",
                  label = paste("ctrl \nn = ",nrow(.), sep = ""))

  length <- rbind(case_length, ctrl_length)

  length %>% ggplot2::ggplot(ggplot2::aes(x = factor(label), y = length, fill = factor(label))) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.25) +
    cowplot::theme_cowplot() +
    ggpubr::stat_compare_means(method = "wilcox") +
    ggplot2::guides(fill = FALSE) +
    ggplot2::labs(x = "", y = "length")
}
#' @export
#' @describeIn Plots GC content of all case and ctrl sequences as violin plots.
GC_plot <- function(case, ctrl){
  case_GC <- RNAreachR::get_GC(case) %>%
    dplyr::mutate(group = "case",
                  label = paste("case \nn = ", nrow(.), sep = ""))

  ctrl_GC <- RNAreachR::get_GC(ctrl) %>%
    dplyr::mutate(group = "ctrl",
                  label = paste("ctrl \nn = ",nrow(.), sep = ""))

  GC <- rbind(case_GC, ctrl_GC)

  GC %>% ggplot2::ggplot(ggplot2::aes(x = factor(label), y = GC, fill = factor(label))) +
    ggplot2::geom_violin() +
    ggplot2::geom_boxplot(width = 0.25) +
    cowplot::theme_cowplot() +
    ggpubr::stat_compare_means(method = "wilcox") +
    ggplot2::guides(fill = FALSE) +
    ggplot2::labs(x = "", y = "Percent GC Content")
}
