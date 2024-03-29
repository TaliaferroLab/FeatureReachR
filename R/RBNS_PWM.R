#' RNA Bind and Seq PWM list
#'
#' A list of position weight Matrices (PWM) of recognition motifs for RNA
#' binding proteins (RBPs) as determined by RNA bind and seq RBNS experiments.
#' RNA bind and seq is an In vitro technique to define RBP binding sites. These
#' motifs were deined by RBNS experiments performed with 78 purified human RBPs
#' and random synthetic RNAs.
#'
#' @format A list of named matrices with 131 entries. The name is unique
#'   refering to both the RBP and its specific motif, the matrix columns
#'   represent positions, and the rows are the probability of occurance for each
#'   base "A", "C", "G", and "T". Each column sums to 1. \describe{
#'   \item{PWM}{probability matrix of a RBP motif} }
#'
#' @source Dominguez, D., Freese, P., Alexis, M. S., Su, A., Hochman, M.,
#'   Palden, T., Bazile, C., Lambert, N. J., Van Nostrand, E. L., Pratt, G. A.,
#'   Yeo, G. W., Graveley, B. R., & Burge, C. B. (2018). Sequence, Structure,
#'   and Context Preferences of Human RNA Binding Proteins. Molecular cell,
#'   70(5), 854–867.e9. https://doi.org/10.1016/j.molcel.2018.05.001
#'   \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6062212/}
"RBNS_PWM"
