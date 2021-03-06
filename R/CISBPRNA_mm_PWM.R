#' CISBPRNA Mus muscululs PWM list
#'
#' A list of position weight Matrices (PWM) of recognition motifs for RNA
#' binding proteins (RBPs) as documented by CISBPRNA for Mus musculus.
#'
#' @format A list of named matrices with 188 entries. The name is unique
#'   refering to both the RBP and its specific motif, the matrix columns
#'   represent positions, and the rows are the probability of occurance for each
#'   base "A", "C", "G", and "T". Each column sums to 1.
#'   \describe{
#'       \item{PWM}{probability matrix of a RBP motif}
#' }
#'
#' @source \url{http://cisbp-rna.ccbr.utoronto.ca/bulk.php}
#'     Download Mus musculus species archive
"CISBPRNA_mm_PWM"
