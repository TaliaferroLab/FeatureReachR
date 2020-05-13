#' RNA Bind and Seq PWM list
#'
#' A list of position weight Matrices (PWM) of recognition motifs for RNA
#' binding proteins (RBPs) as determined by RNA bind and seq RBNS experiments.
#' RNA bind and seq is an In vitro technique to define RBP binding sites.
#' These assays are performed wiht purified RBP and random synthetic RNAs.
#'
#' @format A list of named matrices with 131 entries. The name is unique
#'   refering to both the RBP and its specific motif, the matrix columns
#'   represent positions, and the rows are the probability of occurance for each
#'   base "A", "C", "G", and "T". Each column sums to 1.
#'   \describe{
#'       \item{PWM}{probability matrix of a RBP motif} }
#'
#' @source Matthew Taliaferro's friend D.Dominguez
"RBNS_PWM"
