#' custom PWM list examlpe
#'
#' A list of position weight Matrices (PWM) of recognition motifs for RNA
#' binding proteins (RBPs) to be used as an example
#'
#' @format A list of named matrices with 5 entries. The name is unique refering
#'   to the RBP, the matrix columns represent positions, and the rows are the
#'   probability of occurance for each base "A", "C", "G", and "T". Each column
#'   sums to 1.
#'    \describe{ \item{PWM}{probability matrix of a RBP motif} }
#'
#' @source a subset of RBNS_PWM
"custom_PWM_list"
