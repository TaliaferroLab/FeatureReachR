#' Group Kmers and Represent them as Motif Logos
#'
#' Kmers of interest can be collapsed or grouped by their levenshtein distance
#' (>2) to simplify enriched/depleted kmer outputs. These functions help visualize
#' how kmers are grouped as well as provide PWM which can be visualized using
#' \code{ggseqlogo}.
#'
#' @param kmers character list of kmers of the length.
#' @param title string for plot titles defaults to no title.
#' @return a tree, seqotif plots or a named list of PWMs
#' @seealso
#' \code{kmer_compare}, \code{\link[ggseqlogo]{ggseqlogo}}
#' @examples
#' enriched_sixmers <- c("AAGGAA", "ACACAC", "AGAAGG", "AGAGAG", "AGAGGG",
#' "AGGAAG", "AGGAGG", "AGGGAG", "CACACA", "GAAGGA", "GAGAAG", "GAGAGA",
#' "GAGGAG", "GAGGGA", "GAGGGG", "GGAAGG", "GGAGGA", "GGAGGG", "GGGAGG")
#' kmer2logo(enriched_sixmers, "My Plot Title")
#' kmer2tree(enriched_sixmers, "My Tree Title")
#' kmer_PWMs <- kmer2PWM(enriched_sixmers)
#' kmer_PWMs
#'
#' @describeIn plots the tree and all seq logos together
#' @export
#'
kmer2logo <- function(kmers, title = ""){
  d  <- adist(kmers) #this is Levenshtein distance
  rownames(d) <- kmers
  hc <- hclust(as.dist(d))
  plot(hc, hang = -1, main = title, xlab = "", sub = "")

  if (max(hc$height) <= 2){
    print("All kmers are within a Levenshtein distance of 2")
    Biostrings::DNAStringSet(kmers) %>%
      Biostrings::pairwiseAlignment(., .[[1]], type = "overlap") %>%
      Biostrings::consensusMatrix(endgapCode = "N") %>%
      t() %>%
      dplyr::as_tibble(rownames = "pos") %>%
      dplyr::mutate(A = A + 0.25*N,
                    C = C + 0.25*N,
                    G = G + 0.25*N,
                    T = T + 0.25*N) %>%
      dplyr::select(A,C,G,T) %>%
      t() %>%
      prop.table(.,2) %>%
      ggseqlogo::ggseqlogo(., method = "prob", ncol = 3) +
      ggplot2::ggtitle(title)
  } else{
    rh <- rect.hclust(hc, h=2.5) #groups with levenshtein distance < 2
    names(rh) <- as.character(c(1:length(rh)))

    beg_clus <- head(cumsum(c(1, lengths(rh))), -1)
    text(x = beg_clus, y = 3, col = "red", labels = names(rh), font = 2, cex = 2)

    df <- lapply(rh, names) %>%
      tibble::enframe() %>%
      tidyr::unnest(cols = "value") %>%
      dplyr::rename("group" = name, "kmer" = value)

    num_grps <- df %>%
      dplyr::pull(group) %>%
      unique() %>%
      length()

    x <- c(1: num_grps)

    list <- lapply(x, function(x) df %>%
                     dplyr::filter(group == x) %>%
                     dplyr::pull(kmer) %>%
                     as.character() %>%
                     Biostrings::DNAStringSet() %>%
                     Biostrings::pairwiseAlignment(., .[[1]], type = "overlap") %>%
                     Biostrings::consensusMatrix(endgapCode = "N") %>%
                     t() %>%
                     dplyr::as_tibble(rownames = "pos") %>%
                     dplyr::mutate(A = A + 0.25*N,
                                   C = C + 0.25*N,
                                   G = G + 0.25*N,
                                   T = T + 0.25*N) %>%
                     dplyr::select(A,C,G,T) %>%
                     t() %>%
                     prop.table(.,2))

    ggseqlogo::ggseqlogo(list, method = "prob", ncol = 3) + ggplot2::ggtitle(title)
  }
}
#'
#' @describeIn plots the tree of kmers with groups highlighted.
#' @export
#'
kmer2tree <- function(kmers, title = ""){
  d  <- adist(kmers) #this is Levenshtein distance
  rownames(d) <- kmers
  hc <- hclust(as.dist(d))
  plot(hc, hang = -1, main = title, xlab = "", sub = "")

  rh <- rect.hclust(hc, h=2.5) #groups with levenshtein distance < 2
  names(rh) <- as.character(c(1:length(rh)))

  beg_clus <- head(cumsum(c(1, lengths(rh))), -1)
  text(x = beg_clus, y = 3, col = "red", labels = names(rh), font = 2, cex = 2)
}
#'
#' @describeIn creates a list of PWMs used to create seq logos. The groups match
#'   with groups plotted on trees.
#' @export
#'
kmer2PWM <- function(kmers){
  d  <- adist(kmers) # this is Levenshtein distance
  rownames(d) <- kmers
  hc <- hclust(as.dist(d))

  if (max(hc$height) <= 2){
    print("All kmers are within a Levenshtein distance of 2")
    list <- Biostrings::DNAStringSet(kmers) %>%
      Biostrings::pairwiseAlignment(., .[[1]], type = "overlap") %>%
      Biostrings::consensusMatrix(endgapCode = "N") %>%
      t() %>%
      dplyr::as_tibble(rownames = "pos") %>%
      dplyr::mutate(A = A + 0.25*N,
                    C = C + 0.25*N,
                    G = G + 0.25*N,
                    T = T + 0.25*N) %>%
      dplyr::select(A,C,G,T) %>%
      t() %>%
      prop.table(.,2)
  } else{
    rh <- rect.hclust(hc, h=2.5, border = "#00000000") # Levenshtein distance < 2
    names(rh) <- as.character(c(1:length(rh)))

    df <- lapply(rh, names) %>%
      tibble::enframe() %>%
      tidyr::unnest(cols = "value") %>%
      dplyr::rename("group" = name, "kmer" = value)

    num_grps <- df %>%
      dplyr::pull(group) %>%
      unique() %>%
      length()
    x <- c(1: num_grps)

    list <- lapply(x, function(x) df %>%
                     dplyr::filter(group == x) %>%
                     dplyr::pull(kmer) %>%
                     as.character() %>%
                     Biostrings::DNAStringSet() %>%
                     Biostrings::pairwiseAlignment(., .[[1]], type = "overlap") %>%
                     Biostrings::consensusMatrix(endgapCode = "N") %>%
                     t() %>%
                     dplyr::as_tibble(rownames = "pos") %>%
                     dplyr::mutate(A = A + 0.25*N,
                                   C = C + 0.25*N,
                                   G = G + 0.25*N,
                                   T = T + 0.25*N) %>% # shifts treated as N (equal A, C, G, T)
                     dplyr::select(A,C,G,T) %>%
                     t() %>%
                     prop.table(.,2))
  }
  return(list)
}

