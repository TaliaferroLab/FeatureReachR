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
#' \code{kmer_compare}, \code{\link[ggseqlogo]{ggseqlogo}}, \code{\link[msa]{msa}}
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

  num_clust <- length(unique(cutree(hc,h = 2.5)))
  num_max_members <- cutree(hc,h = 2.5) %>%
    as_tibble(rownames = "kmer") %>%
    group_by(value) %>%
    summarize(n=n(), .groups = "drop") %>%
    pull(n) %>%
    max()

  if (num_clust == 1){
    print("All kmers are within a Levenshtein distance of 2")
    Biostrings::DNAStringSet(kmers) %>%
      msa::msa() %>%
      Biostrings::consensusMatrix() %>%
      t() %>%
      dplyr::as_tibble(rownames = "pos") %>%
      dplyr::rename("gap" = `-`) %>%
      dplyr::mutate(A = A + 0.25*gap,
                    C = C + 0.25*gap,
                    G = G + 0.25*gap,
                    T = T + 0.25*gap) %>%
      dplyr::select(A,C,G,T) %>%
      t() %>%
      prop.table(.,2) %>%
      ggseqlogo::ggseqlogo(., method = "prob", ncol = 3) +
      ggplot2::ggtitle(title)
  } else if (num_max_members == 1) {
    print("No kmers are within a Levenshtein distance of 2, each kmer is its own cluster")
    Biostrings::DNAStringSet(kmers)
    Matrix_list <- lapply(c(1:length(kmers)), function(x) consensusMatrix(Biostrings::DNAStringSet(kmers)[x])[1:4,])
    names(Matrix_list) <- kmers
    ggseqlogo::ggseqlogo(Matrix_list, method = "prob", ncol = 3) + ggplot2::ggtitle(title)

  } else if (num_clust > 1 & num_max_members > 1) {
    plot(hc, hang = -1, main = title, xlab = "", sub = "")
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

    kmer_grp_list <- lapply(x, function(x) df %>%
                     dplyr::filter(group == x) %>%
                     dplyr::pull(kmer) %>%
                     as.character() %>%
                     Biostrings::DNAStringSet())



   Matrix_list <- lapply(kmer_grp_list, function(kmer_grp)
                  if (length(kmer_grp) > 1){
                     kmer_grp %>%
                       msa::msa() %>%
                     Biostrings::consensusMatrix() %>%
                     t() %>%
                     dplyr::as_tibble(rownames = "pos") %>%
                     dplyr::rename("gap" = `-`) %>%
                     dplyr::mutate(A = A + 0.25*gap,
                                   C = C + 0.25*gap,
                                   G = G + 0.25*gap,
                                   T = T + 0.25*gap) %>%
                     dplyr::select(A,C,G,T) %>%
                     t() %>%
                     prop.table(.,2)
                     } else {kmer_grp %>%
                         Biostrings::consensusMatrix() %>%
                         t()%>%
                         dplyr::as_tibble(rownames = "pos") %>%
                         dplyr::select(A,C,G,T) %>%
                         t() %>%
                         prop.table(.,2)
                       })

    ggseqlogo::ggseqlogo(Matrix_list, method = "prob", ncol = 3) + ggplot2::ggtitle(title)
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

  num_clust <- length(unique(cutree(hc,h = 2.5)))
  num_max_members <- cutree(hc,h = 2.5) %>%
    as_tibble(rownames = "kmer") %>%
    group_by(value) %>%
    summarize(n=n(), .groups = "drop") %>%
    pull(n) %>%
    max()

  if (num_clust == 1) {
    warning("all kmers are within the same custer and cannot be plotted as a tree")
  } else if(num_max_members == 1) {
    warning("each kmer is in its own cluster and cannot be plotted as a tree")
  } else if (num_clust > 1 & num_max_members > 1) {
    plot(hc, hang = -1, main = title, xlab = "", sub = "")

    rh <- rect.hclust(hc, h=2.5) #groups with levenshtein distance < 2
    names(rh) <- as.character(c(1:length(rh)))

    beg_clus <- head(cumsum(c(1, lengths(rh))), -1)
    text(x = beg_clus, y = 3, col = "red", labels = names(rh), font = 2, cex = 2)
  } else {
    warning("kmers not suitable for plotting tree or more kmers must be included.")
  }
}
#'
#' @describeIn creates a list of PWMs used to create seq logos. The groups match
#'   with groups plotted on trees.
#' @export
#'
kmer2PWM <- function(kmers){
  d  <- adist(kmers) #this is Levenshtein distance
  rownames(d) <- kmers
  hc <- hclust(as.dist(d))

  num_clust <- length(unique(cutree(hc,h = 2.5)))
  num_max_members <- cutree(hc,h = 2.5) %>%
    as_tibble(rownames = "kmer") %>%
    group_by(value) %>%
    summarize(n=n(), .groups = "drop") %>%
    pull(n) %>%
    max()

  if (num_clust == 1){
    print("All kmers are within a Levenshtein distance of 2")
    Matrix <- Biostrings::DNAStringSet(kmers) %>%
      msa::msa() %>%
      Biostrings::consensusMatrix() %>%
      t() %>%
      dplyr::as_tibble(rownames = "pos") %>%
      dplyr::rename("gap" = `-`) %>%
      dplyr::mutate(A = A + 0.25*gap,
                    C = C + 0.25*gap,
                    G = G + 0.25*gap,
                    T = T + 0.25*gap) %>%
      dplyr::select(A,C,G,T) %>%
      t() %>%
      prop.table(.,2)
    return(Matrix)

  } else if (num_max_members == 1) {
    print("No kmers are within a Levenshtein distance of 2, each kmer is its own cluster")
    Biostrings::DNAStringSet(kmers)
    Matrix_list <- lapply(c(1:length(kmers)), function(x) consensusMatrix(Biostrings::DNAStringSet(kmers)[x])[1:4,])
    names(Matrix_list) <- kmers
    return(Matrix_list)

    } else if (num_clust > 1 & num_max_members > 1) {
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

    kmer_grp_list <- lapply(x, function(x) df %>%
                              dplyr::filter(group == x) %>%
                              dplyr::pull(kmer) %>%
                              as.character() %>%
                              Biostrings::DNAStringSet())



    Matrix_list <- lapply(kmer_grp_list, function(kmer_grp)
      if (length(kmer_grp) > 1){
        kmer_grp %>%
          msa::msa() %>%
          Biostrings::consensusMatrix() %>%
          t() %>%
          dplyr::as_tibble(rownames = "pos") %>%
          dplyr::rename("gap" = `-`) %>%
          dplyr::mutate(A = A + 0.25*gap,
                        C = C + 0.25*gap,
                        G = G + 0.25*gap,
                        T = T + 0.25*gap) %>%
          dplyr::select(A,C,G,T) %>%
          t() %>%
          prop.table(.,2)
      } else {kmer_grp %>%
          Biostrings::consensusMatrix() %>%
          t()%>%
          dplyr::as_tibble(rownames = "pos") %>%
          dplyr::select(A,C,G,T) %>%
          t() %>%
          prop.table(.,2)
      })
  }
  return(Matrix_list)
}

