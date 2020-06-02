## code to prepare kmer motif relation datasets goes here

relate_kmer_PWM <- function(k, PWM_list){

  ##make DNAStringSets for every kmer
  x <- c("T", "A", "C", "G")

  kmers <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer) %>%
    Biostrings::DNAStringSet()

  names(kmers) <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer)

  #Check PWM_list
  print("Checking PWM data", quote = FALSE)

  if (any(as.numeric(lapply(PWM_list, nrow)) != 4)){
    stop("Error: Ensure all PWM matricies have exactly 4 rows (\"A\", \"C\", \"G\" and \"T\")")
  }

  #this function makes every PWM the same length as k or shorter by tiling across longer PWM matrices
  mattiles <- function(Mat, k) {

    if(ncol(Mat[[1]]) > k){
      x <- c(1:(ncol(Mat[[1]])-(k-1)))-1
      list <- lapply(x, function(x) Mat[[1]][,(1+x):(k+x)])
      names(list) <- lapply(x, function(x) paste(names(Mat), "_", x + 1, sep = ""))
      list
    }
    else
      Mat
  }

  x <- 1:length(PWM_list)
  tiled_motifs <- unlist(lapply(x, function(x) mattiles(PWM_list[x], k)), recursive = FALSE)


  print("Counting motif occurances...", quote = FALSE)

  counts_list <- lapply(tiled_motifs, function(x)
    lapply(kmers, function(y)
      Biostrings::countPWM(x, y)) %>% unlist() %>% dplyr::as_tibble(rownames = "kmer"))

  print("Counting motif occurances complete.", quote = FALSE)

  motif_by_kmer <- dplyr::bind_rows(counts_list, .id = "motif") %>%
    tidyr::spread(kmer, value = value) %>%
    dplyr::select(motif, dplyr::everything())

  motif_by_kmer <- motif_by_kmer %>% dplyr::mutate_if(is.numeric, ~1 * (. > 0))

  return(motif_by_kmer)
}


fourmer_mm_cisbpRNA <- relate_Kmer_PWM(4, CISBPRNA_mm_PWM)
fivemer_mm_cisbpRNA <- relate_Kmer_PWM(5, CISBPRNA_mm_PWM)
ixmer_mm_cisbpRNA <- relate_Kmer_PWM(6, CISBPRNA_mm_PWM)
sevenmer_mm_cisbpRNA <- relate_Kmer_PWM(7, CISBPRNA_mm_PWM)

fourmer_hs_cisbpRNA <- relate_Kmer_PWM(4, CISBPRNA_hs_PWM)
fivemer_hs_cisbpRNA <- relate_Kmer_PWM(5, CISBPRNA_hs_PWM)
sixmer_hs_cisbpRNA <- relate_Kmer_PWM(6, CISBPRNA_hs_PWM)
sevenmer_hs_cisbpRNA <- relate_Kmer_PWM(7, CISBPRNA_hs_PWM)

fourmer_RBNS <- relate_Kmer_PWM(4, RBNS_PWM)
fivemer_RBNS <- relate_Kmer_PWM(5, RBNS_PWM)
sixmer_RBNS <- relate_Kmer_PWM(6, RBNS_PWM)
sevenmer_RBNS <- relate_Kmer_PWM(7, RBNS_PWM)


usethis::use_data(fourmer_mm_cisbpRNA, fivemer_mm_cisbpRNA, sixmer_mm_cisbpRNA, sevenmer_mm_cisbpRNA, fourmer_hs_cisbpRNA, fivemer_hs_cisbpRNA, sixmer_hs_cisbpRNA, sevenmer_hs_cisbpRNA, fourmer_RBNS, fivemer_RBNS, sixmer_RBNS, sevenmer_RBNS, internal = TRUE, overwrite = TRUE)


