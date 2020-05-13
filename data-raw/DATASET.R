## code to prepare kmer motif relation datasets goes here

relate_Kmer_RBNS <- function(k, motif_path){

  ##make DNAStringSets for every kmer
  x <- c("T", "A", "C", "G")

  kmers <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer) %>%
    Biostrings::DNAStringSet()

  names(kmers) <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer)

  #get paths to each motif pwm
  print("getting Motif and RBP data...", quote = FALSE)

  motif_paths <- list.files(path = motif_path, full.names = TRUE)
  motif_info <- file.info(path = motif_paths)
  motifs <- motif_info %>%
    dplyr::as_tibble(rownames = "PATH") %>%
    dplyr::mutate(motif = stringr::str_match(PATH, "PWMs/(.*?).PWM")[,2]) %>%
    dplyr::select(PATH, motif)

  print("Getting PWM data", quote = FALSE)

  motifs <- motifs %>%
    mutate(PWM = lapply(PATH, function(x) t(read.table(x, skip = 1, row.names = 1, header = TRUE, col.names = c("pos","A", "C", "G", "T")))))


  PWM_list <- motifs$PWM
  names(PWM_list) <- motifs$motif

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

  motif_by_kmer <- motif_by_kmer %>% dplyr::as_tibble((. > 0) + 0) %>% dplyr::mutate(motif = .$motif)
  return(motif_by_kmer)
}
fourmer_RBNS <- relate_Kmer_RBNS(4, "mydata/RBNSdat/RBNS_PWMs")
fivemer_RBNS <- relate_Kmer_RBNS(5, "mydata/RBNSdat/RBNS_PWMs")
sixmer_RBNS <- relate_Kmer_RBNS(6, "mydata/RBNSdat/RBNS_PWMs")
sevenmer_RBNS <- relate_Kmer_RBNS(7, "mydata/RBNSdat/RBNS_PWMs")

usethis::use_data(fourmer_RBNS, fivemer_RBNS, sixmer_RBNS, sevenmer_RBNS, internal = TRUE, overwrite = TRUE)

## code to prepare cisbp_RBNS dataset goes here

relate_Kmer_CisBPRNA <- function(k, motif_path, RBPinfo){

  ##make DNAStringSets for every kmer
  x <- c("T", "A", "C", "G")

  kmers <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer) %>%
    Biostrings::DNAStringSet()

  names(kmers) <- do.call(expand.grid, rep(list(x), k)) %>%
    tidyr::unite(kmer, sep = "") %>%
    dplyr::pull(., kmer)

  #get paths to each motif pwm
  print("getting Motif and RBP data...", quote = FALSE)

  motif_paths <- list.files(path = motif_path, full.names = TRUE)
  motif_info <- file.info(motif_paths)
  motif_info <- motif_info[motif_info$size != 0, ]
  motifs <- motif_info %>%
    dplyr::as_tibble(rownames = "PATH") %>%
    dplyr:: mutate(motif = stringr::str_match(PATH, "motifs/(.*?).txt")[,2]) %>%
    dplyr::select(PATH, motif)

  #merge motif paths with RBP info

  RBP_info <- read.table(RBPinfo, header = TRUE, sep = "\t")
  RBP_info <- RBP_info %>%
    dplyr::as_tibble() %>%
    dplyr::select(Motif_ID, RBP_Name) %>%
    dplyr::filter(Motif_ID != ".") %>%
    dplyr::group_by(Motif_ID) %>%
    dplyr::summarise(RBP_name = dplyr::first(RBP_Name))

  motifs <- dplyr::left_join(motifs, RBP_info, by = c("motif" = "Motif_ID"))

  print("Getting PWM data", quote = FALSE)

  motifs <- motifs %>%
    dplyr::mutate(PWM = lapply(PATH, function(x) t(read.table(x, header = TRUE, row.names = 1, col.names = c("pos", "A", "C", "G", "T")))))


  PWM_list <- motifs$PWM
  names(PWM_list) <- paste(motifs$RBP_name, "_", motifs$motif, sep = "")

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

  motif_by_kmer <- motif_by_kmer %>% dplyr::as_tibble((. > 0) + 0) %>% dplyr::mutate(motif = .$motif)
  return(motif_by_kmer)
}
fourmer_cisbpRNA <- relate_Kmer_CisBPRNA(4, "mydata/CisBPRNAdat/mm/pwms_all_motifs", "mydata/CisBPRNAdat/mm/RBP_Information_all_motifs.txt")
fivemer_cisbpRNA <- relate_Kmer_CisBPRNA(5, "mydata/CisBPRNAdat/mm/pwms_all_motifs", "mydata/CisBPRNAdat/mm/RBP_Information_all_motifs.txt")
sixmer_cisbpRNA <- relate_Kmer_CisBPRNA(6, "mydata/CisBPRNAdat/mm/pwms_all_motifs", "mydata/CisBPRNAdat/mm/RBP_Information_all_motifs.txt")
sevenmer_cisbpRNA <- relate_Kmer_CisBPRNA(7, "mydata/CisBPRNAdat/mm/pwms_all_motifs", "mydata/CisBPRNAdat/mm/RBP_Information_all_motifs.txt")

usethis::use_data(fourmer_cisbpRNA, fivemer_cisbpRNA, sixmer_cisbpRNA, sevenmer_cisbpRNA, fourmer_RBNS, fivemer_RBNS, sixmer_RBNS, sevenmer_RBNS, internal = TRUE, overwrite = TRUE)
