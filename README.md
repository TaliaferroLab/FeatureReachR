# FeatureReachR
The FeatureReachR packages enables users to take a deep dive into human or mouse data sets to extract biologically relevant sequence features in experimental versus control genes. FeatureReachR allows analysis of a diverse set of features with the convenience of R/RStudio. Biological insight is streamlined with this package, and prevents users from the tasks of wrangling public data sets themselves, such as RNA binding protein occupancy (e.g. CLIP-seq), conservation (phastcons), and protein binding motif databases. The package makes use of a variety of custom functions and existing tools. Our goal was not to reinvent the wheel, but rather provide a convenient platform to query sequence data in a thorough manner. See vignette(FeatureReachR) for more details.
## Installation
devtools::install_github("TaliaferroLab/FeatureReachR")
## Input Data Format
The core function accepts a list of mouse or human ensembl IDs, and returns the sequence (.fa) and annotation (.gff) for those genes. There are functions to conveniently convert gene lists to transcript lists. Users can specify if the entirety of each gene is of interest, or if only specific regions are instead wanted (whole transcript, 5'UTR, 3'UTR, CDS, or upstream promoter region). These filtered .fa and .gff files serve as the input for the downstream functions. 

It's important here to keep in mind that the type of RNAs wanted for analyses. For example, if a total RNA-seq (a.k.a. rRNA depletion) experiment is informing your gene lists of interest and the 3'UTR is selected, the core function by necessity will filter to only include mRNAs. 
## Features Queried by FeatureReachR
### RNA Sequence Content (e.g. motif, kmer, GC, length, etc)
FeatureReachR easily compares simple sequence features between experimental and control gene sets. Functions to determine the length, GC content, kmer content, or RNA or DNA binding motif presence output the raw data/counts as well as the appropriate statistics for comparison between the gene sets. ggplot wrappers use these outputs to make quick visualizations. This allows researchers to explore their data more throuroughly with less hassle.
### RNA Binding Protein Occupancy
ENCODE CLIP-seq data from the human cell lines HepG2 and K562 was wrangled to allow calculation of RBP occupancy enrichment in the experimental versus control gene sets. 
### Conservation
Average Phastcons scores from mouse and human were pre-calculated for all 5'UTRs, coding regions, and 3'UTRs. Users may extract these scores for their genes of interest, or use a custom set of genomic intervals. Due to the large size of the phastcons files, it is not recommended that users query more than 100kbp for custom calculations.

