# RNAreachR
The RNAreachR packages enables users to take a deep dive into their human or mouse RNA-seq data sets to extract biologically relevant features in experimental versus control genes. RNAreachR allows analysis of a diverse set of features with the convenience of R/RStudio. Biological insight is streamlined with this package, and prevents users from the tasks of wrangling public data sets themselves, such as RNA binding protein occupancy (e.g. CLIP-seq), GTEx, and RNA sequence motif databases. The package makes use of a variety of custom functions and existing tools. Our goal was not to reinvent the wheel, but rather provide a convenient platform to query RNA-seq data in a thorough manner. See vignette(RNAreachR) for more details.
## Installation
## Input Data Format
The core function accepts a list of ensembl IDs from either mouse or human, and returns the sequence (.fa) and annotation (.gff) for those genes. Users can specify if the entirety of each gene is of interest, or if only specific RNA regions are instead wanted (5'UTR, 3'UTR, or CDS). These filtered .fa and .gff files serve as the input for the downstream functions. 

It's important here to keep in mind that the type of RNA-seq libraries and sequencing that was performed influence the ability to query specific RNA regions or features. For example, if total RNA-seq (a.k.a. rRNA depletion) is performed and the 3'UTR is selected, the core function by necessity will filter to only include mRNAs. 
## Features Queried by RNAreachR
### RNA Sequence Content (e.g. motif, kmer, GC, length, etc)
### RNA Binding Protein Occupancy
ENCODE CLIP-seq data from the human cell lines HepG2 and K562 was wrangled to allow calculation of RBP occupancy enrichment in the experimental versus control gene sets. 
### Conservation
Average Phastcons scores from mouse and human were pre-calculated for all 5'UTRs, coding regions, and 3'UTRs. Users may extract these scores for their genes of interest, or use a custom set of genomic intervals. Due to the large size of the phastcons files, it is not recommended that users query more than 100kbp for custom calculations.
### RNA Structure
The average free energy, or overall stability of RNA structure, is calculated across the defined intervals using mfold.
