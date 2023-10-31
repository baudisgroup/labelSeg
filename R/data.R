#' A data frame consisting of the sample ID, the location of a genomic locus, number of markers, and 
#' an associated numeric value (logR) representing the relative copy change compared to the reference. 
#'
#' @return segmented data
#'
"example.seg"
#'
#' A vector that stores the number of 1MB-size genomic bins for each chromosome from the hg18 genome. 
#' It is used to normalize the segment lengths to the corresponding chromosomes.  
#'
#' @return the numebr of genomic bins in the hg18 genome
#'
"hg18_chr_bins"
#'
#' A vector that stores the number of 1MB-size genomic bins for each chromosome from the hg19 genome. 
#' It is used to normalize the segment lengths to the corresponding chromosomes.  
#'
#' @return the numebr of genomic bins in the hg19 genome
#'
"hg19_chr_bins"
#'
#' A vector that stores the number of 1MB-size genomic bins for each chromosome from the hg38 genome. 
#' It is used to normalize the segment lengths to the corresponding chromosomes.  
#'
#' @return the numebr of genomic bins in the hg38 genome
#'
"hg38_chr_bins"