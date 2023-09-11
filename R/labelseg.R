#' Annotate segment profiles
#'
#' This function assigns segments SCNA labels to call different levels of SCNA.
#' It estimates calling thresholds based on logR and segment length of individual samples.
#'
#' @param data CNA segment data with logR values
#' @param genome Genome version of input data. Available options are "hg38","hg19" and "hg18". Default is "hg38".
#' @param baseshift Direction of baseline shift. Available options are  "n" (no shift), "h" (shift higher), and "l" (shift lower). Default is "n".
#' @param method Clustering method used in this algorithm. Available options are "dbscan", "hdbscan", and "optics". Default is "dbscan".
#' @param len_sep Logic value to determine whether separate clustering by segment length is applied. Default is TRUE.
#' @param minPts Integer to specify the minimum size of clusters. Default is 1.
#' @param lb_dup Lower bound for calling low-level duplication. Default is 0.15.
#' @param hb_dup Higher bound for calling low-level duplication. Default is 0.7.
#' @param lb_del Lower bound for calling low-level deletion. Default is 0.15.
#' @param hb_del Higher bound for calling low-level deletion. Default is 1.5.
#' @return Segment data with labels. In terms of SCNA labels, "-2" represents high-level deletion. "-1" represents
#' low-level deletion. "+2" represents high-level duplication. "+1" represents low-level duplication.
#' "0" represents no copy changes.
#' @export
#' @examples
#' segfile <- system.file("extdata", "example.seg", package = "labelSeg")
#' segdata <- read.table(segfile,header=T,sep='\t')
#' labeled.data <- labelseg(data=segdata)

labelseg <- function(data, genome ='hg38',baseshift = 'n', method = 'dbscan', len_sep = TRUE, minPts = 1, lb_dup = 0.15, hb_dup = 0.7, lb_del=0.15, hb_del = 1.5){
    if (!genome %in% c("hg38","hg19","hg18")) stop("Invalid genome version")
    if (len_sep  == FALSE & method != 'hdbscan') stop("The `method` parameter is invalid. When `len_sep` is set to FALSE, only \"hdbscan\" method is available")
    if (minPts == 1 & method == 'hdbscan') stop("The `minPts` parameter is invalid. `minPts` should be greater than 1 when using the \"hdbscan\" method")
    if (!method %in% c("dbscan","optics","hdbscan")) stop("The `method` parameter is invalid. (Avaible options: \"dbscan\", \"hdbscan\", \"optics\")")

    data <- as.data.frame(data)
    sample_num <- length(unique(data[,1]))
    if (sample_num == 0) stop("This file is empty") 
    sample_name <- unique(data[,1])

    result <- list()
    for (i in seq_len(sample_num)){
        inddata <- data[data[,1]== sample_name[i],]
        if (len_sep == TRUE){
            result[[i]] <- labelindseg(data=inddata,genome=genome,baseshift=baseshift,method=method,minPts=minPts,lb_dup=lb_dup,hb_dup=hb_dup,lb_del=lb_del,hb_del=hb_del)
        } else{
            result[[i]] <- labelindseg_special(data=inddata,genome=genome,baseshift=baseshift,minPts=minPts,lb_dup=lb_dup,hb_dup=hb_dup,lb_del=lb_del,hb_del=hb_del)
        }
    }
    result <- do.call(rbind,result)
  
    return(result)
}


