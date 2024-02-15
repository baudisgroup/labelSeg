#' Annotate segment profiles
#'
#' This function assigns segments SCNA labels to call different levels of SCNA.
#' It estimates calling thresholds based on logR and segment length of individual samples.
#'
#' @param data CNA segment data with logR values
#' @param genome Genome version of input data. Available options are "hg38","hg19" and "hg18". Default is "hg38".
#' @param baseshift Direction of baseline shift. Available options are  "n" (no shift), "h" (shift higher), and "l" (shift lower). Default is "n".
#' @param method Clustering method used in this algorithm. Available options are "dbscan", "hdbscan", and "optics". Default is "dbscan".
#' @param len.sep Logic value to determine whether separate clustering by segment length is applied. Default is TRUE.
#' @param minPts Integer to specify the minimum size of clusters. Default is 1.
#' @param lb.dup Lower bound for calling low-level duplication. Default is 0.15.
#' @param hb.dup Higher bound for calling low-level duplication. Default is 0.7.
#' @param lb.del Lower bound for calling low-level deletion. Default is 0.15.
#' @param hb.del Higher bound for calling low-level deletion. Default is 1.5.
#' @param shiftnum Number of target baseline clusters to shift. Only used if the parameter `baseshift` is "h" or "l".
#' @return Segment data with labels. In terms of SCNA labels, "-2" represents high-level focal deletion. "-1" represents
#' low-level deletion. "+2" represents high-level focal duplication. "+1" represents low-level duplication.
#' "0" represents no copy changes.
#' @export
#' @examples
#' data(example.seg)
#' labeled.seg <- labelseg(data=example.seg)

labelseg <- function(data, genome ='hg38',baseshift = 'n', method = 'dbscan', len.sep = TRUE, minPts = 1, lb.dup = 0.15, hb.dup = 0.7, lb.del=0.15, hb.del = 1.5, shiftnum=1){
    # check input
    genome <- match.arg(genome, c("hg38","hg19","hg18"))
    method <- match.arg(method, c("dbscan","optics","hdbscan"))
    if (len.sep == FALSE & method != 'hdbscan') stop("The `method` parameter is invalid. When `len.sep` is set to FALSE, only \"hdbscan\" method is available")
    if (minPts == 1 & method == 'hdbscan') stop("The `minPts` parameter is invalid. `minPts` should be greater than 1 when using the \"hdbscan\" method")

    data <- as.data.frame(data)
    sample.num <- length(unique(data[,1]))
    if (sample.num == 0) stop("This file is empty") 
    sample.name <- unique(data[,1])

    result <- list()
    for (i in seq_len(sample.num)){
        inddata <- data[data[,1]== sample.name[i],]
        if (len.sep){
            # separate clustering by segment length
            result[[i]] <- labelindseg(data=inddata,genome=genome,baseshift=baseshift,method=method,minPts=minPts,lb.dup=lb.dup,hb.dup=hb.dup,lb.del=lb.del,hb.del=hb.del,shiftnum=shiftnum)
        } else{
            # uniform clustering
            result[[i]] <- labelindseg_special(data=inddata,genome=genome,baseshift=baseshift,minPts=minPts,lb.dup=lb.dup,hb.dup=hb.dup,lb.del=lb.del,hb.del=hb.del,shiftnum=shiftnum)
        }
    }
    result <- do.call(rbind,result)
  
    return(result)
}


