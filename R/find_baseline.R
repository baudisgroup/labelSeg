find_baseline <- function(peak.mean, l.length, baseshift, shiftnum, oribase){
    # if input data is called, find the cluster closest to this calling
    if (!is.null(oribase)){
        neu.idx <- which.min(abs(peak.mean - oribase))
        neu.thre <- peak.mean[neu.idx]
    } else{
        neu.idx <- 0
        length.control <- 0.4

        while(length.control >= 0.2 & neu.idx == 0){
            for (i in order(abs(peak.mean))){
                if ((l.length[i]/sum(l.length)) >= length.control){
                    neu.idx <- i
                    neu.thre <- peak.mean[i]
                    break
                }
            } 
            length.control <- length.control - 0.1
        }

    }

    # if target baseline is not found
    if (neu.idx == 0) neu.thre <- 0

    # shift the baseline (relax size control)
    if (baseshift != "n"){
        length.control <- 0.01  # based on shortest chr length
        recheck <- 1
        candidate.idx <- NULL
        ## shift the baseline higher
        if (baseshift == 'h'){
            if (neu.idx == 0){
                candidate.idx <- which(peak.mean > 0)
            } else if (neu.idx < length(peak.mean)){
                candidate.idx <- c((neu.idx+1):length(peak.mean))
            } 
        ## shift the baseline lower
        } else if (baseshift == 'l'){
            if (neu.idx == 0){
                candidate.idx <- which(peak.mean < 0)
            } else if (neu.idx > 1){
                candidate.idx <- c((neu.idx-1):1)   
            } 
        }
        shiftnum <- min(shiftnum,length(candidate.idx))
        for (j in candidate.idx){
             if ((l.length[j]/sum(l.length)) >= length.control & recheck == shiftnum){
                neu.idx <- j
                neu.thre <- peak.mean[j]
                break    
             }
             recheck <- recheck + 1
        }
    }

    return(neu.thre)
}
