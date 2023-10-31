find_baseline <- function(peak.mean, l.length, baseshift){
    neu.idx <- 0
    length.control <- 0.4
    recheck <- 0
    while(length.control >= 0.2 & neu.idx == 0){
        for (i in order(abs(peak.mean))){
            if ((l.length[i]/sum(l.length)) >= length.control){
                if (neu.idx == 0 & recheck == 0){
                    neu.idx <- i
                    neu.thre <- peak.mean[i]
                }
                if (baseshift == 'n'){
                    break
                }else if (baseshift == 'h'){
                    neu.idx <- 0
                    recheck <- 1
                    if (peak.mean[i] > neu.thre){
                        neu.idx <- i
                        neu.thre <- peak.mean[i]
                        break
                    }
                }else if (baseshift == 'l'){
                    neu.idx <- 0
                    recheck <- 1
                    if (peak.mean[i] < neu.thre){
                        neu.idx <- i
                        neu.thre <- peak.mean[i]
                        break
                    }
                }
            } 
        }
        length.control <- length.control - 0.1
    }
    # if target baseline is not found
    if (neu.idx == 0) neu.thre <- 0
    return(neu.thre)
}
