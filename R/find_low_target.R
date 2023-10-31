find_low_target <- function(clus.val, clus.length, peak.mean, neu.thre, lb, hb, type, small.sd.control){
    calling.idx <- 0
    calling.thre <- NULL
    calling.sd <- NULL
    # order in descending order based on the cumulative segment length
    for (i in order(-clus.length)){
          clus.diff <- switch(type,
                              dup=peak.mean[i]-neu.thre,
                              del=neu.thre-peak.mean[i])
          
          if (clus.diff >= lb & clus.diff < hb){
            calling.idx <- i
            break
        }
    }
    # sd adjustment
    if (calling.idx != 0){
        calling.thre <- peak.mean[calling.idx]
        calling.sd <- pracma::std(clus.val[[calling.idx]])
        adjusted.sd <- abs(calling.thre-neu.thre)/10
        if (is.na(calling.sd)){
          calling.sd <- adjusted.sd 
        } else if (calling.sd > abs(calling.thre-neu.thre)/2){
          calling.sd <- adjusted.sd 
        } ## for optics and dbscan where variance control is applied
        if (small.sd.control){
            if (calling.sd < abs(calling.thre-neu.thre)/10){
                calling.sd <- adjusted.sd 
            }
        }
    }
    return(list(thre=calling.thre,sd=calling.sd))    
}
