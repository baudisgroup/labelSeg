find_high_target <- function(focal.peak.mean, neu.thre, calling.thre, type){
    high.calling.thre <- NULL
    if (!is.null(calling.thre)){
        calling.gap <- switch(type,
                              dup=(calling.thre-neu.thre) * (2.2-0.6*(calling.thre-neu.thre)),
                              del= 2* (neu.thre-calling.thre))
        
    } else{
      # in case of no low-level CNA called
        calling.gap <- switch(type,
                              dup=0.6,
                              del= 1)
        calling.thre <- neu.thre
    }
        
    if (type == 'dup'){
      if (any(focal.peak.mean > calling.thre)){
        # order in ascending order based on feature value
        for (idx in seq(min(which(focal.peak.mean > calling.thre)),length(focal.peak.mean))){
          if ((focal.peak.mean[idx]-neu.thre) >= calling.gap){
            high.calling.thre <- focal.peak.mean[idx]
            break
          }
        }
      }
    } 
    
    if (type == 'del'){
      if (any(focal.peak.mean <  calling.thre)){
        for (idx in seq(max(which(focal.peak.mean < calling.thre)),1,-1)){
          if ((neu.thre-focal.peak.mean[idx]) >= calling.gap){
            high.calling.thre <- focal.peak.mean[idx]
            break
          }
        }
      }
    }

    return(high.calling.thre)
}