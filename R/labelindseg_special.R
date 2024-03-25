labelindseg_special <- function(data, genome, baseshift, minPts, lb.dup, hb.dup, lb.del, hb.del, shiftnum, labeled){
    set.seed(123)
    segment.mean <- as.numeric(data[,6])
    len <- get_relative_chrlen(data, genome)

    # all segment clustering
    clus.res <- clus(segment.mean,len, minPts = minPts, method = 'hdbscan')
  
    clus.val <- clus.res[[1]]
    clus.count <- clus.res[[2]]
    clus.length <- sapply(clus.count, function(x){sum(x)})
    clus.weight <- lapply(clus.count, function(x){x/sum(x)})
    ## In estimating low-level calling thresholds, feature value for each cluster is weighted mean by segment length
    peak.mean <- c()
    for (i in seq_len(length(clus.val))){
        peak.mean <- c(peak.mean,stats::weighted.mean(x=clus.val[[i]],w=clus.weight[[i]]))
    }
    ### this sorting is for result stability when closed clusters have the same length
    clus.val <- clus.val[order(peak.mean)]
    clus.length <- clus.length[order(peak.mean)]
    peak.mean <- peak.mean[order(peak.mean)]

    ## Find baseline cluster
    ### check if the data has been called
    if (labeled){
        orilabel <- data[['label']]
        oribase <- stats::weighted.mean(x=segment.mean[orilabel == '0'],w=len[orilabel == '0']/sum(len[orilabel == '0']))
    } else{
        oribase <- NULL
    }
    neu.thre <- find_baseline(peak.mean, clus.length, baseshift, shiftnum, oribase)

    ## Find low-level clusters
    low.dup.clus <- find_low_target(clus.val, clus.length, peak.mean, neu.thre, lb.dup, hb.dup, 'dup', FALSE)
    gain.thre <- low.dup.clus[[1]]
    gain.sd <- low.dup.clus[[2]]

    low.del.clus <- find_low_target(clus.val, clus.length, peak.mean, neu.thre, lb.del, hb.del, 'del', FALSE)
    loss.thre <- low.del.clus[[1]]
    loss.sd <- low.del.clus[[2]]
    
    ## Find high-level clusters
    ## In estimating high-level calling thresholds, feature value for each cluster is the logR value closest to 0
    focal.peak.mean <- vapply(clus.val, function(x){x[which.min(abs(x))]}, numeric(1))
    focal.peak.mean <- focal.peak.mean[order(focal.peak.mean)]

    gain.2.thre <-  find_high_target(focal.peak.mean, neu.thre, gain.thre, 'dup')
    loss.2.thre <-  find_high_target(focal.peak.mean, neu.thre, loss.thre, 'del')
      
    # only label focal change as high-level CNAs (+2/-2)
    if (!is.null(gain.2.thre) & !is.null(gain.thre)){
      gain.2.thre <- ifelse(any(len >= 0.2),max(gain.2.thre,max(segment.mean[len >= 0.2])+0.01), gain.2.thre)
    }
    if (!is.null(loss.2.thre) & !is.null(loss.thre)){
      loss.2.thre <- ifelse(any(len >= 0.2),min(loss.2.thre,min(segment.mean[len >= 0.2])-0.01), loss.2.thre)
    }

    data <- assign.label(data,gain.thre,gain.sd, gain.2.thre, loss.thre, loss.sd, loss.2.thre)
  
    return(data)  
}