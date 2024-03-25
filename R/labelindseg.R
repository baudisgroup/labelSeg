labelindseg <- function(data, genome, baseshift, method, minPts, lb.dup, hb.dup, lb.del, hb.del, shiftnum, labeled){
    set.seed(123)
    segment.mean <- as.numeric(data[,6])
    len <- get_relative_chrlen(data, genome)

    # Long segment clustering
    large.value <- segment.mean[len >= 0.2]
    large.len <- len[len >= 0.2]

    if (length(large.value) == 0){
        warning(unique(data[,1]),' is excluded because no large segments exist \n')
        return()
    }
    clus.res <- clus(large.value,large.len,var.thre = 0.05,eps = 0.05,minPts = minPts, method = method)
    ## logR value in each cluster
    l.clus.val <- clus.res[[1]]
    ## segment length in each cluster
    l.count <- clus.res[[2]]
    l.length<- sapply(l.count, function(x){sum(x)})
    l.weight <- lapply(l.count, function(x){x/sum(x)})
    ## feature value in long segment clusters is weighted mean by segment length
    peak.mean <- c()
    for (i in seq_len(length(l.clus.val))){
        peak.mean <- c(peak.mean,stats::weighted.mean(x=l.clus.val[[i]],w=l.weight[[i]]))
    }
    ### this sorting is for result stability when closed clusters have the same length
    l.clus.val <- l.clus.val[order(peak.mean)]
    l.length <- l.length[order(peak.mean)]
    peak.mean <- peak.mean[order(peak.mean)]

    ## Find baseline cluster
    ### check if the data has been called
    if (labeled){
        orilabel <- data[['label']]
        oribase <- stats::weighted.mean(x=segment.mean[orilabel == '0'],w=len[orilabel == '0']/sum(len[orilabel == '0']))
    } else{
        oribase <- NULL
    }
    neu.thre <- find_baseline(peak.mean, l.length, baseshift, shiftnum, oribase)
    ## Find low-level target clusters
    ### for optics and dbscan where variance control is applied
    small.sd.control <- TRUE
    if (method == "hdbscan") small.sd.control <- FALSE
    
    low.dup.clus <- find_low_target(l.clus.val, l.length, peak.mean, neu.thre, lb.dup, hb.dup, 'dup', small.sd.control)
    gain.thre <- low.dup.clus[[1]]
    gain.sd <- low.dup.clus[[2]]

    low.del.clus <- find_low_target(l.clus.val, l.length, peak.mean, neu.thre, lb.del, hb.del, 'del', small.sd.control)
    loss.thre <- low.del.clus[[1]]
    loss.sd <- low.del.clus[[2]]

    # Short segment clustering
    small.value <- segment.mean[len < 0.2]
    small.len <- len[len < 0.2]

    if (length(small.value) == 0){
        data <- assign.label(data,gain.thre,gain.sd, NULL, loss.thre, loss.sd, NULL)
        return(data)
    }
   
    clus.res <- clus(small.value,small.len,var.thre = 0.1,eps = 0.1,minPts = minPts, method = method)
    s.clus.val <- clus.res[[1]]
    s.count <- clus.res[[2]]
    s.length <- sapply(s.count, function(x){sum(x)})
    ## feature value in short segment clusters is the logR value closest to 0 
    focal.peak.mean <- vapply(s.clus.val, function(x){x[which.min(abs(x))]}, numeric(1))
    s.clus.val <- s.clus.val[order(focal.peak.mean)]
    s.length <- s.length[order(focal.peak.mean)]
    focal.peak.mean <- focal.peak.mean[order(focal.peak.mean)]
    ## Find low-level clusters in short clustering in case long clustering returns nothing
    if (is.null(gain.thre)){
        low.dup.clus <- find_low_target(s.clus.val, s.length, focal.peak.mean, neu.thre, 0.1, 0.6, 'dup', small.sd.control)
        gain.thre <- low.dup.clus[[1]]
        gain.sd <- low.dup.clus[[2]]
    }

    if (is.null(loss.thre)){
        low.del.clus <- find_low_target(s.clus.val, s.length, focal.peak.mean, neu.thre, 0.1, 1, 'del', small.sd.control)
        loss.thre <- low.del.clus[[1]]
        loss.sd <- low.del.clus[[2]]
    }
    
    ## Find high-level target clusters

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