labelindseg_special <- function(data, genome, baseshift, minPts, lb_dup, hb_dup, lb_del, hb_del){
    set.seed(123)
    segment_mean <- as.numeric(data[,6])
    len <- get_relative_chrlen(data, genome)

    # all segment clustering
    clus.res <- clus(segment_mean,len, minPts = minPts, method = 'hdbscan')
  
    clus_val <- clus.res[[1]]
    clus_count <- clus.res[[2]]
    clus_length <- sapply(clus_count, function(x){sum(x)})
    clus_weight <- lapply(clus_count, function(x){x/sum(x)})
    peak_mean <- c()
    for (i in seq_len(length(clus_val))){
        peak_mean <- c(peak_mean,weighted.mean(x=clus_val[[i]],w=clus_weight[[i]]))
    }
    ### for result stability when closed clusters have the same length
    clus_val <- clus_val[order(peak_mean)]
    clus_length <- clus_length[order(peak_mean)]
    peak_mean <- peak_mean[order(peak_mean)]

    ## Find baseline cluster
    neu_thre <- find_baseline(peak_mean, clus_length, baseshift)

    ## Find low-level clusters
    low_dup_clus <- find_low_dup_target(clus_val, clus_length, peak_mean, neu_thre, lb_dup, hb_dup, FALSE)
    gain_thre <- low_dup_clus[[1]]
    gain_sd <- low_dup_clus[[2]]

    low_del_clus <- find_low_del_target(clus_val, clus_length, peak_mean, neu_thre, lb_del, hb_del, FALSE)
    loss_thre <- low_del_clus[[1]]
    loss_sd <- low_del_clus[[2]]
    
    ## Find high-level clusters
    focal_peak_mean <- vapply(clus_val, function(x){x[which.min(abs(x))]}, numeric(1))
    focal_peak_mean <- focal_peak_mean[order(focal_peak_mean)]

    gain_2_thre <-  find_high_dup_target(focal_peak_mean, neu_thre, gain_thre)
    loss_2_thre <-  find_high_del_target(focal_peak_mean, neu_thre, loss_thre)
      
    # only label focal change as +2
    if (!is.null(gain_2_thre) & !is.null(gain_thre)){
      gain_2_thre <- ifelse(any(len >= 0.2),max(gain_2_thre,max(segment_mean[len >= 0.2])+0.01), gain_2_thre)
    }
    if (!is.null(loss_2_thre) & !is.null(loss_thre)){
      loss_2_thre <- ifelse(any(len >= 0.2),min(loss_2_thre,min(segment_mean[len >= 0.2])-0.01), loss_2_thre)
    }

    data <- assign.label(data,gain_thre,gain_sd, gain_2_thre, loss_thre, loss_sd, loss_2_thre)
  
    return(data)  
}