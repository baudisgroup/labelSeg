labelindseg <- function(data, genome, baseshift, method, minPts, lb_dup, hb_dup, lb_del, hb_del){
    set.seed(123)
    segment_mean <- as.numeric(data[,6])
    len <- get_relative_chrlen(data, genome)

    # Long segment clustering
    large.value <- segment_mean[len >= 0.2]
    large.len <- len[len >= 0.2]

    if (length(large.value) == 0){
        warning(unique(data[,1]),' is excluded because no large segments exist \n')
        return()
    }
    clus.res <- clus(large.value,large.len,var_thre = 0.05,eps = 0.05,minPts = minPts, method = method)

    l_clus_val <- clus.res[[1]]
    l_count <- clus.res[[2]]
    l_length<- sapply(l_count, function(x){sum(x)})
    l_weight <- lapply(l_count, function(x){x/sum(x)})
    peak_mean <- c()
    for (i in seq_len(length(l_clus_val))){
        peak_mean <- c(peak_mean,weighted.mean(x=l_clus_val[[i]],w=l_weight[[i]]))
    }
    ### for result stability when closed clusters have the same length
    l_clus_val <- l_clus_val[order(peak_mean)]
    l_length <- l_length[order(peak_mean)]
    peak_mean <- peak_mean[order(peak_mean)]

    ## Find baseline cluster
    neu_thre <- find_baseline(peak_mean, l_length, baseshift)
    ## Find low-level clusters
    low_dup_clus <- find_low_dup_target(l_clus_val, l_length, peak_mean, neu_thre, lb_dup, hb_dup, TRUE)
    gain_thre <- low_dup_clus[[1]]
    gain_sd <- low_dup_clus[[2]]

    low_del_clus <- find_low_del_target(l_clus_val, l_length, peak_mean, neu_thre, lb_del, hb_del, TRUE)
    loss_thre <- low_del_clus[[1]]
    loss_sd <- low_del_clus[[2]]

    # Short segment clustering
    small.value <- segment_mean[len < 0.2]
    small.len <- len[len < 0.2]

    if (length(small.value) == 0){
        data <- assign.label(data,gain_thre,gain_sd, NULL, loss_thre, loss_sd, NULL)
        return(data)
    }
   
    clus.res <- clus(small.value,small.len,var_thre = 0.1,eps = 0.1,minPts = minPts, method = method)
    s_clus_val <- clus.res[[1]]
    s_count <- clus.res[[2]]
    s_length <- sapply(s_count, function(x){sum(x)})
    focal_peak_mean <- vapply(s_clus_val, function(x){x[which.min(abs(x))]}, numeric(1))
    s_clus_val <- s_clus_val[order(focal_peak_mean)]
    s_length <- s_length[order(focal_peak_mean)]
    focal_peak_mean <- focal_peak_mean[order(focal_peak_mean)]
    ## Find low-level clusters in short clustering in case long clustering returns nothing
    if (is.null(gain_thre)){
        low_dup_clus <- find_low_dup_target(s_clus_val, s_length, focal_peak_mean, neu_thre, 0.1, 0.6, TRUE)
        gain_thre <- low_dup_clus[[1]]
        gain_sd <- low_dup_clus[[2]]
    }

    if (is.null(loss_thre)){
        low_del_clus <- find_low_del_target(s_clus_val, s_length, focal_peak_mean, neu_thre, 0.1, 1, TRUE)
        loss_thre <- low_del_clus[[1]]
        loss_sd <- low_del_clus[[2]]
    }
    
    ## Find high-level clusters

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