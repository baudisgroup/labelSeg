find_high_del_target <- function(focal_peak_mean, neu_thre, loss_thre){
    loss_2_thre <- NULL
    if (!is.null(loss_thre)){
        loss_gap <- 2* (neu_thre-loss_thre)
        if (any(focal_peak_mean <  loss_thre)){
            for (idx in seq(max(which(focal_peak_mean < loss_thre)),1,-1)){
                if ((neu_thre-focal_peak_mean[idx]) >= loss_gap){
                    loss_2_thre <- focal_peak_mean[idx]
                    break
                }
            }
        }
    } else{
        loss_gap <- 1
        if (any(focal_peak_mean < neu_thre)){
            for (idx in seq(max(which(focal_peak_mean < neu_thre)),1,-1)){
                if ((neu_thre-focal_peak_mean[idx]) >= loss_gap){
                    loss_2_thre <- focal_peak_mean[idx]
                    break
                }
            }
        }
    }
    return(loss_2_thre)
}