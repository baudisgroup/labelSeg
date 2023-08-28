find_low_del_target <- function(clus_val, clus_length, peak_mean, neu_thre, lb_del, hb_del, small_sd_control){
    loss_idx <- 0
    loss_thre <- NULL
    loss_sd <- NULL
    for (i in order(-clus_length)){
        if (neu_thre-peak_mean[i] >= lb_del & neu_thre-peak_mean[i] < hb_del){
            loss_idx <- i
            break
        }
    }

    if (loss_idx != 0){
        loss_thre <- peak_mean[loss_idx]
        loss_sd <- pracma::std(clus_val[[loss_idx]])
        adjusted_sd <- (neu_thre-loss_thre)/10
        if (is.na(loss_sd)){
            loss_sd <- adjusted_sd 
        } else if (loss_sd > (neu_thre-loss_thre)/2){
            loss_sd <- adjusted_sd 
        }
        if (small_sd_control){
            if (loss_sd < (neu_thre-loss_thre)/10){
                loss_sd <- adjusted_sd 
            }
        }
    }
    return(list(thre=loss_thre,sd=loss_sd))    
}
