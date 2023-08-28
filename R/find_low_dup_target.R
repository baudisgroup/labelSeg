find_low_dup_target <- function(clus_val, clus_length, peak_mean, neu_thre, lb_dup, hb_dup, small_sd_control){
    gain_idx <- 0
    gain_thre <- NULL
    gain_sd <- NULL
    for (i in order(-clus_length)){
        if (peak_mean[i]-neu_thre >= lb_dup & peak_mean[i]-neu_thre < hb_dup){
            gain_idx <- i
            break
        }
    }

    if (gain_idx != 0){
        gain_thre <- peak_mean[gain_idx]
        gain_sd <- pracma::std(clus_val[[gain_idx]])
        # control sd
        adjusted_sd <- (gain_thre-neu_thre)/10
        if (is.na(gain_sd)){
            gain_sd <- adjusted_sd
        } else if (gain_sd > (gain_thre-neu_thre)/2){
                gain_sd <- adjusted_sd
        }

        if (small_sd_control){
            if (gain_sd < (gain_thre-neu_thre)/10){
                gain_sd <- adjusted_sd
            } 
        }
    }
      
    return(list(thre=gain_thre,sd=gain_sd))
}