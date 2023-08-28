find_high_dup_target <- function(focal_peak_mean, neu_thre, gain_thre){
    gain_2_thre <- NULL
    if (!is.null(gain_thre)){
        gain_gap <- (gain_thre-neu_thre) * (2.2-0.6*(gain_thre-neu_thre))
        if (any(focal_peak_mean > gain_thre)){
            for (idx in seq(min(which(focal_peak_mean > gain_thre)),length(focal_peak_mean))){
                if ((focal_peak_mean[idx]-neu_thre) >= gain_gap){
                    gain_2_thre <- focal_peak_mean[idx]
                    break
                }
            }
        }
    } else{
        gain_gap <- 0.6
        if (any(focal_peak_mean > neu_thre)){
            for (idx in seq(min(which(focal_peak_mean > neu_thre)),length(focal_peak_mean))){      
                if ((focal_peak_mean[idx]-neu_thre) >= gain_gap){               
                    gain_2_thre <- focal_peak_mean[idx]
                    break
                }
            }
        }
    }
    return(gain_2_thre)
}