find_baseline <- function(peak_mean, l_length, baseshift){
    neu_idx <- 0
    length_control <- 0.4
    recheck <- 0
    while(length_control >= 0.2 & neu_idx == 0){
        for (i in order(abs(peak_mean))){
            if ((l_length[i]/sum(l_length)) >= length_control){
                if (neu_idx == 0 & recheck == 0){
                    neu_idx <- i
                    neu_thre <- peak_mean[i]
                }
                if (baseshift == 'n'){
                    break
                }else if (baseshift == 'h'){
                    neu_idx <- 0
                    recheck <- 1
                    if (peak_mean[i] > neu_thre){
                        neu_idx <- i
                        neu_thre <- peak_mean[i]
                        break
                    }
                }else if (baseshift == 'l'){
                    neu_idx <- 0
                    recheck <- 1
                    if (peak_mean[i] < neu_thre){
                        neu_idx <- i
                        neu_thre <- peak_mean[i]
                        break
                    }
                }
            } 
        }
        length_control <- length_control - 0.1
    }
    # if target baseline is not found
    if (neu_idx == 0) neu_thre <- 0
    return(neu_thre)
}
