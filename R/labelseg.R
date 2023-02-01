#' Annotate segment profiles
#'
#' This function assigns segments SCNA labels to call different levels of SCNA.
#' It estimates calling thresholds from logR distribution of individual samples.
#'
#' @param data CNA segment data with logR values
#' @param genome genome version of input data
#' @return segment data with labels
#' @export

labelseg <- function(data, genome ='hg38'){
  data <- as.data.frame(data)
  sample_num <- length(unique(data[,1]))
  if (sample_num == 0){
    cat(paste('This file is empty','\n'))
    return()
  }
  if (sample_num > 1){
    sample_name <- unique(data[,1])
    result <- list()
    for (i in seq(sample_num)){
      inddata <- data[data[,1]== sample_name[i],]
      result[[i]] <- labelindseg(data = inddata,genome=genome)
    }
    result <- do.call(rbind,result)
  } else{
    result <- labelindseg(data = data,genome=genome)
  }

  return(result)
}


labelindseg <- function(data, genome){
  set.seed(123)
  chr.bin <- get(paste0(genome,'_bins'))
  data <- as.data.frame(data)
  Start <- data[,3]
  End <- data[,4]
  chr <- data[,2]
  chr[which(chr == 'X')] <- 23
  chr[which(chr == 'Y')] <- 24
  chr <- as.numeric(chr)
  chr_len <- sapply(chr, function(x){chr.bin[x]})
  Segment_mean <- data[,6]
  len <- ceiling((End-Start)/1000000)/chr_len

  # Long segment clustering
  large.value <- Segment_mean[len >= 0.2]
  large.len <- len[len >= 0.2]

  if (length(large.value) == 0){
    cat(paste(unique(data[,1]),'is excluded because no large segments exist','\n'))
    return()
  }
  clus.res <- clus(large.value,large.len,var_thre = 0.05,eps = 0.05)

  l_db_lst <- clus.res[[1]]
  l_db_count_lst <- clus.res[[2]]
  l_db_length_lst <- sapply(l_db_count_lst, function(x){sum(x)})
  l_db_count_lst <- lapply(l_db_count_lst, function(x){x/sum(x)})
  peak_mean <- c()
  for (i in seq(length(l_db_lst))){
    peak_mean <- c(peak_mean,weighted.mean(x=l_db_lst[[i]],w=l_db_count_lst[[i]]))
  }

  l_db_length_lst <- l_db_length_lst[order(peak_mean)]
  l_db_lst <- l_db_lst[order(peak_mean)]
  peak_mean <- peak_mean[order(peak_mean)]

  ## Find baseline cluster
  neu_idx <- 0
  length_control <- 0.4
  while(length_control >= 0.2 & neu_idx == 0){
    for (i in order(abs(peak_mean))){
      if ((l_db_length_lst[i]/sum(l_db_length_lst)) >= length_control){
        if (neu_idx == 0){
          neu_idx <- i
          neu_thre <- peak_mean[i]
        }
      }
    }
    length_control <- length_control - 0.1
  }

  if (neu_idx == 0){
    neu_thre <- 0
  }

  ## Find low-level clusters
  gain_idx <- 0
  gain_thre <- 0
  gain_2_idx <- 0
  loss_idx <- 0
  loss_thre <- 0
  loss_2_idx <- 0

  for (i in c(1:length(order(-l_db_length_lst)))){
    i <- order(-l_db_length_lst)[i]
    if ((peak_mean[i]-neu_thre >= 0.15 & peak_mean[i]-neu_thre < 0.7)){
      gain_idx <- i
      break
    }
    }


  ### control sd
  if (gain_idx != 0){
    gain_thre <- peak_mean[gain_idx]
    gain_sd <- pracma::std(l_db_lst[[gain_idx]])
    if (is.na(gain_sd)){
      gain_sd <- (gain_thre-neu_thre)/10
    } else{
      if (gain_sd < (gain_thre-neu_thre)/10){
        gain_sd <-  (gain_thre-neu_thre)/10
        }
      else if (gain_sd > (gain_thre-neu_thre)/2){
        gain_sd <-  (gain_thre-neu_thre)/10
      }
    }
    }

  for (i in c(1:length(order(-l_db_length_lst)))){
    i <- order(-l_db_length_lst)[i]
    if ((neu_thre-peak_mean[i] >= 0.15 & neu_thre-peak_mean[i] < 1.5)){
      loss_idx <- i
      break
    }
    }


  if (loss_idx != 0){
    loss_thre <- peak_mean[loss_idx]
    loss_sd <- pracma::std(l_db_lst[[loss_idx]])
    if (is.na(loss_sd)){
      loss_sd <- (neu_thre-loss_thre)/10
    } else{
      if (loss_sd < (neu_thre-loss_thre)/10){
        loss_sd <-  (neu_thre-loss_thre)/10
      }
      else if (loss_sd > (neu_thre-loss_thre)/2){
        loss_sd <-  (neu_thre-loss_thre)/10
      }
    }
    }

  # Short segment clustering
  small.value <- Segment_mean[len < 0.2]
  small.len <- len[len < 0.2]

  if (length(small.value) > 0){
    clus.res <- clus(small.value,small.len,var_thre = 0.1,eps = 0.1)

    s_db_lst <- clus.res[[1]]
    s_db_count_lst <- clus.res[[2]]
    s_db_length_lst <- sapply(s_db_count_lst, function(x){sum(x)})

    focal_peak_mean <- c()
    for (i in seq(length(s_db_lst))){
      focal_peak_mean <- c(focal_peak_mean,s_db_lst[[i]][which.min(abs(s_db_lst[[i]]))])
    }
    s_db_lst <- s_db_lst[order(focal_peak_mean)]
    s_db_length_lst  <- s_db_length_lst[order(focal_peak_mean)]
    focal_peak_mean <- focal_peak_mean[order(focal_peak_mean)]
  }

  ## Find high-level clusters
  ### if low-level dup clusters exist
  if (gain_idx != 0 & length(small.value) > 0){
    gain_gap <- (gain_thre-neu_thre) * (2.2-0.6*(gain_thre-neu_thre))
    if (length(which(focal_peak_mean > gain_thre)) != 0){
      for (idx in seq(min(which(focal_peak_mean > gain_thre)),length(focal_peak_mean),1)){
        if ((focal_peak_mean[idx]-neu_thre) >= gain_gap){
          gain_2_idx <- idx
          gain_2_thre <- focal_peak_mean[idx]
          break
        }
      }
    }
  ### if low-level dup clusters are not found
  } else{
    if (length(small.value) > 0){
      for (idx in order(-s_db_length_lst)){
        if ((focal_peak_mean[idx]-neu_thre) > 0.1 & (focal_peak_mean[idx]-neu_thre) <= 0.6){
          gain_idx <- idx
          gain_thre <- focal_peak_mean[idx]
          gain_sd <- pracma::std(s_db_lst[[idx]])
          if (is.na(gain_sd)){
            gain_sd <- (gain_thre-neu_thre)/10
          } else{
            if (gain_sd < (gain_thre-neu_thre)/10){
              gain_sd <-  (gain_thre-neu_thre)/10
            }
            else if (gain_sd > (gain_thre-neu_thre)/2){
              gain_sd <-  (gain_thre-neu_thre)/10
            }
          }
          gain_gap <- (gain_thre-neu_thre) * (2.2-0.6*(gain_thre-neu_thre))
          break
          }
        }
      }
    }

  if (gain_2_idx == 0 & length(small.value) > 0){
    for (idx in seq(which.min(abs(focal_peak_mean)),length(focal_peak_mean))){
      if (gain_idx == 0){gain_gap <- 0.6}
      if ((focal_peak_mean[idx]-neu_thre) > gain_gap){
        gain_2_idx <- idx
        gain_2_thre <- focal_peak_mean[idx]
        break
      }
    }
  }

  ### if low-level del clusters exist
  if (loss_idx != 0 & length(small.value) > 0){
    loss_gap <- 2* (neu_thre-loss_thre)
    if (length(which(focal_peak_mean <  loss_thre)) != 0){
      for (idx in seq(max(which(focal_peak_mean <  loss_thre)),1,-1)){
        if ((neu_thre-focal_peak_mean[idx]) >= loss_gap){
          loss_2_idx <- idx
          loss_2_thre <- focal_peak_mean[idx]
          break
        }
      }
    }
    ### if low-level del clusters are not found
  } else{
    if (length(small.value) > 0){
    for (idx in order(-s_db_length_lst)){
      if ((neu_thre-focal_peak_mean[idx]) > 0.1 & (neu_thre-focal_peak_mean[idx]) < 1){
        loss_idx <- idx
        loss_thre <- focal_peak_mean[idx]
        loss_sd <- pracma::std(s_db_lst[[idx]])
        if (is.na(loss_sd)){
          loss_sd <- (neu_thre-loss_thre)/10
        } else{
          if (loss_sd < (neu_thre-loss_thre)/10){
            loss_sd <-  (neu_thre-loss_thre)/10
          }
          else if (loss_sd > (neu_thre-loss_thre)/2){
            loss_sd <-  (neu_thre-loss_thre)/10
          }
        }
        loss_gap <- 2* (neu_thre-loss_thre)
        break
        }
      }
  }}

  if (loss_2_idx == 0 & length(small.value) > 0){
    for (idx in seq(which.min(abs(focal_peak_mean)),1,-1)){
      if (loss_idx == 0){loss_gap <- 1}
      if ((neu_thre-focal_peak_mean[idx]) >= loss_gap){
        loss_2_idx <- idx
        loss_2_thre <- focal_peak_mean[idx]
        break
      }
    }
  }

  # only label focal change as +2
  if (gain_2_idx != 0 & gain_idx != 0){
    gain_2_thre <- ifelse(any(len >= 0.2),max(gain_2_thre,max(Segment_mean[len >= 0.2])+0.01), gain_2_thre)
  }
  if (loss_2_idx != 0 & loss_idx != 0){
    loss_2_thre <- ifelse(any(len >= 0.2),min(loss_2_thre,min(Segment_mean[len >= 0.2])-0.01), loss_2_thre)
  }


  if (loss_idx == 0){
    loss_thre <- NULL
    loss_sd <- NULL
  }
  if (loss_2_idx == 0){
    loss_2_thre <- NULL
  }
  if (gain_idx == 0){
    gain_thre <- NULL
    gain_sd <- NULL
  }
  if (gain_2_idx == 0){
    gain_2_thre <- NULL
  }

  data <- assign.label(data,gain_thre,gain_sd, gain_2_thre, loss_thre, loss_sd, loss_2_thre)


return(data)

}
