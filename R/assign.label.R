assign.label <- function(data,gain_thre,gain_sd, gain_2_thre,loss_thre, loss_sd, loss_2_thre){
  Segment_mean <- data[,6]
  data$label <- rep(NA, dim(data)[1])
  
  if (!(is.null(gain_thre)) & !(is.null(loss_thre))){
    data$label[Segment_mean < gain_thre-2*gain_sd & Segment_mean > loss_thre+2*loss_sd] <- '0'
    if (!(is.null(gain_2_thre))){
      data$label[Segment_mean >= gain_thre-2*gain_sd & Segment_mean < gain_2_thre] <- '+1'
      data$label[Segment_mean >= gain_2_thre] <- '+2'
    } else{
      data$label[Segment_mean >= gain_thre-2*gain_sd] <- '+1'
    }
    
    if (!(is.null(loss_2_thre))){
      data$label[Segment_mean <= loss_thre+2*loss_sd & Segment_mean > loss_2_thre] <- '-1'
      data$label[Segment_mean <= loss_2_thre] <- '-2'
    } else{
      data$label[Segment_mean <= loss_thre+2*loss_sd] <- '-1'
    }
    
    
    
  } else if (!(is.null(gain_thre)) & is.null(loss_thre)){
    
    if (!(is.null(gain_2_thre))){
      data$label[Segment_mean >= gain_thre-2*gain_sd & Segment_mean < gain_2_thre] <- '+1'
      data$label[Segment_mean >= gain_2_thre] <- '+2'
    } else{
      data$label[Segment_mean >= gain_thre-2*gain_sd] <- '+1'
    }
    
    if (!(is.null(loss_2_thre))){
      data$label[Segment_mean < gain_thre-2*gain_sd & Segment_mean > loss_2_thre] <- '0'
      data$label[Segment_mean <= loss_2_thre] <- '-2'
    } else{
      data$label[Segment_mean < gain_thre-2*gain_sd] <- '0'
    }
    
    
    
  } else if (is.null(gain_thre) & !(is.null(loss_thre))){
    
    if (!(is.null(gain_2_thre))){
      data$label[Segment_mean < gain_2_thre & Segment_mean > loss_thre+2*loss_sd] <- '0'
      data$label[Segment_mean >= gain_2_thre] <- '+2'
    } else{
      data$label[Segment_mean > loss_thre+2*loss_sd] <- '0'
    }
    
    if (!(is.null(loss_2_thre))){
      data$label[Segment_mean <= loss_thre+2*loss_sd & Segment_mean > loss_2_thre] <- '-1'
      data$label[Segment_mean <= loss_2_thre] <- '-2'
    } else{
      data$label[Segment_mean <= loss_thre+2*loss_sd] <- '-1'
    }
    
    
  } else if (is.null(gain_thre) & is.null(loss_thre)){
    if (!(is.null(gain_2_thre))){
      data$label[Segment_mean >= gain_2_thre] <- '+2'
      if (!(is.null(loss_2_thre))){
        data$label[Segment_mean <= loss_2_thre] <- '-2'
        data$label[Segment_mean < gain_2_thre & Segment_mean > loss_2_thre] <- '0'
      } else{
        data$label[Segment_mean < gain_2_thre] <- '0'
      }
    } else if (!(is.null(loss_2_thre))){
      data$label[Segment_mean > loss_2_thre] <- '0'
      data$label[Segment_mean <= loss_2_thre] <- '-2'
    } else{
      data$label <- '0'
    }}
  return(data)
}


