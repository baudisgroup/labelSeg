get_relative_chrlen <- function(data, genome){
    chr.bin <- get(paste0(genome,'_chr_bins'))
    Start <- as.numeric(data[,3])
    End <- as.numeric(data[,4])
    chr <- data[,2]
    chr[which(chr == 'X')] <- 23
    chr[which(chr == 'Y')] <- 24
    chr <- as.numeric(chr)
    chr_len <- sapply(chr, function(x){chr.bin[x]})
    # segments with len < 1MB also considered as 1MB to avoid length bias in focal cnv 
    len <- ceiling((End-Start)/1000000)/chr_len
    return(len)
}




assign.label <- function(data,gain_thre,gain_sd, gain_2_thre,loss_thre, loss_sd, loss_2_thre){
    segment_mean <- as.numeric(data[,6])
    
    data$label <- rep("0", dim(data)[1])

    if (!is.null(gain_thre)){
        data$label[segment_mean >= gain_thre-2*gain_sd] <- '+1'
    }

    if (!is.null(gain_2_thre)){
        data$label[segment_mean >= gain_2_thre] <- '+2'
    }

    if (!is.null(loss_thre)){
        data$label[segment_mean <= loss_thre+2*loss_sd] <- '-1'
    }

    if (!is.null(loss_2_thre)){
        data$label[segment_mean <= loss_2_thre] <- '-2'
    }

    return(data)
}







