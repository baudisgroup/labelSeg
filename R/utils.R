get_relative_chrlen <- function(data, genome){
    data(sysdata, envir=environment())
    chr.bin <- get(paste0(genome,'_chr_bins'))
    Start <- as.numeric(data[,3])
    End <- as.numeric(data[,4])
    chr <- data[,2]
    chr[which(chr == 'X')] <- 23
    chr[which(chr == 'Y')] <- 24
    chr <- as.numeric(chr)
    chr.len <- sapply(chr, function(x){chr.bin[x]})
    # segments with len < 1MB also considered as 1MB to avoid length bias in focal cnv 
    len <- ceiling((End-Start)/1000000)/chr.len
    return(len)
}




assign.label <- function(data,gain.thre,gain.sd, gain.2.thre,loss.thre, loss.sd, loss.2.thre){
    segment.mean <- as.numeric(data[,6])
    
    data$label <- rep("0", dim(data)[1])

    if (!is.null(gain.thre)){
        data$label[segment.mean >= gain.thre-2*gain.sd] <- '+1'
    }

    if (!is.null(gain.2.thre)){
        data$label[segment.mean >= gain.2.thre] <- '+2'
    }

    if (!is.null(loss.thre)){
        data$label[segment.mean <= loss.thre+2*loss.sd] <- '-1'
    }

    if (!is.null(loss.2.thre)){
        data$label[segment.mean <= loss.2.thre] <- '-2'
    }

    return(data)
}







