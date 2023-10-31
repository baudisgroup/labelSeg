clus <- function(value, len, var.thre, eps, minPts, method){
    if (method == 'dbscan' | method == 'optics'){
        change.eps <- TRUE
        while(change.eps){
            if (method == 'dbscan'){
                # clustering using dbscan
                db <- dbscan::dbscan(as.matrix(value),minPts = minPts,eps = eps)
            } else {
                # clustering using optics
                db <- dbscan::optics(as.matrix(value),minPts = minPts)
                db <- dbscan::extractDBSCAN(db,eps_cl = eps)
            }
            # treat noisy points as independent clusters
            if (any(db$cluster == 0)){
                noise.pt.idx <- which(db$cluster == 0)
                db$cluster[noise.pt.idx] <- seq(max(unique(db$cluster))+1,max(unique(db$cluster))+length(noise.pt.idx))
            }
            db.lst <- list()
            db.count.lst <- list()
            for (i in unique(db$cluster)){
                db.lst[[i]] <- value[db$cluster == i]
                db.count.lst[[i]] <- len[db$cluster == i]
            }
            # self-adaptive choosing for eps by controlling sd of each cluster
            db.lst.sd <- sapply(db.lst, function(x) {if(length(x)>1){pracma::std(x)}else{0}})
            if (any(db.lst.sd >= var.thre)){
                change.eps <- TRUE
                eps <- eps-0.01
            }else{
                change.eps <- FALSE
            }
        }
          # clustering using hdbscan
    } else if (method == 'hdbscan'){
        db <- dbscan::hdbscan(as.matrix(value),minPts = minPts)
        if (any(db$cluster == 0)){
            noise.pt.idx <- which(db$cluster == 0)
            db$cluster[noise.pt.idx] <- seq(max(unique(db$cluster))+1,max(unique(db$cluster))+length(noise.pt.idx))
        }
        db.lst <- list()
        db.count.lst <- list()
        for (i in unique(db$cluster)){
            db.lst[[i]] <- value[db$cluster == i]
            db.count.lst[[i]] <- len[db$cluster == i]
        }
    }

    clus.res <- list()
    clus.res[[1]] <- db.lst
    clus.res[[2]] <- db.count.lst
    return(clus.res)
}

