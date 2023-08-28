clus <- function(value, len, var_thre, eps, minPts, method){
    if (method == 'dbscan' | method == 'optics'){
        change_eps <- TRUE
        while(change_eps){
            if (method == 'dbscan'){
                db <- dbscan::dbscan(as.matrix(value),minPts = minPts,eps = eps)
            } else {
                db <- dbscan::optics(as.matrix(value),minPts = minPts)
                db <- dbscan::extractDBSCAN(db,eps_cl = eps)
            }
            if (any(db$cluster == 0)){
                noise_pt_idx <- which(db$cluster == 0)
                db$cluster[noise_pt_idx] <- seq(max(unique(db$cluster))+1,max(unique(db$cluster))+length(noise_pt_idx))
            }
            db_lst <- list()
            db_count_lst <- list()
            for (i in unique(db$cluster)){
                db_lst[[i]] <- value[db$cluster == i]
                db_count_lst[[i]] <- len[db$cluster == i]
            }
            db_lst_sd <- sapply(db_lst, function(x) {if(length(x)>1){pracma::std(x)}else{0}})
            if (any(db_lst_sd >= var_thre)){
                change_eps <- TRUE
                eps <- eps-0.01
            }else{
                change_eps <- FALSE
            }
        }
    } else if (method == 'hdbscan'){
        db <- dbscan::hdbscan(as.matrix(value),minPts = minPts)
        if (any(db$cluster == 0)){
            noise_pt_idx <- which(db$cluster == 0)
            db$cluster[noise_pt_idx] <- seq(max(unique(db$cluster))+1,max(unique(db$cluster))+length(noise_pt_idx))
        }
        db_lst <- list()
        db_count_lst <- list()
        for (i in unique(db$cluster)){
            db_lst[[i]] <- value[db$cluster == i]
            db_count_lst[[i]] <- len[db$cluster == i]
        }
    }

    clus.res <- list()
    clus.res[[1]] <- db_lst
    clus.res[[2]] <- db_count_lst
    return(clus.res)
}

