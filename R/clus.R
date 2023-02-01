clus <- function(value,len,var_thre,eps){
  change_eps <- TRUE
  while(change_eps){
    db <- dbscan::dbscan(as.matrix(value),minPts = 1,eps = eps)
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
    }}
  
  clus.res <- list()
  clus.res[[1]] <- db_lst
  clus.res[[2]] <- db_count_lst
  return(clus.res)
}

