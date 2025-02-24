imls_to_mat <- function(imls,start = 1,end = time_points, na.matrix=NULL){
  if(!is.null(na.matrix)){
    imls <- lapply(start:end, function(x) as.matrix(imls[[x]])*na.matrix)
  }else{
    imls <- lapply(start:end, function(x) as.matrix(imls[[x]]))
  }
  imls <- array(unlist(imls),c(dim(imls[[1]]),end-start+1))
  
  return(imls)
}