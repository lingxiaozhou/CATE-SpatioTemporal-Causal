#' Calculating the estimand.
#' 
#' Integrating the intensity over the specified set B.
Integrate <- function(intensity, B) {
  
  r <- matrix(NA, nrow = length(intensity), ncol = length(B))
  for (bb in 1 : length(B)) {
    r[, bb] <- sapply(intensity, function(x) integral.im(x, domain = B[[bb]]))
  }
  colnames(r) <- names(B)
  
  return(r)
  
}