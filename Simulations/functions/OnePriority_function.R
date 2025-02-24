OnePriority <- function(x, priority_coef = 2,dimyx = c(128,128)) {
  
  this_dist <- distmap(x,dimyx = dimyx)
  dims <- dim(this_dist$v)
  if (length(unique(as.numeric(this_dist$v))) == 1) {
    this_dist$v <- matrix(Inf, nrow = dims[1], ncol = dims[2])
  }
  this_priority <- exp(- priority_coef * this_dist)
  
  return(list(distance = this_dist, priority = this_priority))
  
}
