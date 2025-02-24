GetPriority <- function(point_patterns, priority_coef = 2, col_name = NULL,dimyx = c(128,128)) {
  
  if (!is.null(col_name)) {
    wh_col <- which(names(point_patterns) == col_name)
    names(point_patterns)[wh_col] <- 'Points'
  }
  
  all_res <- lapply(point_patterns$Points, OnePriority,
                    priority_coef = priority_coef,dimyx = dimyx)
  
  point_patterns$Dist <- lapply(all_res, function(x) x[[1]])
  point_patterns$Priority <- lapply(all_res, function(x) x[[2]])
  
  return(point_patterns)
}
