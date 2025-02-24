#' Generating time varying confounders
#' 
GenConf <- function(time_points, roads, Xt_int, priority_coef = 2,dimyx = c(128,128)) {
  
  Xt <- rpoispp(lambda = exp(Xt_int[1] + Xt_int[2] * roads$Priority[[1]]),
                nsim = time_points)
  Xt <- hyperframe(Points = Xt)
  Xt <- GetPriority(Xt, priority_coef = priority_coef,dimyx = dimyx)
  
  return(Xt)
}