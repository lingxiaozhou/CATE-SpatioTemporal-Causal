specify_knots <- function(n,first_point = 0.2,equal_spaced = TRUE){
  
  if(equal_spaced){
    return(1:n/(n+1))
  }else{
    res <- (1-first_point)^(1:(n+1))
    res <- cumsum(res)
    res <- res/res[n+1]
    return(res[-(n+1)])
  }
  

}
