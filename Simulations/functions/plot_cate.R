plot_cate <- function(est,V.points,points,ylim = NULL, main="", xlab="Xt", ylab = "CATE", ratio = FALSE,meaneffect = NULL,V.points.ratio = NULL,xrange = NULL){

  
  if(ratio == FALSE){
    m <- est
    V <- diag(V.points)
  }else{
    m <- est/abs(meaneffect)
    V <- diag(V.points.ratio)
  }
  
  ub <- m+qnorm(c(0.975))*sqrt(V)
  lb <- m+qnorm(c(0.025))*sqrt(V)
  df <- cbind.data.frame(ub,lb,m,x = points)
  
  if(is.null(xrange) == FALSE){
   df <- df[df$x >= xrange[1],]
   df <- df[df$x <= xrange[2],] 
  }
  
  b <- max(df$ub)
  a <- min(df$lb)
  
  if(is.null(ylim)==FALSE){
    b <- max(b,ylim[2])
    a <- min(a,ylim[1])
  }
  
  main <- ifelse(ratio==TRUE, paste0(main," total effect =", round(meaneffect*8455,2)),main)
  ylab <- ifelse(ratio==TRUE,"Ratio",'CATE')
  
  
  plot(df$x, df$m,ylim = c(a,b),type = "l",xlab = xlab, ylab = ylab,main = main)
  #make polygon where coordinates start with lower limit and
  # then upper limit in reverse order
  polygon(c(df$x,rev(df$x)),c(df$lb,rev(df$ub)),col = "grey75", border = FALSE)
  lines(df$x, df$m, lwd = 2, col = "black")
  abline(h=0,col = "red",lty = "dashed")
  if(ratio==TRUE){
    abline(h=meaneffect,col = "blue",lty = "dashed")
  }


}
