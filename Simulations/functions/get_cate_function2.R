#' Calculating the conditional average treatment effects
#' This version is for using less memory!
#' 
#' Function that calculates the conditional average treatment effects: Compute the estimates for tau, using different varaince bound including standadization weights by mean weights ^(1/M), mean M_1_weights ^M and mean weights
#' 
#' @param surface1 The 3D array of weighted outcome surface according to counterfactural distribution 1. Third dimension corresponding to time points
#' @param surface2 The 3D array of weighted outcome surface according to counterfactural distribution 2. Similar to surface1
#' @param E.mat The 3D array of the covariates. Simliar to surface1
#' @param E.basis The basis matrix for the spline of covariate
#' @param base.t The number of base for time in the spline
#' @param C The standardization constant for the pixel area
#' @param points a vector of values of E.mat for prediction. 
#' @param Hajek logical. If true, then we use the Hajek version estimator
#' @param used_time a vector of used time indeces for the effect modifier. Typically, it will be 2:T-M+1
#' @param separate_model logical. If true, fit separate models for each time period (and do not include t in each model)
#' @param weights_M1 IPW weights for single time period intervention
#' @param M Lag of intervention
#' 
get_cate2 <- function(surface1, surface2,
                      E.basis,t.basis,point.basis,used_time=NULL, 
                      C=1,Hajek = TRUE,weights = NULL, separate_model=TRUE,
                      confidence_region=TRUE,perform_pred = FALSE,M=1,weights_M1=NULL){
  total_pixel <- 8451
  if(Hajek==TRUE & is.null(weights)){
    cat("no weights!\n")
  }
  if(is.null(weights_M1)){
    weights_M1 <- weights
  }
  
  test.c <- NULL
  reject <- NULL
  p.value <- NULL
  pred <- NULL
  beta_ver2 <- NULL
  beta_ver3 <- NULL
  
  npixel <- dim(surface2)[1]
  time_points <- dim(surface2)[3]
  
  # construct the basis
  if(separate_model==TRUE){
    tensor_product_basis <- cbind(1,E.basis)
  }else{
    tensor_product_basis <- cbind(1,E.basis)%.%cbind(1,t.basis)
  }
  
  if(!is.null(used_time)){
    keep_index <- ((used_time[1]-1)*total_pixel+1):(max(used_time)*total_pixel)
    tensor_product_basis <- tensor_product_basis[keep_index,]
    if(is.null(dim(tensor_product_basis))){
      tensor_product_basis <- matrix(tensor_product_basis,ncol = 1)
    }
  }
  
  
  est <- (na.omit(c(surface2))-na.omit(c(surface1)))/C
  
  # form a dataframe for the estimates and covariates
  df <- cbind.data.frame(E = tensor_product_basis,estimates1 = na.omit(c(surface1))/C,estimates2 = na.omit(c(surface2))/C,estimates = est,time = sort(rep(1:time_points,total_pixel)))
  rm(est)
  rm(tensor_product_basis)
  
  if(perform_pred){
    pred <- rep(NA,nrow(df)) 
  }
  df <- na.omit(df)
  meaneffect <- mean(df$estimates)
  
  if(separate_model == FALSE){
    df <- df[,-ncol(df)]
    # Fit the spline
    cat("fit the model\n")
    Z <- as.matrix(df[,-c(ncol(df)-0:2)])
    Q <- solve(1/time_points*t(Z)%*%Z)
    beta1 <-1/time_points*Q%*%t(Z)%*%df$estimates1
    beta2 <-1/time_points*Q%*%t(Z)%*%df$estimates2
    beta <- beta2-beta1
    if(Hajek==TRUE){
      beta <- beta2/mean(weights[2,])-beta1/mean(weights[1,])
      beta_ver2 <- beta2/(mean(weights_M1[2,])^M)-beta1/(mean(weights_M1[1,])^M)
      beta_ver3 <- beta2/(mean(weights[2,])^(1/M))-beta1/(mean(weights[1,])^(1/M))
      
    }
    
    # compute the variance and fitted values
    cat("get the asymptotic variance\n")
    if(Hajek==FALSE){
      tmp <- array(0,dim = c(time_points,ncol(Z)))
      for (i in 1:ncol(Z)) {
        tmp[,i] <- colSums(matrix((df$estimates2-df$estimates1)*Z[,i],nrow=nrow(df)/time_points))
      }
      # rm(Z)
      s2 <- t(tmp)%*%tmp/time_points
      V <- Q%*%s2%*%t(Q)/(time_points)
    }else{
      
      tmp <- array(0,dim = c(time_points,2*ncol(Z)))
      for (i in 1:(2*ncol(Z))) {
        if(i<=ncol(Z)){
          tmp[,i] <- colSums(matrix(df$estimates2*Z[,i],nrow=nrow(df)/time_points))
          
        }else{
          tmp[,i] <- colSums(matrix(df$estimates1*Z[,i-ncol(Z)],nrow=nrow(df)/time_points))
        }
      }
      # rm(Z)
      tmp <- cbind(tmp,weights[2,],weights[1,])
      k <- length(beta)
      s2 <- t(tmp)%*%tmp/time_points
      J <- cbind(diag(k),-diag(k),-t(Z)%*%df$estimates2/time_points,t(Z)%*%df$estimates1/time_points)
      V <- Q%*%J%*%s2%*%t(J)%*%t(Q)/(time_points) 
    }
    
    
    
    cat("get the predicted value\n    ") 
    V.points <- diag(length(points))
    est <- rep(0,length(points))
    
    
    for (j in 1:length(points)) {
      
      # Z.new <- matrix(rep(points.basis[j,],time_points),ncol = length(points.basis[j,]),byrow = TRUE)
      # Z.new <- Z.new %.% cbind(1,predict(t.basis,newx = 1:time_points))
      Z.new <- matrix(points.basis[j,],nrow = 1)%.%matrix(c(1,colMeans(t.basis)),nrow = 1)
      est[j] <- Z.new%*%beta
      V.points[j,j] <- Z.new%*%V%*%t(Z.new)
    }
    
    return(list(points = points, est = est, V.points = V.points,meaneffect = meaneffect))
    
  }else{
    valid_time_points <- 0
    cat("start fit the model\n")
    p <- ncol(E.basis)+1 # number of covariates
    if(Hajek==TRUE){
      beta_arr <- array(NA,c(time_points,2*p+2))
      Vt <- array(0,c(2*p+2,2*p+2))
      Vt2 <- array(0,c(2*p+2,2*p+2))
    }else{
      beta_arr <- array(NA,c(time_points,p))
      Vt <- array(0,c(p,p))
    }
    
    for (tt in 1:time_points) {
      df_tt <- df[df$time==tt,-ncol(df)]
      if(perform_pred){
        index_tt <- 1:total_pixel+(tt-1)*total_pixel
      }
      
      # Fit the spline
      tryCatch({
        Z <- as.matrix(df_tt[,-c(ncol(df_tt)-0:2)])
        Q <- solve(t(Z)%*%Z)
        if(Hajek==TRUE){
          beta_arr[tt,] <- c(Q%*%t(Z)%*%df_tt$estimates1,Q%*%t(Z)%*%df_tt$estimates2,weights[,tt])
          Vt <- Vt+beta_arr[tt,]%*%t(beta_arr[tt,])
          factor_tmp <- c(rep(mean(weights[1,]),p),rep(mean(weights[2,]),p),1,1)
          Vt2 <- Vt2+(beta_arr[tt,]/factor_tmp)%*%(t(beta_arr[tt,])/factor_tmp)
          if(perform_pred){
            pred[index_tt] <- tensor_product_basis[index_tt,]%*%(beta_arr[tt,(p+1):(2*p)]/mean(weights[2,])-beta_arr[tt,1:p]/mean(weights[1,]))
          }
        }else{
          beta_arr[tt,] <- Q%*%t(Z)%*%(df_tt$estimates2-df_tt$estimates1)
          Vt <- Vt + beta_arr[tt,]%*%t(beta_arr[tt,])
          if(perform_pred){
            pred[index_tt] <- tensor_product_basis[index_tt,]%*%beta_arr[tt,]
          }
        }
        valid_time_points <- valid_time_points+1
      }, error=function(e){
        cat("ERROR :",conditionMessage(e),"time_point is ",tt, "\n")
      })
      
    }
    Vt <- Vt/valid_time_points
    Vt2 <- Vt2/valid_time_points
    if(Hajek==TRUE){
      beta <- colMeans(beta_arr,na.rm=TRUE)
      beta1 <- beta[1:p]
      beta2 <- beta[(p+1):(2*p)]
      # Jocobian matrix 
      
      estJ <- cbind(-diag(p)/mean(weights[1,]),
                    diag(p)/mean(weights[2,]),
                    beta1/mean(weights[1,])^2,
                    -beta2/mean(weights[2,])^2)
      
      J <- cbind(-diag(p)/mean(weights[1,]),
                 diag(p)/mean(weights[2,]),
                 beta1/mean(weights[1,])^((M+1)/M),
                 -beta2/mean(weights[2,])^((M+1)/M))
      
      J2 <- cbind(-diag(p)/mean(weights[1,]),
                  diag(p)/mean(weights[2,]),
                  beta1/(mean(weights[1,])*mean(weights_M1[1,])^M),
                  -beta2/(mean(weights[2,])*mean(weights_M1[2,])^M))
      
      V_ad_hoc <- J2%*%Vt%*%t(J2)
      V <- J%*%Vt%*%t(J) # variance for beta_hat
      estV <- estJ%*%Vt%*%t(estJ) # variance for beta_hat (used in the paper)
      
      beta <- beta2/mean(weights[2,])-beta1/mean(weights[1,])  # beta_hat
      beta_ver2 <- beta2/(mean(weights_M1[2,])^M)-beta1/(mean(weights_M1[1,])^M)
      beta_ver3 <- beta2/(mean(weights[2,])^(1/M))-beta1/(mean(weights[1,])^(1/M))
      if(confidence_region){
        tryCatch({
          test.c <- t(beta[-1])%*%solve(estV[-1,-1])%*%beta[-1]*time_points
          p.value <- pchisq(test.c,df = p-1,lower.tail = FALSE)
          reject <- test.c>qchisq(0.95,df = p-1)},
          error = function(e) {
            message("Error occurred: ", conditionMessage(e))
          }
        )
        
      }
    }else{
      V <- Vt
      estV <- NULL
      V_ad_hoc <- NULL
      beta <- colMeans(beta_arr,na.rm = TRUE)
    }
    
    
    
    cat("get the predicted value and variance \n") 
    
    V.points <- points.basis%*%V%*%t(points.basis)/time_points   # NOTE we divide by T here!
    if(is.null(estV)){
      V.points.est <- NULL
    }else{
      V.points.est <- points.basis%*%estV%*%t(points.basis)/time_points
    }
    
    if(is.null(V_ad_hoc)){
      V.points.adhoc <- NULL
    }else{
      V.points.adhoc <- points.basis%*%V_ad_hoc%*%t(points.basis)/time_points
    }
    
    est <- c(points.basis%*%beta)
    est2 <- c(points.basis%*%beta_ver2)
    est3 <- c(points.basis%*%beta_ver3)
    
    
    
    return(list(points = points, est = est, est2 = est2, est3=est3, V.points = V.points,
                meaneffect = meaneffect,valid_time_points = valid_time_points, 
                beta = beta, beta_ver2 = beta_ver2,beta_ver3 = beta_ver3,
                V.points.est = V.points.est,test.c = test.c,reject = reject,estV = estV/time_points,p.value = p.value,pred=pred,beta_arr = beta_arr,estJ=estJ,V.points.adhoc = V.points.adhoc))
  }
  
  
  
}






# example
# res <- get_cate(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],E.mat,E.basis,nbaset=3,used_time=used_time, C=1, points = points,Hajek = TRUE,weights = estimates$weights[c(interv1,interv2),], separate_model=FALSE)



