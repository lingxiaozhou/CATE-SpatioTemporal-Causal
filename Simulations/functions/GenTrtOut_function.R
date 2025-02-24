#' Function that generates treatment and outcome point patterns.
#' 
#' @param roads Image confounder that does not vary over time. There can be one
#' or two of them.
#' @param Xt Time varying confounder.
#' @param Xt2 Second time varying confounder. Can be left NULL.
#' @param trt_coefs Coefficients for the treatment model intensity. The order
#' is for intercept, roads, Xt at current lag only, previous treatments, and
#' previous outcomes (additional coefficient for prev_trt and interaction terms)
#' @param trt_lags Numeric vector of length 2 including the lags of treatments
#' and outcomes that can affect the treatment at time t.
#' @param out_coefs Coefficients in the outcome intensity. The order of the
#' coefficients is the same as in trt_coefs.
#' @param out_lags Numeric vector of length 2 including the lags of treatments
#' and outcomes that can affect the outcome at time t.
#' @param save_hist If True, then save the hist covariates. Default is false
#' @param separate_prev_trt If True, then indicate that we will previous treatment at each time separately
#' so for example W_{t-3:t} will be a vector of 4 in the outcome model
#' @param interaction a value from {0,1,2}. If 0, then no interaction. If 1, we add the term road2*prev_trt
#' If 2, we add the term prev_out*prev_trt (in this case we also set separate_prev_trt = TRUE)
#' @param linear if true, the intensity of outcome is a linear function. default is FALSE 
#' 
GenTrtOut <- function(roads, Xt, Xt2 = NULL, trt_coefs, trt_lags, out_coefs,
                      out_lags, save_hist = FALSE, separate_prev_trt = FALSE, interaction=0, linear = FALSE,dimyx = c(128,128)) {
  if(interaction==2|interaction==3){
    separate_prev_trt <- TRUE
  }
  
  correct_num_coefs_trt <- 4 + nrow(roads) + !is.null(Xt2)
  correct_num_coefs_out <- (4 + nrow(roads) + !is.null(Xt2)) + as.numeric(interaction==1)+
    as.numeric(separate_prev_trt)*(out_lags[1]-1) + 
    as.numeric(interaction==2)*out_lags[1] +
    as.numeric(interaction==3)*out_lags[1] 
    
  if(length(trt_coefs) != correct_num_coefs_trt | length(out_coefs) != correct_num_coefs_out){
    stop('Wrong number of coefficients.')
  }
  
  time_points <- length(Xt$Points)
  window <- Xt$Points[[1]]$window
  
  Wt <- NULL
  Yt <- NULL
  
  for (tt in 1 : time_points) {
    
    # Creating history of treatments and outcomes for generating treatment.
    hist <- GetHistory(tt = tt, Wt = Wt, Yt = Yt, center_variable = 'W',
                       lags = trt_lags, window = window,dimyx = dimyx)

    # Generating the treatment.
    
    mean_struc <-
      trt_coefs[1] + trt_coefs[2] * roads$Priority[[1]] +
      trt_coefs[3] * Xt$Priority[[tt]] +
      trt_coefs[4] * hist$Wt_hist$Priority[[1]] +
      trt_coefs[5] * hist$Yt_hist$Priority[[1]]
    
    if (nrow(roads) == 2 & !is.null(Xt2)) {
      mean_struc <- mean_struc +
        trt_coefs[6] * roads$Priority[[2]] +
        trt_coefs[7] * Xt2$Priority[[tt]]
    }

    
    Wt[[tt]] <- rpoispp(exp(mean_struc), nsim = 1)
    
    
    # Creating history of treatments and outcomes for generating the outcome.
    hist <- GetHistory(tt = tt, Wt = Wt, Yt = Yt, center_variable = 'Y',
                       lags = out_lags, window = window, all = !separate_prev_trt,dimyx = dimyx)
    
    # Generating the outcome.
    
    mean_struc <-
      out_coefs[1] + out_coefs[2] * roads$Priority[[1]] +
      out_coefs[3] * Xt$Priority[[tt]] +
      out_coefs[4] * hist$Wt_hist$Priority[[1]] +
      out_coefs[5] * hist$Yt_hist$Priority[[1]]
    
    if (nrow(roads) == 2 & !is.null(Xt2)) {
      mean_struc <- mean_struc +
        out_coefs[6] * roads$Priority[[2]] +
        out_coefs[7] * Xt2$Priority[[tt]]
    }
    
    upto <- min(out_lags[1], tt)
    if(separate_prev_trt==TRUE & upto!=1){
      for (this_lag in 1:(upto-1)) {
        mean_struc <- mean_struc + out_coefs[this_lag+7] * hist$Wt_hist$Priority[[this_lag+1]]
      }
    }
    
    if(interaction==1){
        mean_struc <- mean_struc + 
          out_coefs[8+separate_prev_trt*(out_lags[1]-1)]* roads$Priority[[2]] * hist$Wt_hist$Priority[[1]]
    }
    
    if(interaction==2){
      for (this_lag in 1:upto) {
        if(tt>4 | this_lag < tt){
          mean_struc <- mean_struc + 
            out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * hist$Yt_hist$Priority[[this_lag]]
        }
        
      }
    }
    
    if(interaction==3){
      for (this_lag in 1:upto) {
        if(tt>4 | this_lag < tt){
          mean_struc <- mean_struc + 
            out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * Xt2$Priority[[tt-this_lag+1]]
        }
        
      }
    }
    
    
    if(linear==FALSE){
      Yt[[tt]] <- rpoispp(exp(mean_struc), nsim = 1)
    }else{
      Yt[[tt]] <- rpoispp(mean_struc, nsim = 1)
    }
    
    
  }
  
  return(list(Wt = hyperframe(Points = Wt),
              Yt = hyperframe(Points = Yt)))
}