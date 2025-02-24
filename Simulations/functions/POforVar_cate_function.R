#' Calculating the estimator many times under the observed path.
#' 
#' Function that generates potential outcomes under the natural progression,
#' and calculates the estimate we would have got at every time point.
#' 
#' @param roads The time-invariant confounder that might be used in the
#' specification of the intervention.
#' @param Xt The point patterns of the time-varying confounder.
#' @param Xt2 Second time-varying confounder. Can be left NULL.
#' @param Wt The observed treatment point patterns.
#' @param Yt The observed outcome point patterns.
#' @param trt_coefs Coefficients for the treatment model intensity. The order
#' is for intercept, roads, Xt at current lag only, previous treatments, and
#' previous outcomes.
#' @param trt_lags Numeric vector of length 2 including the lags of treatments
#' and outcomes that can affect the treatment at time t.
#' @param out_coefs The true coefficients used to generate potential outcomes.
#' The order is for intercept, roads, Xt at current lag only, previous
#' treatments, and previous outcomes.
#' @param out_lags The true lags of treatment and outcome that are used to
#' generate potential outcomes.
#' @param intervention Array of dimensions: number of interventions, 2, M.
#' Entry ijk corresponds to the intercept (if j = 1) and coefficient of roads
#' (if j = 2) for the kth interventional time point of intervention i. Entry
#' k corresponds to k lags ago (for example, k = 1 corresponds to current time
#' point).
#' @param interv_cov The time-invariant covariate that can be used in the
#' specification of the intervention. If set to NULL, the first covariate from
#' roads is used. Can be left NULL if the weights are provided.
#' @param ps_coefs The coefficients used in the propensity score. If left NULL
#' the true ones are used.
#' @param ps_lags The lags of previous treatments and previous outcomes that
#' should be used for covariates based on the propensity score.
#' @param reps Number of times the simulations will be performed.
#' @param spat_kernel The spatial kernel used in estimation. Defaults to
#' 'gaussian'. Other options should be in agreement with the spatstat package.
#' @param spat_adjust Controling the bandwidth of the spatial kernel.
#' @param edge_cor Whether edge correction should be performed. Defaults to
#' FALSE.
#' @param separate_prev_trt If True, then indicate that we will previous treatment at each time separately
#' so for example W_{t-3:t} will be a vector of 4 in the outcome model
#' @param interaction a value from {0,1,2}. If 0, then no interaction. If 1, we add the term road2*prev_trt and
#' the interaction of prev_out*prev_trt
#' If 2, we add the term Xt2*prev_trt (in this case we also set separate_prev_trt = TRUE)
#' @param E.basis The covariate matrix of effect modifiers corresponding to the selected spline basis
#' @param t.basis The covatriate matrix of time corresponding to the selected spline basis
#' @param separate_model If True, assume separate model for each time period
#' 
#' @examples
#' po_for_var <- POforVar(roads = roads, Xt = Xt, Wt = Wt, Yt = Yt,
#'                        trt_coefs = trt_coefs, trt_lags = trt_lags,
#'                        out_coefs = out_coefs, out_lags = out_lags,
#'                        intervention = intervention, ps_coefs = ps_coefs,
#'                        ps_lags = ps_lags, reps = reps,
#'                        spat_kernel = spat_kernel, spat_adjust = spat_adjust,
#'                        edge_cor = edge_cor)
#' true_var <- GetTrueVar(estimates = po_for_var$estimates, B = B)
#' 
POforVar <- function(roads, Xt, Xt2 = NULL, Wt, Yt, trt_coefs, trt_lags,
                     out_coefs, out_lags, intervention, interv_cov = NULL,
                     ps_coefs = NULL, ps_lags = NULL, reps = 1000,
                     spat_kernel = 'gaussian', spat_adjust = 0.3,
                     edge_cor = FALSE, linear = FALSE,
                     interaction = 0,E.basis = E.basis,t.basis = t.basis,separate_model = TRUE,dimyx = NULL) {
  separate_prev_trt <-  (interaction==2|interaction==3)
  if (is.null(ps_coefs)) {
    ps_coefs <- trt_coefs
  }
  if (is.null(ps_lags)) {
    ps_lags <- trt_lags
  }
  if (is.null(interv_cov)) {
    interv_cov <- roads$Priority[[1]]
  }
  
  time_points <- length(Wt$Points)
  window <- Xt$Points[[1]]$window
  num_interv <- dim(intervention)[1]
  M <- dim(intervention)[3]
  N <- nrow(points.basis)
  p <- ifelse(separate_model==TRUE,ncol(E.basis)+1,(ncol(E.basis)+1)*(ncol(t.basis)+1))  # number of covariates
  npixel <- (nrow(E.basis)/time_points)^0.5  # number of pixels
  if(is.null(dimyx)){
    dimyx <- c(npixel,npixel)
  }
  
  # ----- General strategy ------ #
  # For each time point, and each of the reps:
  # (1) generate M treatments and outcomes from the observed path
  # (2) calculate the estimate
  
  progress <- floor(seq(1, time_points, length.out = 11))[2 : 11]
  
  estimates <- array(NA, dim = c(num_interv, time_points - M + 1, reps, p))
  dimnames(estimates) <- list(interv = 1 : num_interv, time = M : time_points,
                              rep = 1 : reps, beta = 1:p)
  
  weights <- array(NA, dim = c(num_interv, time_points - M + 1, reps))
  dimnames(weights) <- list(interv = 1 : num_interv, time = M : time_points,
                            rep = 1 : reps)
  

  
  for (tt in M : time_points) {
    
    if (tt %in% progress) {
      cat(tt / time_points * 100, '% completed. \n')
    }
    
    # Observed treatments and outcomes before the interventions took place.
    prev_trts <- NULL
    prev_outs <- NULL
    if (tt > M) {
      prev_trts <- Wt[1 : (tt - M)]$Points
      prev_outs <- Yt[1 : (tt - M)]$Points
    }
    
    Zt <- cbind(1,E.basis[(tt-M)*npixel**2+1:(npixel**2),])
    Zt <- na.omit(Zt)
    
    #-----------------------------------if we consider separate models for each time period-----------------------------------------------
    if(separate_model==TRUE){
      valid_time_points <- 0
      tryCatch({
        Q <- solve(t(Zt)%*%Zt)
        
        for (ss in 1 : reps) {
          
          these_outs <- prev_outs
          these_trts <- prev_trts
          
          # Generating treatments and potential outcomes for M time points:
          for (tt_ss in (tt - M + 1) : tt) {
            
            # Treatment.
            hist <- GetHistory(tt = tt_ss, Wt = these_trts, Yt = these_outs,
                               center_variable = 'W', lags = trt_lags,
                               window = window,dimyx = dimyx)
            
            this_log_lambda <-
              trt_coefs[1] + trt_coefs[2] * roads$Priority[[1]] +
              trt_coefs[3] * Xt$Priority[[tt_ss]] +
              trt_coefs[4] * hist$Wt_hist$Priority[[1]] +
              trt_coefs[5] * hist$Yt_hist$Priority[[1]]
            
            if (nrow(roads) == 2 & !is.null(Xt2)) {
              this_log_lambda <- this_log_lambda +
                trt_coefs[6] * roads$Priority[[2]] +
                trt_coefs[7] * Xt2$Priority[[tt_ss]]
            }
            
            this_lambda <- exp(this_log_lambda)
            these_trts <- append(these_trts, list(rpoispp(this_lambda, nsim = 1)))
            
            # Outcome.
            hist <- GetHistory(tt = tt_ss, Wt = these_trts, Yt = these_outs,
                               center_variable = 'Y', lags = out_lags,
                               window = window,all=!separate_prev_trt,dimyx = dimyx)
            
            
            this_log_lambda <-
              out_coefs[1] + out_coefs[2] * roads$Priority[[1]] +
              out_coefs[3] * Xt$Priority[[tt_ss]] +
              out_coefs[4] * hist$Wt_hist$Priority[[1]] +
              out_coefs[5] * hist$Yt_hist$Priority[[1]]
            
            if (nrow(roads) == 2) {
              this_log_lambda <- this_log_lambda +
                out_coefs[6] * roads$Priority[[2]] +
                out_coefs[7] * Xt2$Priority[[tt_ss]]
            }
            
            upto <- min(out_lags[1], tt_ss)
            if(separate_prev_trt==TRUE & upto>1){
              for (this_lag in 1:(upto-1)) {
                this_log_lambda <- this_log_lambda + out_coefs[this_lag+7] * hist$Wt_hist$Priority[[this_lag+1]]
              }
            }
            
            if(interaction==1){
              this_log_lambda <- this_log_lambda + 
                out_coefs[8+separate_prev_trt*(out_lags[1]-1)]* roads$Priority[[2]] * hist$Wt_hist$Priority[[1]]
              
            }
            
            if(interaction==2){
              for (this_lag in 1:upto) {
                if(tt_ss>4 | this_lag < tt_ss){
                  this_log_lambda <- this_log_lambda + 
                    out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * hist$Yt_hist$Priority[[this_lag]]
                }
              }
            }
            
            if(interaction==3){
              for (this_lag in 1:upto) {
                if(tt_ss>4 | this_lag < tt_ss){
                  this_log_lambda <- this_log_lambda + 
                    out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * Xt2$Priority[[tt_ss-this_lag+1]]
                }
              }
            }
            
            if(linear){
              this_lambda <- this_log_lambda
            }else{
              this_lambda <- exp(this_log_lambda)
            }

            these_outs <- append(these_outs, list(rpoispp(this_lambda, nsim = 1)))
          }
          
          
          est <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2,
                               Wt = hyperframe(Points = these_trts),
                               Yt = hyperframe(Points = these_outs),
                               ps_coefs = ps_coefs, ps_lags = ps_lags,
                               truncation_level = NULL,
                               intervention = intervention,
                               interv_cov = interv_cov,
                               times = tt,
                               spat_kernel = spat_kernel,
                               spat_adjust = spat_adjust, edge_cor = edge_cor)
          
          weights[, tt - M + 1, ss] <- est$weights
          for (vv in 1:num_interv) {
            tmp_est <- c(est$surfaces[[vv]][,,1])
            tmp_est <- tmp_est[!is.na(tmp_est)]
            
            estimates[vv, tt - M + 1, ss,] <- Q%*%t(Zt)%*%tmp_est
 
          }
        }
      
        valid_time_points <- valid_time_points+1
        
        }, error=function(e){
        cat("ERROR :",conditionMessage(e), "\n")
      })

    }else{
      #---------------------if we only fit one model------------------------------------
      for (ss in 1 : reps) {
        
        these_outs <- prev_outs
        these_trts <- prev_trts
        
        # Generating treatments and potential outcomes for M time points:
        for (tt_ss in (tt - M + 1) : tt) {
          
          # Treatment.
          hist <- GetHistory(tt = tt_ss, Wt = these_trts, Yt = these_outs,
                             center_variable = 'W', lags = trt_lags,
                             window = window,dimyx = dimyx)
          
          this_log_lambda <-
            trt_coefs[1] + trt_coefs[2] * roads$Priority[[1]] +
            trt_coefs[3] * Xt$Priority[[tt_ss]] +
            trt_coefs[4] * hist$Wt_hist$Priority[[1]] +
            trt_coefs[5] * hist$Yt_hist$Priority[[1]]
          
          if (nrow(roads) == 2 & !is.null(Xt2)) {
            this_log_lambda <- this_log_lambda +
              trt_coefs[6] * roads$Priority[[2]] +
              trt_coefs[7] * Xt2$Priority[[tt_ss]]
          }
          
          this_lambda <- exp(this_log_lambda)
          these_trts <- append(these_trts, list(rpoispp(this_lambda, nsim = 1)))
          
          # Outcome.
          hist <- GetHistory(tt = tt_ss, Wt = these_trts, Yt = these_outs,
                             center_variable = 'Y', lags = out_lags,
                             window = window,all=!separate_prev_trt,dimyx = dimyx)
          
          
          this_log_lambda <-
            out_coefs[1] + out_coefs[2] * roads$Priority[[1]] +
            out_coefs[3] * Xt$Priority[[tt_ss]] +
            out_coefs[4] * hist$Wt_hist$Priority[[1]] +
            out_coefs[5] * hist$Yt_hist$Priority[[1]]
          
          if (nrow(roads) == 2) {
            this_log_lambda <- this_log_lambda +
              out_coefs[6] * roads$Priority[[2]] +
              out_coefs[7] * Xt2$Priority[[tt_ss]]
          }
          
          upto <- min(out_lags[1], tt_ss)
          if(separate_prev_trt==TRUE & upto>1){
            for (this_lag in 1:(upto-1)) {
              this_log_lambda <- this_log_lambda + out_coefs[this_lag+7] * hist$Wt_hist$Priority[[this_lag+1]]
            }
          }
          
          if(interaction==1){
            this_log_lambda <- this_log_lambda + 
              out_coefs[8+separate_prev_trt*(out_lags[1]-1)]* roads$Priority[[2]] * hist$Wt_hist$Priority[[1]]
            
          }
          
          if(interaction==2){
            for (this_lag in 1:upto) {
              if(tt_ss>4 | this_lag < tt_ss){
                this_log_lambda <- this_log_lambda + 
                  out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * hist$Yt_hist$Priority[[this_lag]]
              }
            }
          }
          
          if(interaction==3){
            for (this_lag in 1:upto) {
              if(tt_ss>4 | this_lag < tt_ss){
                this_log_lambda <- this_log_lambda + 
                  out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * Xt2$Priority[[tt_ss-this_lag+1]]
              }
            }
          }
          if(linear){
            this_lambda <- this_log_lambda
          }else{
            this_lambda <- exp(this_log_lambda)
          }
          
          these_outs <- append(these_outs, list(rpoispp(this_lambda, nsim = 1)))
        }
        
        
        est <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2,
                             Wt = hyperframe(Points = these_trts),
                             Yt = hyperframe(Points = these_outs),
                             ps_coefs = ps_coefs, ps_lags = ps_lags,
                             truncation_level = NULL,
                             intervention = intervention,
                             interv_cov = interv_cov,
                             times = tt,
                             spat_kernel = spat_kernel,
                             spat_adjust = spat_adjust, edge_cor = edge_cor)
        
        weights[, tt - M + 1, ss] <- est$weights
        for (vv in 1:num_interv) {
          tmp_est <- c(est$surfaces[[vv]][,,1])
          tmp_est <- tmp_est[!is.na(tmp_est)]
          
          estimates[vv, tt - M + 1, ss,] <- t(Zt)%*%tmp_est
          
        }
      }
      
    }
    
    

  }
  
  return(list(estimates = estimates,weights = weights))
}

