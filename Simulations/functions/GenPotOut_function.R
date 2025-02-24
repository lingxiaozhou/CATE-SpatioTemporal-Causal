#' Generating potential outcomes under the intervention.
#' 
#' Function that generates potential outcomes under the stochastic intervention
#' at M time points which can be used to calculate the estimand of interest.
#' 
#' @param intervention Array of dimensions: number of interventions, 2, M.
#' Entry ijk corresponds to the intercept (if j = 1) and coefficient of roads
#' (if j = 2) for the kth interventional time point of intervention i. Entry
#' k corresponds to k lags ago (for example, k = 1 corresponds to current time
#' point).
#' @param interv_cov The time-invariant covariate that can be used in the
#' specification of the intervention. If set to NULL, the first covariate from
#' roads is used.
#' @param roads The time-invariant confounder. Can be of length 1 or 2.
#' @param Xt The point patterns of the time-varying confounder.
#' @param Xt2 Second time-varying confounder. Can be left NULL.
#' @param Wt The observed treatment point patterns.
#' @param Yt The observed outcome point patterns.
#' @param out_coefs The true coefficients used to generate potential outcomes.
#' The order is for intercept, roads, Xt at current lag only, previous
#' treatments, and previous outcomes.
#' @param out_lags The true lags of treatment and outcome that are used to
#' generate potential outcomes.
#' @param reps Number of times the simulations will be performed in order to
#' approximately calculate the estimand.
#' @param short Logical. If short is equal to TRUE, only the average surface
#' is returned. If short is equal to FALSE, the individually generated
#' potential outcomes and average surface at each time point are also returned.
#' @param separate_prev_trt If True, then indicate that we will previous treatment at each time separately
#' so for example W_{t-3:t} will be a vector of 4 in the outcome model
#' @param interaction a value from {0,1,2,3}. If 0, then no interaction. If 1, we add the term road2*prev_trt and
#' the interaction of prev_out*prev_trt
#' If 2, we add the term Xt2*prev_trt (in this case we also set separate_prev_trt = TRUE)
#' @param linear if true, the outcome intensity is a linear function. Default is false
#' 
GenPotOut <- function(intervention, interv_cov = NULL, roads, Xt, Xt2 = NULL,
                      Wt, Yt, out_coefs, out_lags, reps = 1000, short = TRUE, 
                      separate_prev_trt = FALSE, interaction = 0, linear =FALSE,dimyx = NULL) {
  
  time_points <- length(Wt$Points)
  num_interv <- dim(intervention)[1]
  M <- dim(intervention)[3]
  window <- Xt$Points[[1]]$window
  progress <- floor(seq(1, time_points, length.out = 11))[2 : 11]
  
  if (is.null(interv_cov)) {
    interv_cov <- roads$Priority[[1]]
  }

  # List where we will save all the potential outcomes.
  pixels <- dim(density(Wt$Points[[1]])$v)[1]
  if(is.null(dimyx)){
    dimyx <- c(pixels,pixels)
  }
  
  av_surface <- NULL
  for (ii in 1 : num_interv) {
    av_surface[[ii]] <- array(0, dim = c(pixels, pixels))
  }
  
  if (!short) {
    pot_out <- NULL
    surfaces <- NULL
    for (ii in 1 : num_interv) {
      pot_out[[ii]] <- list()
      surfaces[[ii]] <- list()
    }
  }
  
  # If all outcome model coefficients (except of intercept & roads) are zero,
  # we do not need to generate anything.
  
  if ((nrow(roads) == 1 & all(out_coefs[- c(1, 2)] == 0)) |
      (nrow(roads) == 2 & all(out_coefs[- c(1, 2, 5)] == 0))) {
    
    this_lambda <- out_coefs[1] + out_coefs[2] * roads$Priority[[1]]
    if (nrow(roads) == 2) {
      this_lambda <- this_lambda + out_coefs[5] * roads$Priority[[2]]
    }
    if(!linear){
      this_lambda <- exp(this_lambda) 
    }else{
      this_lambda <- this_lambda
    } 
    
    for (ii in 1 : num_interv) {
      av_surface[[ii]] <- this_lambda
      
      if (!short) {
        for (tt in M : time_points) {
          pot_out[[ii]][[tt - M + 1]] <- list(ss1 = this_lambda)
          surfaces[[ii]][[tt - M + 1]] <- this_lambda
        }
      }
    }
    
  } else {
    
    # Generating a number of treatments from which I will sample from later.
    # These treatments are from the intervention distribution.
    trts <- NULL
    for (ii in 1 : num_interv) {
      trts[[ii]] <- list()
      for (mm in 1 : M) {
        these_coefs <- intervention[ii, , mm]
        log_lam <- these_coefs[1] + these_coefs[2] * interv_cov
        trts[[ii]][[mm]] <- rpoispp(exp(log_lam), nsim = 200 * reps)
      }
    }
    
    
    # ----------------- #
    # The following code samples possible treatments for time point tt and
    # calculates the potential outcome point process intensity.
    # If the outcome process depends on itself at previous time points,
    # intermediate outcomes have to be generated.
    #
    # ----- General strategy ------ #
    # For each time point, and each of the reps:
    # (1) take M treatments
    # (2) generate potential outcome.
    # ----------------- #
    
    gen_all_pot_out <- (out_coefs[5] != 0 & M > 1)
    if (gen_all_pot_out) {
      cat('Intermediate potential outcomes will be generated. \n')
    } else {
      cat('Generating intermediate potential outcomes is avoided. \n')
    }
    
    
    for (tt in M : time_points) {
      
      if (tt %in% progress) {
        cat(tt / time_points * 100, '% completed. \n')
      }
      
      # List of potential outcomes and estimated surfaces for time point tt.
      if (!short) {
        pot_out_tt <- NULL
        for (ii in 1 : num_interv) {
          pot_out_tt[[ii]] <- list()
        }
      }

      for (ii in 1 : num_interv) {
        
        pot_out_tt_ii <- list()
        
        # Observed treatments and outcomes before the interention took place.
        prev_trts <- NULL
        prev_outs <- NULL
        if (tt > M) {
          prev_trts <- Wt[1 : (tt - M)]$Points
          prev_outs <- Yt[1 : (tt - M)]$Points
        }
        
        
        for (ss in 1 : reps) {
          
          # Treatments are the observed before tt - M and M randomly generated.
          # We need to assign treatment M first, then M - 1, up to treatment 1.
          these_trts <- prev_trts
          for (mm in 1 : M) {
            wh_trt <- sample(3 * reps, 1)
            these_trts <- c(these_trts, trts[[ii]][[M - mm + 1]][wh_trt])
          }
          names(these_trts) <- 1 : length(these_trts)

          # If we do not need to generate intermediates, we use the observed.
          these_outs <- Yt$Points[1 : tt]
          
          # Generating the intermediate outcomes if necessary.
          if (gen_all_pot_out & tt > 1) {
            
            # Outcomes at non-interventional time points:
            these_outs <- prev_outs
            
            # Outcomes for M - 1 interventional time points:
            for (tt_ss in (tt - M + 1) : (tt - 1)) {
              
              hist <- GetHistory(tt = tt_ss, Wt = these_trts, Yt = these_outs,
                                 center_variable = 'Y', lags = out_lags,
                                 window = window, all = !separate_prev_trt,dimyx = dimyx)
              
              # change the code here!
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
              
              if(linear==FALSE){
                this_out <- rpoispp(exp(this_log_lambda), nsim = 1)
              }else{
                this_out <- rpoispp(this_log_lambda, nsim = 1)
              }
              
              these_outs <- append(these_outs, list(this_out))
              
            }
          }
          
          hist <- GetHistory(tt = tt, Wt = these_trts, Yt = these_outs,
                             center_variable = 'Y', lags = out_lags,
                             window = window, all = !separate_prev_trt,dimyx = dimyx)
          
          this_log_lambda <-
            out_coefs[1] + out_coefs[2] * roads$Priority[[1]] +
            out_coefs[3] * Xt$Priority[[tt]] +
            out_coefs[4] * hist$Wt_hist$Priority[[1]] +
            out_coefs[5] * hist$Yt_hist$Priority[[1]]
          
          if (nrow(roads) == 2) {
            this_log_lambda <- this_log_lambda +
              out_coefs[6] * roads$Priority[[2]] +
              out_coefs[7] * Xt2$Priority[[tt]]
          }
          
          upto <- min(out_lags[1], tt)
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
              if(tt>4 | this_lag < tt){
                this_log_lambda <- this_log_lambda + 
                  out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * hist$Yt_hist$Priority[[this_lag]]
              }
            }
          }
          
          if(interaction==3){
            for (this_lag in 1:upto) {
              if(tt>4 | this_lag < tt){
                this_log_lambda <- this_log_lambda + 
                  out_coefs[this_lag+6+out_lags[1]]* hist$Wt_hist$Priority[[this_lag]] * Xt2$Priority[[tt-this_lag+1]]
              }
            }
          }

          # Keeping the last outcome (outcome at time point tt).
          if(!linear){
            pot_out_tt_ii[[ss]] <- exp(this_log_lambda)
          }else{
            pot_out_tt_ii[[ss]] <- this_log_lambda
          }
          
          
        }
        
        # Now we have all samples for intervention ii and time tt.
        # We can get the mean surface.
        
        mat_im <- sapply(pot_out_tt_ii, function(x) as.matrix.im(x))
        mat_im <- array(mat_im, dim = c(pixels, pixels, dim(mat_im)[2]))
        surface_tt_ii <- apply(mat_im, c(1, 2), mean)
        av_surface[[ii]] <- av_surface[[ii]] + surface_tt_ii
        
        if (!short) {
          pot_out[[ii]][[tt - M + 1]] <- pot_out_tt_ii
          surfaces[[ii]][[tt - M + 1]] <- as.im(surface_tt_ii, W = window)
        }

      }
    }
  }
  
  if ('matrix' %in% class(av_surface[[1]])) {
    for (ii in 1 : num_interv) {
      av_surface[[ii]] <- av_surface[[ii]] / (length(M : time_points))
      av_surface[[ii]] <- as.im(av_surface[[ii]], W = window)
    }    
  }

  if (short) {
    return(list(av_surface = av_surface))
  }
  
  for (ii in 1 : num_interv) {
    names(pot_out[[ii]]) <- paste('tt', M : time_points)
    names(surfaces[[ii]]) <- paste('tt', M : time_points)
  }
  return(list(surfaces = surfaces, av_surface = av_surface))
}

