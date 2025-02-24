#' 
#' @param interv_cov The time-invariant covariate that can be used in the
#' specification of the intervention. If set to NULL, the first covariate from
#' roads is used. Can be left NULL if the weights are provided.
#' 
GetLogWeights <- function(roads, Xt, Xt2 = NULL, Wt, Yt, ps_coefs, ps_lags,
                          intervention, interv_cov = NULL, times = NULL) {
  
  # Reminder: The interventional time points correspond to current time point
  # for m = 1, previous time point for m = 2, etc. So the last entry in
  # intervention corresponds to the treatment at M time points ago.
  
  time_points <- length(Xt$Points)
  window <- Xt$Points[[1]]$window
  M <- dim(intervention)[3]
  num_interv <- dim(intervention)[1]

  # The times for which the log weights will be calculated.
  if (is.null(times)) {
    times <- M : time_points
  }
  
  if (is.null(interv_cov)) {
    interv_cov <- roads$Priority[[1]]
  }
  
  interv_intens <- NULL
  for (ii in 1 : num_interv) {
    interv_intens[[ii]] <- list()
    for (mm in 1 : M) {
      this_interv <- intervention[ii, , mm]
      log_lam <- this_interv[1] + this_interv[2] * interv_cov
      interv_intens[[ii]][[mm]] <- exp(log_lam)
    }
  }
  
  int_interv <- lapply(interv_intens, function(x)
    sapply(x, function(y) integral(y, domain = window)))
  
  log_weights <- array(NA, dim = c(num_interv, length(times), M))
  dimnames(log_weights) <- list(interv = 1 : num_interv,
                                time = times, prev_lags_mm = 1 : M)
  
  
  for (tt in 1 : length(times)) {
    
    this_time <- times[tt]
    
    for (mm in 1 : M) {
      
      # Treatment at mm - 1 time points ago:
      this_trt <- Wt$Points[[this_time - mm + 1]]
      
      # History for treatment mm - 1 steps before.
      hist <- GetHistory(tt = this_time - (mm - 1), Wt = Wt$Points,
                         Yt = Yt$Points, center_variable = 'W', lags = ps_lags,
                         window = window)
      
      ps_log_intens <- ps_coefs[1] + ps_coefs[2] * roads$Priority[[1]] +
                         ps_coefs[3] * Xt$Priority[[this_time - mm + 1]] +
                         ps_coefs[4] * hist$Wt_hist$Priority[[1]] +
                         ps_coefs[5] * hist$Yt_hist$Priority[[1]]
      
      if (nrow(roads) == 2) {
        ps_log_intens <- ps_log_intens + ps_coefs[6] * roads$Priority[[2]] +
          ps_coefs[7] * Xt2$Priority[[this_time - mm + 1]]
      }
      
      ps_intens <- exp(ps_log_intens)
      
      int_ps_intens <- integral(ps_intens, domain = window)
      ps_at_points <- interp.im(ps_intens, x = this_trt)
      
      # If interpolation did not work, shifting points sligtly towards center
      w <- 0.999
      window_center <- centroid.owin(window)
      while(any(is.na(ps_at_points))) {
        wh_na <- which(is.na(ps_at_points))
        x_locs <- w * this_trt$x[wh_na] + (1 - w) * window_center$x
        y_locs <- w * this_trt$y[wh_na] + (1 - w) * window_center$y
        ps_at_points[wh_na] <- interp.im(ps_intens, x = x_locs, y = y_locs)
        w <- w - 0.001
      }
      
      for (ii in 1 : num_interv) {
      
        # Intensity at mm - 1 time points ago for intervention ii:
        interv_at_points <- interp.im(interv_intens[[ii]][[mm]], x = this_trt)
        
        # If interpolation did not work, shifting points sligtly towards center
        w <- 0.999
        while(any(is.na(interv_at_points))) {
          wh_na <- which(is.na(interv_at_points))
          x_locs <- w * this_trt$x[wh_na] + (1 - w) * window_center$x
          y_locs <- w * this_trt$y[wh_na] + (1 - w) * window_center$y
          interv_at_points[wh_na] <- interp.im(interv_intens[[ii]][[mm]],
                                               x = x_locs, y = y_locs)
          w <- w - 0.001
        }
        
        lw <- - int_interv[[ii]][mm] + int_ps_intens
        lw <- lw + sum(log(interv_at_points)) - sum(log(ps_at_points))
        log_weights[ii, tt, mm] <- lw  
        
      }
    }
  }
  
  return(log_weights)
}