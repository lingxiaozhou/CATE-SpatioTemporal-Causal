#' Calculating unadjusted estimator
#' 
#' Estimator that assumes that the observed treatments are assigned as
#' homogeneous poisson process with intensity equal to the avearge number of
#' observed treatment points.
CalcUnadjusted <- function(Wt, Yt, roads, Xt, intervention, interv_cov = NULL,
                           times = NULL, truncation_level = NULL,
                           spat_kernel = 'gaussian', spat_adjust = 0.5,
                           spat_sigma = 1 / 8, edge_cor = FALSE,
                           smoothed_outcome = NULL) {

  if (is.null(times)) {
    time_points <- length(Wt$Points)
    M <- dim(intervention)[3]
    times <- M : time_points
  }
  
  if (is.null(interv_cov)) {
    interv_cov <- roads$Priority[[1]]
  }
  
  window <- Wt$Points[[1]]$window
  area_window <- area(window)
  obs_num_points <- mean(with(Wt[times], Points$n))
  
  these_coefs <- c(log(obs_num_points) - log(area_window), 0, 0, 0, 0)
  
  unadjusted <- CalcEstimates(Wt = Wt, Yt = Yt, roads = roads[1, ],
                              Xt = Xt, Xt2 = NULL, ps_coefs = these_coefs,
                              ps_lags = c(1, 1),
                              intervention = intervention,
                              interv_cov = interv_cov,
                              times = times,
                              truncation_level = truncation_level,
                              spat_kernel = spat_kernel,
                              spat_adjust = spat_adjust,
                              spat_sigma = spat_sigma,
                              edge_cor = edge_cor, weights = NULL,
                              return_smoothed_outcome = FALSE,
                              smoothed_outcome = smoothed_outcome)

  unadjusted <- list(surfaces = unadjusted$surfaces,
                     av_surface = unadjusted$av_surface)
  
  return(unadjusted)
}