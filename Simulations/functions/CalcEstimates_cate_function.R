#' Calculating the weighted surfaces.
#' 
#' Function that calculates the weighted surface at every time point, and
#' their average over time.
#' 
#' @param Wt The list of observed treatment point patterns.
#' @param Yt The list of observed outcome point patterns.
#' @param roads The hyperframe of time invariant covariates. Has to include
#' a column named Priority that will be used. Can have 1 or 2 rows. Can be
#' left NULL if the weights are provided.
#' @param Xt The first time varying covariate. Needs to include a column
#' that is named Priority, and this column will be used. Can be left NULL if
#' the weights are provided.
#' @param Xt2 Second time varying covariate. Can be left NULL if roads has
#' one row only. Can be left NULL if the weights are provided.
#' @param ps_coefs The coefficients in the propensity score model. Should be
#' of length 5 or 7 depending on whether we have 2 or 4 covariates. The order
#' of the coefficients is intercept, first row of roads, Xt, previous
#' treatment, previous outcome, second road, and Xt2. Can be left NULL if the
#' weights are provided.
#' @param ps_lags Vector of length 2. Specifies how many previous lags of
#' treatment and outcome are used in the propensity score. Can be left NULL if
#' the weights are provided.
#' @param intervention Array of dimension 3. The first dimension is the
#' number of unique interventions we consider. The second dimension includes
#' two values, the intercept and coefficient of covariate in the intervention,
#' the third covariate correspond to the lag. Interventions can vary across
#' lags.
#' @param interv_cov The time-invariant covariate that can be used in the
#' specification of the intervention. If set to NULL, the first covariate from
#' roads is used. Can be left NULL if the weights are provided.
#' @param times
#' @param truncation_level
#' @param spat_kernel
#' @param spat_adjust
#' @param spat_sigma
#' @param edge_cor
#' @param weights
#' @param return_smoothed_outcome
#' @param smoothed_outcome
#' 
#' 
CalcEstimates <- 
  function(Wt, Yt = NULL, roads = NULL, Xt = NULL, Xt2 = NULL,
                          ps_coefs = NULL, ps_lags = NULL, intervention,
                          interv_cov = NULL, times = NULL,
                          truncation_level = NULL,
                          spat_kernel = 'gaussian', spat_adjust = 0.5,
                          spat_sigma = 1 / 8, edge_cor = FALSE, weights = NULL,
                          return_smoothed_outcome = FALSE,
                          smoothed_outcome = NULL) {
  
  
  # -------- PART 1: The setup ----------- #
  
  if (all(class(Wt) != 'hyperframe')) {
    Wt <- hyperframe(Points = Wt)
  }
  if (!is.null(Yt) & all(class(Yt) != 'hyperframe')) {
    Yt <- hyperframe(Points = Yt)
  }
  
  time_points <- length(Wt$Points)
  window <- Wt$Points[[1]]$window
  M <- dim(intervention)[3]

  if (is.null(times)) {
    start_time <- max(ps_lags + 1, M)
    times <- start_time : time_points
  }
  
  if (is.null(interv_cov)) {
    interv_cov <- roads$Priority[[1]]
  }
  
  # Items that will be returned:
  
  surfaces <- NULL
  surfaces_haj <- NULL
  
  
  
  # --------- PART 2: The weights ------------- #

  # If not provided, getting the weights we will use in the estimation.
  
  if (is.null(weights)) {
    log_weights <- GetLogWeights(roads = roads, Xt = Xt, Xt2 = Xt2, Wt = Wt,
                                 Yt = Yt, ps_coefs = ps_coefs,
                                 ps_lags = ps_lags,
                                 intervention = intervention,
                                 interv_cov = interv_cov, times = times)
    weights <- exp(apply(log_weights, c(1, 2), sum))
  }
  num_interv <- dim(weights)[1]
    
  if (!is.null(truncation_level)) {
    truncate_at <- apply(weights, 1, quantile, probs = truncation_level)
    for (ii in 1 : num_interv) {
      for (tt in 1 : ncol(weights)) {
        weights[ii, tt] <- min(weights[ii, tt], truncate_at[ii])
      }
    }
  }
  
  
  # --------- PART 3: The outcome surfaces ---------- #
  
  # smoothing the outcome at the time points we are considering.
  if (is.null(smoothed_outcome)) {
    for (tt in 1 : length(times)) {
      # smooooth <- density.ppp(Yt$Points[[times[tt]]], edge = edge_cor,
      #                         diggle = TRUE, kernel = spat_kernel,
      #                         adjust = spat_adjust, sigma = spat_sigma)
      smoothed_outcome[[tt]] <- pixellate.ppp(Yt$Points[[times[tt]]]) 
      # smoothed_outcome[[tt]] <- smooooth
    }
  }
  
  mat_im <- sapply(smoothed_outcome, function(x) as.matrix.im(x))
  pixels <- sqrt(dim(mat_im)[1])
  mat_im <- array(mat_im, dim = c(pixels, pixels, dim(mat_im)[2]))
  
  for (ii in 1 : dim(weights)[1]) {
    # Weighting each time point by the corresponding weights:
    mat_im_ii <- sweep(mat_im, MARGIN = 3, STATS = weights[ii, ], FUN = '*')
    surfaces[[ii]] <- mat_im_ii
    surfaces_haj[[ii]] <- surfaces[[ii]] / mean(weights[ii, ])

  }
  
  if (return_smoothed_outcome) {
    return(list(surfaces = surfaces, surfaces_haj = surfaces_haj,
                weights = weights, smoothed_outcome = smoothed_outcome))
  }
  
  return(list(surfaces = surfaces, surfaces_haj = surfaces_haj,
              weights = weights))
}