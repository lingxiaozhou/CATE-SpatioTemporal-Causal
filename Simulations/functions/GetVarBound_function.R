GetVarBound <- function(smoothed_outcome = NULL, out_in_B = NULL, weights,
                        B = NULL) {
  
  if (any((is.null(smoothed_outcome) | is.null(B))) & is.null(out_in_B)) {
    stop('Specify (smoothed_outcome and B) or (out_in_B).')
  }
  
  num_interv <- dim(weights)[1]
  num_B <- length(B)

  if (is.null(out_in_B)) {  # Integrate the smoothed outcome over sets B.
    
    num_time_points <- length(smoothed_outcome)
    if (num_time_points != ncol(weights)) {
      stop('Dimensions of smoothed_outcome and weights do not agree.')
    }
    
    out_in_B <- matrix(NA, nrow = num_time_points, ncol = num_B)
    for (bb in 1 : num_B) {
      out_in_B[, bb] <- sapply(smoothed_outcome, function(s) {
        integral(s, domain = B[[bb]])})
    }
  }
  
  num_time_points <- nrow(out_in_B)
  num_B <- ncol(out_in_B)
  
  all_est <- array(NA, dim = c(num_interv, num_time_points, num_B))
  for (ii in 1 : num_interv) {
    for (bb in 1 : num_B) {
      all_est[ii, , bb] <- out_in_B[, bb] * weights[ii, ]
    }
  }
  
  bound <- apply(all_est, c(1, 3), function(x) mean(x ^ 2) / num_time_points)
  colnames(bound) <- names(B)
  rownames(bound) <- paste0('interv', 1 : num_interv)
  
  bound_t <- array(NA, dim = c(num_interv, num_interv, num_B))
  for (ii in 1 : (num_interv - 1)) {
    for (jj in (ii + 1) : num_interv) {
      est <- all_est[ii, , , drop = FALSE] - all_est[jj, , , drop = FALSE]
      bound_t[ii, jj, ] <- apply(est, 3, function(x) mean(x ^ 2) / num_time_points)
    }
  }
  
  mean_weight <- apply(weights, 1, mean)
  bound_haj <- sweep(bound, 1, mean_weight ^ 2, FUN = '/')
  
  all_est <- sweep(all_est, MARGIN = 1, mean_weight, FUN = '/')
  bound_t_haj <- array(NA, dim = c(num_interv, num_interv, num_B))
  for (ii in 1 : (num_interv - 1)) {
    for (jj in (ii + 1) : num_interv) {
      est <- all_est[ii, , , drop = FALSE] - all_est[jj, , , drop = FALSE]
      bound_t_haj[ii, jj, ] <- apply(est, 3, function(x) mean(x ^ 2) / num_time_points)
    }
  }
  
  return(list(bound_po = bound, bound_tau = bound_t,
              bound_po_haj = bound_haj, bound_tau_haj = bound_t_haj))
  
}
