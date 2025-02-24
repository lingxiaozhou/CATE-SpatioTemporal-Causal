#' Creating the history of treatment and outcome point patterns
#'
#' Getting the treatment and outcome point patterns from previous time points
#' that affect current time points.
#' 
#' @param tt The time point whose history we want
#' @param Wt The treatment point patterns.
#' @param Yt The outcome point patterns.
#' @param center_variable Which variable's history we want. Options are 'W' for
#' treatment and 'Y' for outcome. If 'Y' is used, the time point tt is one of
#' the time points included in the history.
#' @param lags The lags of previous treatment and outcome point patterns that
#' matter.
#' @param first_only Logical, defaults to FALSE. If TRUE, only the history of
#' the first variable (W) will be returned.
#' @param Logical,defaults to TRUE. If true, the superimpose all points in the history
#' 
#' @return A list of two hyperframes, one for history of treatments and one for
#' history of outcomes.
#' 
GetHistory <- function(tt, Wt, Yt, center_variable = c('W', 'Y'), lags,
                       priority_coefs = c(2, 2), window = NULL,
                       first_only = FALSE, all= TRUE,dimyx = c(128,128)) {

  if (is.null(window)) {
    window <- owin(c(0, 1), c(0, 1))
  }
  
  # up_to_time is length 2 corresponding to previous treatment and outcome.
  up_to_time <- sapply(lags, function(x) max(1, tt - x))
  if (center_variable == 'Y') {
    up_to_time[1] <- max(1, tt - lags[1] + 1)
  }
  
  
  # ------- PART A: Creating the treatment history. --------- #
  
  # If the center variable is W, we get previous treatments over the time
  # period up_to_time[1] until tt - 1.
  # If the center is Y, the time period is up to tt.
  Wt_hist <- as.ppp(X = matrix(0, nrow = 0, ncol = 2), W = window)
  start_time <- tt - 1
  if (center_variable == 'Y') {
    start_time <- tt
  }
  
  # If the time period is one time point, then this is the time point we use,
  # otherwise, we superimpose all points.
  if (start_time == up_to_time[1]) {
    Wt_hist <- Wt[[start_time]]
  } else if (start_time > up_to_time[1]) {
    if(all == TRUE){
      Wt_hist <- do.call('superimpose', Wt[start_time : up_to_time[1]])
      Wt_hist$marks <- NULL
    }else{
      Wt_hist <- Wt[start_time : up_to_time[1]]
    }

  }
  
  Wt_hist <- GetPriority(hyperframe(Points = Wt_hist),
                         priority_coef = priority_coefs[1],dimyx = dimyx)
  
  if (first_only) {
    return(list(Wt_hist = Wt_hist))
  }
  
  # ------- PART B: Creating the outcome history. --------- #
  
  Yt_hist <- as.ppp(X = matrix(0, nrow = 0, ncol = 2), W = window)
  if (tt > 1) {
    Yt_hist <- Yt[(tt - 1)]
    if(all==TRUE){
      if (up_to_time[2] < tt - 1) {
        Yt_hist <- do.call('superimpose', Yt[(tt - 1) : up_to_time[2]])
        Yt_hist$marks <- NULL
      }
    }else{
      Yt_hist <- Yt[(tt - 1) : up_to_time[2]]
    }

  }
  
  Yt_hist <- GetPriority(hyperframe(Points = Yt_hist),
                         priority_coef = priority_coefs[2],dimyx = dimyx)

  return(list(Wt_hist = Wt_hist, Yt_hist = Yt_hist))
  
}