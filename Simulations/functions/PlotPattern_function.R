PlotPattern <- function(pattern, exclude = NULL) {
  
  if (is.null(exclude)) {
    x <- do.call('superimpose', pattern$Points)
  } else {
    x <- do.call('superimpose', pattern$Points[- exclude])
  }
  x$marks <- NULL
  
  
  num_points <- with(pattern, Points$n)
  if (!is.null(exclude)) {
    num_points <- num_points[- exclude]
  }
  
  par(mfrow = c(1, 2), mar = c(3, 2, 3, 0))
  hist(num_points, main = paste0('Number - mean ', round(mean(num_points), 1)))
  plot(x, pch = 16, cex = 0.2, main = 'Location')
  
  
}