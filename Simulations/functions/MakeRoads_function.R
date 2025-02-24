MakeRoads <- function(window, second_road = FALSE)  {

  window_length <- 1
  if (is.null(window)) {
    window <- owin(c(0, 1), c(0, 1))
  }
  
  x <- seq(0, window_length / 2, length.out = 3000)
  roads <- cbind(x, sqrt((window_length / 2) ^ 2 - x ^ 2))
  
  x <- seq(0, window_length, length.out = 3000)
  ysq <- (window_length / 2) ^ 2 - (x - 0.7) ^ 2
  roads <- rbind(roads, cbind(x[ysq > 0], 1 - sqrt(ysq[ysq > 0])))
  
  x <- seq(0, window_length / 2, length.out = 3000)
  roads <- rbind(roads, cbind(x, 0.5 + 0.7 * x))
  
  x <- seq(0, 1, length.out = 3000)
  roads <- rbind(roads, cbind(x, 2.5 * x - 1.2))
  
  roads <- na.omit(roads)
  roads <- roads[apply(roads, 1, max) < 1 & apply(roads, 1, min) > 0, ]
  
  roads <- spatstat::as.ppp(roads, window)

  roads <- hyperframe(Points = roads)
  roads <- GetPriority(point_patterns = roads)
  roads$Priority <- roads$Priority / with(roads, integral(Priority, domain = window))
  
  
  if (second_road) {
    x <- seq(0, window_length / 3, length.out = 1000) * c(- 1, 1) + window_length / 2
    y <- sqrt(1 - x ^ 2)
    x2 <- c(x, x, y)
    y2 <- c(y, 1 - y, x)
    roads2 <- hyperframe(Points = as.ppp(cbind(x = x2, y = y2), window))
    roads2 <- GetPriority(point_patterns = roads2, priority_coef = 3)
    roads <- rbind(roads, roads2)
  }
  
  
  return(roads)
  
}