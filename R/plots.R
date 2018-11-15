#' Plot of credibility intervals for each time interval, plotted as timeseries
#' 
#' @description The function plots the credibility intervals for each
#' time interval and returns a ggplot object, which is further customizable.
#' 
#' @param object model of class \code{\link{TemporalIso}}
#' @param prop double between 0 and 1: length of credibility interval. The default value is 80 percent.
#' @param plotShifts boolean if shifts shall be marked or not. Defaults to False.
#' @param ... arguments handed to \code{\link{getShiftTime}}
#' 
#' @return a \link[ggplot2]{ggplot} object.
#' 
#' @export
plotTime <- function(object, prop = 0.8, plotShifts = FALSE, ...){
  stopifnot(prop < 1)
  
  x <- getPlotData(object, prop = prop)
  breaks <- 1:length(object@time)
  p <- ggplot(x, aes_string(x = "time")) + geom_line(aes_string(y = "median")) +
    geom_line(aes_string(y = "lower"), size = 0.05) +
    geom_line(aes_string(y = "upper"), size = 0.05) +
    geom_ribbon(aes_string(ymin = "lower", ymax = "upper"), linetype = 2, alpha = 0.1) +
    geom_point(aes_string(x = "time", y = "median")) +
    scale_x_continuous(breaks = breaks,
                       labels = as.character(object@time)) +
    labs(title = paste0(prop * 100, "%-Credibility-Interval for isotopic values over time"),
         x = "Time", y = "Estimation") +
    theme(panel.grid.major.x = element_line(size = 0.1))
  
  if (plotShifts){
    index <- getShiftIndex(object, ...)
    p <- p + geom_vline(xintercept = breaks[which(index)] + 0.5, col = "darkgrey")
  }
  p
}

getPlotData <- function(object, prop = 0.8, time = NULL){
  lLim <- (1 - prop) / 2
  uLim <- 1 - lLim
  dat <- rstan::extract(object)$interval
  if (is.null(time)) time <- 1:ncol(dat)
  out <- as.data.frame(
    t(apply(dat, 2, quantile, probs = c(lLim, 0.5, uLim)))
  )
  names(out) <- c("lower", "median", "upper")
  out$time <- time
  return(out)
}
