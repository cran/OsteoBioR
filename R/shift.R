getShiftIndex <- function(object,
                          absolute = TRUE,
                          threshold = NULL,
                          probability = 0.5){

  if (is.null(threshold)){
    if (absolute){
      threshold <- 1.5
    }
    else{
      threshold <- 0.15
    }
    
  }
  
  samples <- as.data.frame(rstan::extract(object)$interval)
  names(samples) <- object@time
  diff <- t(apply(samples, 1, diff))
  if (!absolute){
    diff <- abs(diff / samples[, -NCOL(samples)])
  }
  q1 <- apply(diff, 2, quantile, probs = probability)
  q2 <- apply(diff, 2, quantile, probs = 1 - probability)
  out <- apply(rbind(q1, q2), 2, function(x) min(abs(x)))
  outsign <- sign(apply(rbind(q1, q2), 2, prod))
  
  out <- (out > threshold & outsign == 1)
  return(out)
}


#' Calculate time points of shifted isotopic values
#' 
#' @description The function calculates points in time where large changes happen in the isotopic values. It uses the posterior
#' distributions to estimate the probability of changes to be (absolutely or relatively) large.
#' 
#' @param object model of class \code{\link{TemporalIso}}
#' @param absolute boolean. If the calculation shall be based on absolute or relative differences. Defaults to TRUE.
#' @param threshold numeric. The threshold for a shift to be considered "large". Defaults to 1.5 for absolute isotopic values and
#' 15 percent for relative changes.
#' @param probability the probability for the differences to be larger than the threshold. Defaults to 50 percent.
#' 
#' @return a data.frame containing the interval starts (intStart) and ends (intEnd) of changes.
#' 
#' @export
getShiftTime <- function(object,
                         absolute = TRUE,
                         threshold = NULL,
                         probability = 0.5){
  index <- getShiftIndex(object = object,
                         absolute = absolute,
                         threshold = threshold,
                         probability = probability)
  
  out <- data.frame(
    intStart = object@time[c(index, FALSE)],
    intEnd = object@time[c(FALSE, index)]
  )
  return(out)
  
}
