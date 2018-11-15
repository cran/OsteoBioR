#' Function to estimate temporal changes in isotopic values
#' 
#' @description Given renewal rates for different isotopic probes over time, a model
#' estimates a normal distribution of the isotopic values for each time interval. Out of
#' the given rates of renewal, first the influence of each time interval on the final
#' isotopic value is calculated. This proportion is used as a regressor in a fully Bayesian
#' linear model. The variance is estimated as rbf-Kernel-Matrix, where the hyperparameters are
#' chosen such that the correlation between time intervals close to each other is about 0.5.
#' The estimation is implemented using \code{link[rstan]{rstan}}.
#' 
#' @param renewalRates A dataframe specifying the renewal rates of different probes for
#' each time interval. The renewalRates should be between 0 and 100 (percentages). The dataframe should include a column specifying a time-index (e.g. 1, 2, 3, ...)
#' as well as columns for the different bones. The renewal rates should contain the times of origin,
#' containing 100.
#' @param timeVars A character string specifying the name of the column indicating the time.
#' @param boneVars A vector of character strings indicating the relevant variables containing the
#' renewal rates of bones and teeth. If not given, the variables with bones are all variables of the dataframe apart from the time variables.
#' @param isoMean A numeric number indicating the mean of the isotopic values measured.
#' @param isoSigma A numeric, positive number indicating the standard deviation of the isotopic values measured.
#' @param mc A boolean indicating if multiple cores should be used. If \code{TRUE}, which is the default, 4 cores are used.
#' @param adapt_delta A numeric value between 0 and 1 controlling the samplers behavior. Defaults to 0.9999. See \code{\link[rstan]{stan}} for more details.
#' @param max_treedepth A numeric, positive value controling the NUTS sampler. Defaults to 25. See \code{\link[rstan]{stan}} for more details.
#' @param mu_df Hyperparameter for the mean of the interval estimates: degrees of freedom of a non-standardized t-Student distribution. Defaults to 1.
#' @param mu_mean Hyperparameter for the mean of the interval estimates: mean of a non-standardized t-Student distribution. Defaults to 0.
#' @param mu_sd Hyperparameter for the mean of the interval estimates: standard deviation of a non-standardized t-Student distribution. Defaults to 10.
#' @param rho_mean Hyperparameter for the length scale of the rbf-kernel: mean of a normal. Defaults to 1.
#' @param rho_sd Hyperparameter for the length scale of the rbf-kernel: standard deviation of a normal. Defaults to 0.25.
#' @param alpha_mean Hyperparameter for the signal variance of the rbf-kernel: mean of a normal. Defaults to 2.
#' @param alpha_sd Hyperparameter for the signal variance of the rbf-kernel: standard deviation of a normal. Defaults to 0.5.
#' @param chains Number of chains for mcmc. Defaults to 4
#' @param iter Number of iterations per chain for mcmc. Defaults to 2000
#' @param ... Arguments passed to rstand \link[rstan]{sampling}
#'
#' @return A fitted object of class \link{TemporalIso}.
#' 
#' @examples
#' \dontrun{
#' data <- data.frame(
#' intStart = 0:5,
#' intEnd = 1:6,
#' bone1 = c(100, 50, 20, 10, 5, 2),
#' bone2 = c(100, 10, 5, 1, 1, 1),
#' tooth1 = c(0, 100, 0, 0, 0, 0),
#' tooth2 = c(0, 0, 100, 0, 0, 0)
#' )
#' y_mean <- c(-10, -7, -12, -9)
#' y_sigma <- c(2, 1.5, 2.5, 2.5)
#' fit <- estimateIntervals(renewalRates = data,
#'                          timeVars = "intStart",
#'                          boneVars = c("bone1", "bone2", "tooth1", "tooth2"),
#'                          isoMean = y_mean,
#'                          isoSigma = y_sigma)
#' print(fit)
#' plot(fit)
#' plotTime(fit)
#' 
#' # get estimates for specific time points
#' estimateTimePoint(fit, time = seq(0,5, by = 0.5))
#' 
#' # shift point estimation
#' plotTime(fit, plotShifts = TRUE, threshold = 0.5)
#' getShiftTime(fit, threshold = 0.5)
#' 
#' #Staying time estimation
#' estimatedStayTimes <- getSiteStayTimes(object = fit,
#' siteMeans = c(-8, -10),
#' siteSigma = c(1, 1.5))
#' }
#' 
#' @seealso \link[rstan]{sampling}
#'
#' @export
#' 
#' 
estimateIntervals <- function(renewalRates,
                              timeVars,
                              boneVars = NULL,
                              isoMean,
                              isoSigma,
                              mu_df = 1,
                              mu_mean = 0,
                              mu_sd = 10,
                              rho_mean = 1,
                              rho_sd = 0.25,
                              alpha_mean = 2,
                              alpha_sd = 0.5,
                              mc = TRUE,
                              adapt_delta = 0.9999,
                              max_treedepth = 25,
                              chains = 4,
                              iter = 2000,
                              ...) {
  
  stopifnot(
    isoSigma > 0,
    max(renewalRates[boneVars]) <= 100
  )

  if (is.null(boneVars))
  {
    boneVars <- names(renewalRates)[!names(renewalRates) %in% timeVars]
  }
  if (length(timeVars) > 1) timeVars <- timeVars[1]
  
  N <- ncol(renewalRates[boneVars])
  T <- nrow(renewalRates[timeVars])
  y_mean <- isoMean
  y_sigma <- isoSigma
  time <- renewalRates[[timeVars]]
  t <- time / min(abs(diff(time))) # normalization: nearest neighbour is 1 away
  x <- t(as.matrix(apply(renewalRates[boneVars], 2, calcInfluence)))  
  
  cores <- getOption("mc.cores", if (mc) min(4, chains) else 1)
  model <- suppressWarnings(sampling(stanmodels$linRegGP,
                                     chains = chains,
                                     iter = iter,
                                     cores = cores,
                                     verbose = FALSE,
                                     refresh = 0,
                                     control = list(adapt_delta = adapt_delta,
                                                    max_treedepth = max_treedepth),
                                     ...))
  
  return(new("TemporalIso",
             model,
             regressor = x,
             time = time))
}

#' Function to estimate isotopic value for specific time point(s)
#' 
#' @description Once the isotopic values in time course are estimated along
#' with estimateIntervals(), this function provides an interface to extract the mean, standard deviation and quantiles for a specific point in time
#' @param object An object of class \link{TemporalIso} fitted with estimateIntervals()
#' @param time A vector of numerics, indicating the point in time to extract estimates
#' @param intervalProb A numeric between 0 and 1 indicating the coverage of the credible interval
#' 
#' @return A data.frame with time points and their means, standard deviations and credible intervals
#' 
#' @export
estimateTimePoint <- function(object, time, intervalProb = 0.95){
  if (any(time < min(object@time) | time > max(object@time))){
    stop("Provided time must be inside fitted intervals")
  }
  if (intervalProb <= 0 | intervalProb >= 1){
    stop("Value of intervalProb must be between 0 and 1")
  }
  
  estimates <- rstan::extract(object)$interval
  mean <- approx(object@time, colMeans(estimates), xout = time)$y
  sd <- approx(object@time, apply(estimates, 2, sd), xout = time)$y
  intLower <- approx(object@time, apply(estimates, 2,
                                        quantile, (1 - intervalProb) / 2), xout = time)$y
  intUpper <- approx(object@time, apply(estimates, 2,
                                        quantile, 1 - ( (1 - intervalProb) / 2)), xout = time)$y
  return(data.frame(time = time, mean = mean, sd = sd, intLower = intLower, intUpper = intUpper))
}

calcInfluence <- function(x){
  if (max(x) > 1) x <- x / 100
  out <- rep(0, length(x))
  temp <- 1
  for (i in length(x):1){
    out[i] <- x[i] * temp
    temp <- temp - out[i]
  }
  return(out)
}
