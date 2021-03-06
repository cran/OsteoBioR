#' S4 class for temporal estimation of isotopic values
#' @description Inherits from \code{\link[rstan]{stanfit}}
#' @slot regressor regressor used for estimation (class matrix)
#' @slot time timestamps (class vector)
#' @export
TemporalIso <- setClass("TemporalIso",
                              slots = list(regressor = "matrix",
                                           time = "vector"),
                              contains = "stanfit")

#' Plots for TemporalIso objects
#' 
#' @description Inherited from \link[rstan]{rstan-plotting-functions}, the default plot shows
#' posterior uncertainty intervals and point estimates of isotopic values for each interval. The
#' plot method can also regress to standard rstan methods.
#'
#' @param x TemporalIso: An instance of class TemporalIso,
#' generated from \code{\link{estimateIntervals}}
#' 
#' @param plotfun character: a character string specifying which plotting function is to be used. The
#' default is "isotopic", whis is a \code{\link[rstan]{stan_plot}} only for the interval values.
#' "stan_plot", "stan_trace", "stan_dens" etc regress to the functionality of rstan.
#' 
#' @param ... Optional arguments to the plotting functions.
#' 
#' @return a ggplot object which can be further customized using the ggplot2 package.
#' 
#' @seealso \link[rstan]{rstan-plotting-functions}, \link[rstan]{rstan_gg_options}
#' 
#' @export
setMethod("plot",
          signature(x = "TemporalIso", y = "missing"),
          function(x, ..., plotfun = "isotopic") {
            if (plotfun == "isotopic"){
              suppressMessages(
                callNextMethod(x, pars = c("interval"), plotfun = "stan_plot", ...) +
                ggtitle("Isotopic Values", subtitle = "Credibility Intervals") +
                scale_y_continuous(breaks = 1:length(x@time), labels = rev(as.character(x@time)))
              )
            }
            else {
              callNextMethod()
            }
          })
