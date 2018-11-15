#' The 'OsteoBioR' package.
#' 
#' @description A DESCRIPTION OF THE PACKAGE
#' 
#' @docType package
#' @name OsteoBioR-package
#' @aliases OsteoBioR
#' @useDynLib OsteoBioR, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#' @importFrom ggplot2 ggplot geom_line geom_point geom_ribbon labs scale_x_continuous theme ggtitle
#' aes_string element_line scale_y_continuous geom_vline
#' @importFrom stats median quantile dnorm approx
#' @importFrom utils write.csv
#' 
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version NA. http://mc-stan.org
#' 
#' 
globalVariables(".")
NULL
