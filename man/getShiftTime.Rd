% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift.R
\name{getShiftTime}
\alias{getShiftTime}
\title{Calculate time points of shifted isotopic values}
\usage{
getShiftTime(object, absolute = TRUE, threshold = NULL,
  probability = 0.5)
}
\arguments{
\item{object}{model of class \code{\link{TemporalIso}}}

\item{absolute}{boolean. If the calculation shall be based on absolute or relative differences. Defaults to TRUE.}

\item{threshold}{numeric. The threshold for a shift to be considered "large". Defaults to 1.5 for absolute isotopic values and
15 percent for relative changes.}

\item{probability}{the probability for the differences to be larger than the threshold. Defaults to 50 percent.}
}
\value{
a data.frame containing the interval starts (intStart) and ends (intEnd) of changes.
}
\description{
The function calculates points in time where large changes happen in the isotopic values. It uses the posterior
distributions to estimate the probability of changes to be (absolutely or relatively) large.
}
