\name{preProcess}
\alias{preProcess}
\title{Compute a weighted sample from initial observations}
\description{Generates weights from initial sample.}
\usage{preProcess(x, xgrid = NULL)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, not necessarily unique.}
  \item{xgrid}{Parameter that governs the generation of weights: If \code{xgrid = NULL} a new sample
  of unique observations is generated with corresponding vector of weights. If \code{xgrid} is
  a positive number, observations are binned in a grid with grid length \code{xgrid}.
  Finally, an entire vector specifying a user-defined grid can be supplied.}
}
\value{
\item{x}{Vector of unique and sorted observations deduced from the input \code{x} according to the specification 
given by \code{xgrid}.}
\item{w}{Vector of corresponding weights, normalized to sum to one.}
\item{sig}{Standard deviation of the inputed observations. This quantity is needed when computing the smoothed
log-concave density estimator via \code{\link{evaluateLogConDens}}.}
\item{n}{Number of initial observations.}
}
\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
