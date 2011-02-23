\name{ROCx}
\alias{ROCx}
\title{Compute ROC curve at a given x based on log-concave estimates for the constituent distributions}
\description{Computes the value of the ROC curve at \eqn{x} (which may be a vector) based on log-concave density estimates of the constituent distributions.}
\usage{ROCx(x, res0, res1, smooth = FALSE)}
\arguments{
    \item{x}{Vector of numbers in \eqn{[0, 1]} where the ROC curve should be computed at.}
    \item{res0}{\code{dlc} object as a result of a call to \code{\link{logConDens}} for the data of the controls.}
    \item{res1}{\code{dlc} object as a result of a call to \code{\link{logConDens}} for the data of the cases.}
    \item{smooth}{Logical. If \code{TRUE} kernel smoothed log-concave estimate is used.}
}
\value{A real number or a vector of dimension the same as \eqn{x}, the value of the ROC curve at \code{x}.}
\references{
Rufibach, K. (2011).
A smooth ROC curve estimator based on log-concave density estimates.
Preprint.
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html}}
\seealso{Used for the computation of AUC in \code{\link{logConROC}}.}
\keyword{nonparametric}
