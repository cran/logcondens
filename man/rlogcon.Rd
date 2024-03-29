\name{rlogcon}
\Rdversion{1.1}
\alias{rlogcon}
\title{Generate random sample from the log-concave and the smoothed log-concave density estimator}
\description{Generate a random sample from a distribution with density \eqn{\hat f_n} and \eqn{\hat f_n^*},
as described in Duembgen and Rufibach (2009, Section 3).}
\usage{rlogcon(n, x0)}
\arguments{
  \item{n}{Size of random sample to be generated.}
  \item{x0}{Sorted vector of independent and identically distributed numbers, not necessarily unique.}
}
\value{
\item{X}{Random sample from \eqn{\hat f_n}.}
\item{X_star}{Random sample from \eqn{\hat f_n^*}.}
\item{U}{Uniform random sample of size \code{n} used in the generation of \code{X}.}
\item{Z}{Normal random sample of size \code{n} used in the generation of \code{X_star}.}
\item{f}{Computed log-concave density estimator.}
\item{f.smoothed}{List containing smoothed log-concave density estimator, as output by \code{\link{evaluateLogConDens}}.}
\item{x}{Vector of distinct observations generated from \code{x0}.}
\item{w}{Weights corresponding to \code{x}.}
}
\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \doi{https://doi.org/10.18637/jss.v039.i06}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}
}
\examples{

## ===================================================
## Generate random samples as described in Section 3 of
## Duembgen and Rufibach (2009)
## ===================================================
x0 <- rnorm(111)
n <- 22
random <- rlogcon(n, x0)

## sample of size n from the log-concave density estimator
random$X

## sample of size n from the smoothed log-concave density estimator
random$X_star
}
\keyword{htest}
\keyword{nonparametric}
