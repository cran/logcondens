\name{evaluateLogConDens}
\alias{evaluateLogConDens}
\title{Evaluates the Log-Density MLE and Smoothed Estimator at Arbitrary Real Numbers xs}
\description{Based on the output of \code{\link{activeSetLogCon}}, this function computes the values of  

\deqn{\widehat \varphi_m(t)}
\deqn{\widehat f_m(t) = \exp(\widehat \varphi_m(t))}
\deqn{\widehat F_m(t) = \int_{x_1}^t \exp(\widehat \varphi_m(x)) dx}
\deqn{\widehat f_m^*(t) = \exp(\widehat \varphi_m^*(t))}
\deqn{\widehat F_m^*(t) = \int_{x_1}^t \exp(\widehat \varphi_m^*(x)) dx}

at all real number \eqn{t} in \code{xs}. The exact formula for \eqn{\widehat F_m} and \eqn{t \in [x_j,x_{j+1}]} is

\deqn{\widehat F_m(t)  = \widehat F_m(x_j)  + (x_{j+1}-x_j) J\Big(\widehat \varphi_j, \widehat \varphi_{j+1}, \frac{t-x_j}{x_{j+1}-x_j} \Big) }

for the function \eqn{J} introduced in \code{\link{Jfunctions}}. Closed formulas can also be given for \eqn{\widehat f_m^*(t)}
and \eqn{\widehat F_m^*(t)}.
}
\usage{evaluateLogConDens(xs, res, which = 1:5, gam = NULL, print = FALSE)}
\arguments{
  \item{xs}{Vector of real numbers where the functions should be evaluated at.}
  \item{res}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens}.}
  \item{which}{A (sub-)vector of \code{1:5} specifying which of the above quantities should be computed.}
  \item{gam}{Only necessary if \code{smoothed = TRUE}. The standard deviation of the normal kernel. If equal to 
            \code{NULL}, \code{gam} is chosen such that the variances of the original sample \eqn{x_1, \ldots, x_n} 
            and \eqn{\widehat f_n^*} coincide. See \code{\link{logConDens}} for details.}
  \item{print}{Progress in computation of smooth estimates is shown.}
}
\value{Matrix with rows \eqn{(x_{0, i}, \widehat \varphi_m(x_{0, i}), \widehat f_m(x_{0, i}), \widehat F_m(x_{0, i}), \widehat f_m^*(x_{0, i}), \widehat F_m^*(x_{0, i}))}
where \eqn{x_{0,i}} is the \eqn{i}-th entry of \code{xs}.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}
\examples{
## estimate gamma density
set.seed(1977)
x <- rgamma(200, 2, 1)
res <- logConDens(x, smoothed = TRUE, print = FALSE)

## compute function values at an arbitrary point
xs <- (res$x[100] + res$x[101]) / 2
evaluateLogConDens(xs, res)

## only compute function values for non-smooth estimates
evaluateLogConDens(xs, res, which = 1:3)
}
\keyword{htest}
\keyword{nonparametric}
