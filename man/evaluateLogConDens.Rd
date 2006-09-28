\name{evaluateLogConDens}
\alias{evaluateLogConDens}
\title{Evaluates the Log-Density Estimator at an Arbitrary Real Number x0}
\description{Based on the output of the function \code{\link{activeSetLogCon}}, this gives the values of the 
functions 

\deqn{\widehat \varphi_m(t)}
\deqn{\widehat f_m(t) = \exp(\widehat \varphi_m(t))}
\deqn{\widehat F_m(t) = \int_{x_1}^t \exp(\widehat \varphi_m(t)) dt}

at an arbitrary real number \eqn{t = x_0}.}
\usage{evaluateLogConDens(x0, x, f, F, IsKnot)}
\arguments{
  \item{x0}{Real number where the functions should be evaluated at.}
  \item{x}{Vector \eqn{x = (x_1, \ldots, x_m)} of original observations (sorted).}
  \item{f}{Vector \eqn{(\widehat f_m(x_i))_{i=1}^m}, as computed by \code{\link{activeSetLogCon}}.}
  \item{F}{Vector \eqn{(\widehat F_m(x_i))_{i=1}^m}, as computed by \code{\link{activeSetLogCon}}.}
  \item{IsKnot}{0-1 vector giving the kinks of \eqn{\widehat \varphi_m}, as computed by \code{\link{activeSetLogCon}}.}
}
\value{Vector \eqn{(\widehat \varphi_m(x_0), \widehat f_m(x_0), \widehat F_m(x_0))}.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}
\seealso{
This function uses the output of \code{\link{activeSetLogCon}}.

For log-concave density estimation via an iterative convex minorant algorithm see \code{\link{icmaLogCon}}.
}
\examples{
## estimate gamma density
set.seed(1977)
x <- sort(rgamma(200, 2, 1))
res <- activeSetLogCon(x, w = NA, print = FALSE)

## compute function values at an arbitrary point
x0 <- (x[100] + x[101]) / 2
evaluateLogConDens(x0, x, exp(res$phi), res$F, res$IsKnot)
}
\keyword{htest}
\keyword{nonparametric}
