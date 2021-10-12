\name{isoMean}
\alias{isoMean}
\title{Pool-Adjacent Violaters Algorithm: Least Square Fit under Monotonicity Constraint}
\description{
  Fits a vector \eqn{\widehat {\bold{g}}}{\hat g} with nondecreasing components to the data vector 
  \eqn{{\bold{y}}}{y} such that 
  
  \deqn{\sum_{i=1}^n (y_i - \widehat g_i)^2 }{\sum_{i=1}^n (y_i - \hat g_i)^2 }
  
  is minimal (pool - adjacent - violators algorithm). In case a weight vector with positive entries (and the same size as \eqn{{\bold{y}}}{y}) is provided, the function produces an isotonic vector minimizing 
    
  \deqn{\sum_{i=1}^n w_i(y_i - \widehat g_i)^2 .}{\sum_{i=1}^n w_i(y_i - \hat g_i)^2.}
}
\usage{isoMean(y, w)}
\arguments{
  \item{y}{Vector \eqn{(y_1, \ldots, y_n)} of data points.}
  \item{w}{Arbitrary vector \eqn{(w_1, \ldots, w_n)} of weights.}
}
\value{Returns vector \eqn{\widehat {\bold{g}}}{\widehat g}.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}}

\examples{
## simple regression model
n <- 50
x <- sort(runif(n, 0, 1))
y <- x ^ 2 + rnorm(n, 0, 0.2)
s <- seq(0, 1, by = 0.01)
plot(s, s ^ 2, col = 2, type = 'l', xlim = range(c(0, 1, x)), 
    ylim = range(c(0, 1 , y))); rug(x)

## plot pava result
lines(x, isoMean(y, rep(1 / n, n)), type = 's')
}
\keyword{ htest }
\keyword{ nonparametric }
