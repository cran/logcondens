\name{isoMean}
\alias{isoMean}
\title{Pool-Adjacent Violaters Algorithm: Least Square Fit under Monotonicity Constraint}
\description{
  Fits a vector \eqn{\widehat g} with nondecreasing components to the data vector \eqn{y} such that 
  
  \deqn{\sum_{i=1}^n (y_i - \widehat g_i)^2 }
  
  is minimal (pool - adjacent - violators algorithm). In case a weight vector with positive entries (and the same size as \eqn{y}) is provided, the function produces an isotonic vector minimizing 
    
  \deqn{\sum_{i=1}^n w_i(y_i - \widehat g_i)^2 .}
}
\usage{isoMean(y, w)}
\arguments{
  \item{y}{Vector \eqn{(y_1, \ldots, y_n)} of data points.}
  \item{w}{Arbitrary vector \eqn{(w_1, \ldots, w_n)} of weights.}
}
\value{Returns vector \eqn{\widehat g}.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\examples{
## simple regression model
n <- 50
x <- sort(runif(n, 0, 1))
y <- x ^ 2 + rnorm(n, 0, 0.2)
plot(seq(0, 1, by = 0.01), seq(0, 1, by = 0.01) ^ 2, col = 2, type = 'l', xlim = range(c(0, 1, x)), ylim = range(c(0, 1 , y))); rug(x)

## plot pava result
lines(x, isoMean(y, rep(1 / n, n)), type = 's')
}
\keyword{ htest }
\keyword{ nonparametric }
