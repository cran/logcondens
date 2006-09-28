\name{quadDeriv}
\alias{quadDeriv}
\title{Gradient and Diagonal of Hesse Matrix of Quadratic Approximation to Log-Likelihood Function L}
\description{
  Computes gradient and diagonal of the Hesse matrix w.r.t. to \eqn{\eta} of a quadratic approximation to the 
  reparametrized original log-likelihood function 
  
  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{-\infty}^{\infty} \exp(\varphi(t)) dt. }
  
  where \eqn{L} is parametrized via
  
  \deqn{\eta(\varphi) = \Bigl(\varphi_1, \Bigl(\eta_1+ \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr). }
  
  \eqn{\varphi}: vector \eqn{(\varphi(x_i))_{i=1}^m} representing concave, piecewise linear function \eqn{\varphi},\cr \eqn{\eta}: vector representing successive slopes of \eqn{\varphi.}
}
\usage{quadDeriv(dx, w, eta)}
\arguments{
  \item{dx}{Vector \eqn{(0, x_i-x_{i-1})_{i=2}^m.}}
  \item{w}{Vector of weights as in \code{\link{activeSetLogCon}}.}
  \item{eta}{Vector \eqn{\eta.}}
}
\value{\eqn{m \times 2} matrix. First column contains gradient and second column diagonal of Hesse matrix.}
\seealso{\code{\link{quadDeriv}} is used by the function \code{\link{icmaLogCon}}.}

\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\keyword{htest}
\keyword{nonparametric}
