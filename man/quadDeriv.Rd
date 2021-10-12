\name{quadDeriv}
\alias{quadDeriv}
\title{Gradient and Diagonal of Hesse Matrix of Quadratic Approximation to Log-Likelihood Function L}
\description{
  Computes gradient and diagonal of the Hesse matrix w.r.t. to \eqn{\eta} of a quadratic approximation to the 
  reparametrized original log-likelihood function 
  
  \deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{-\infty}^{\infty} \exp(\phi(t)) dt. }{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - int_{-\infty}^{\infty} exp(\phi(t)) dt. }
  
  where \eqn{L} is parametrized via
  
  \deqn{{\bold{\eta}}({\bold{\phi}}) = \Bigl(\phi_1, \Bigl(\eta_1+ \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m\Bigr).}{\eta(\phi) = (\phi_1, (\eta_1 + \sum_{j=2}^i (x_i-x_{i-1}) \eta_i)_{i=2}^m).}
  
  \eqn{{\bold{\phi}}}{\phi}: vector \eqn{(\phi(x_i))_{i=1}^m} representing concave, piecewise linear function \eqn{\phi},\cr \eqn{{\bold{\eta}}}{\eta}: vector representing successive slopes of \eqn{\phi.}
}
\usage{quadDeriv(dx, w, eta)}
\arguments{
  \item{dx}{Vector \eqn{(0, x_i-x_{i-1})_{i=2}^m.}}
  \item{w}{Vector of weights as in \code{\link{activeSetLogCon}}.}
  \item{eta}{Vector \eqn{{\bold{\eta}}.}{\eta.}}
}
\value{\eqn{m \times 2} matrix. First column contains gradient and second column diagonal of Hesse matrix.}
\seealso{\code{\link{quadDeriv}} is used by the function \code{\link{icmaLogCon}}.}

\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}}

\keyword{htest}
\keyword{nonparametric}
