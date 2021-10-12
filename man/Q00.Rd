\name{Q00}
\alias{Q00}
\title{Numerical Routine Q}
\description{This function is used in the computation of \eqn{\widehat f^*}{\hat f^*} and \eqn{\widehat F^*}{\hat F^*}.}
\usage{
Q00(x, a, u, v, gamma, QFhat = FALSE)}
\arguments{
  \item{x}{Number at which to compute \eqn{q} and/or \eqn{Q}.}
  \item{a}{Vector of length \eqn{m} with real entries.}
  \item{u}{Vector of length \eqn{m} with real entries.}
  \item{v}{Vector of length \eqn{m} with real entries.}
  \item{gamma}{Real number. Standard deviation to be used.}
  \item{QFhat}{Logical. Should \eqn{Q} be computed?}
}
\value{The vector(s) \eqn{q} and/or \eqn{Q}.}
\note{Taylor approximation is used if \eqn{a} is small. In addition, as described in Duembgen and Rufibach (2011) at
the end of Appendix C, in extreme situations, e.g. when data sets contain extreme spacings, numerical problems may 
occur in the computation of the function \eqn{q_\gamma} (eq. (7) in Duembgen and Rufibach, 2011). 
For it may happen that the exponent is rather large while the difference of Gaussian CDFs is very 
small. To moderate these problems, we are using the following bounds:

\deqn{\exp(- m^2/2) \bigl( \Phi(\delta) - \Phi(-\delta) \bigr) \ \le \ \Phi(b) - \Phi(a) \ \le \ \exp(- m^2/2) \cosh(m\delta) \bigl( \Phi(\delta) - \Phi(-\delta) \bigr)}{exp(- m ^ 2 / 2) (\Phi(\delta) - \Phi(-\delta)) \le \Phi(b) - \Phi(a) \le exp(- m^2/2) cosh(m\delta) ( \Phi(\delta) - \Phi(-\delta)}

for arbitrary numbers \eqn{a < b} and \eqn{m := (a + b) / 2}, \eqn{\delta := (b - a) / 2}.

However, the function \code{Q00} is not intended to be invoked by the end user.
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}
}
\references{
Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \url{https://www.jstatsoft.org/v39/i06}
}


\keyword{htest}
\keyword{nonparametric}
