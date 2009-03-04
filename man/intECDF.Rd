\name{intECDF}
\alias{intECDF}
\title{Computes the Integrated Empirical Distribution Function at Arbitrary Real Numbers in s}
\description{Computes the value of

\deqn{\bar{I}(t) = \int_{x_1}^t \bar{F}(r) d \, r}

where \eqn{\bar F} is the empirical distribution function of \eqn{x_1,\ldots,x_m}, at all real numbers \eqn{t} in the 
vector \eqn{\bold{s}}. Note that \eqn{t} (so all elements in \eqn{\bold{s}}) must lie in \eqn{[x_1,x_m]}. The exact formula for 
\eqn{\bar I(t)} is

\deqn{\bar I(t)  = \Big(\sum_{i=2}^{i_0}(x_i-x_{i-1})\frac{i-1}{n} \Big) + (t-x_{i_0})\frac{i_0-1}{n}}

where \eqn{i_0 = \max_{i=1,\ldots,m} \{x_i \le t\}}.
}
\usage{intECDF(s, x)}
\arguments{
  \item{s}{Vector of real numbers in \eqn{[x_1,x_m]} where \eqn{\bar{I}} should be evaluated at.}
  \item{x}{Vector \eqn{{\bold{x}} = (x_1, \ldots, x_m)} of original observations (sorted).}
}
\value{Vector of the same length as \eqn{\bold{s}}, containing the values of \eqn{\bar I} at the elements of \eqn{\bold{s}}.}

\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}
\seealso{
This function together with \code{\link{intF}} can be used to check the characterization of the log-concave density
estimator in terms of distribution functions, see Rufibach (2006) and Duembgen and Rufibach (2009).
}

\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Rufibach K. (2006) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
\cr Available at \url{http://www.stub.unibe.ch/download/eldiss/06rufibach_k.pdf}.
}

\examples{
# for an example see the function intF.
}
\keyword{htest}
\keyword{nonparametric}
