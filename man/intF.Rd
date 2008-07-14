\name{intF}
\alias{intF}
\title{Computes the Integral of the estimated CDF at Arbitrary Real Numbers in s}
\description{Based on the output of the function \code{\link{activeSetLogCon}}, this gives values of 

\deqn{\widehat I(t) = \int_{x_1}^t \widehat{F}(r) d r}

at all numbers in \eqn{\bold{s}}. Note that \eqn{t} (so all elements in \eqn{\bold{s}}) must lie in \eqn{[x_1,x_m]}.
The exact formula for \eqn{\widehat I(t)} is

\deqn{\widehat I(t)  = \Bigl(\sum_{i=1}^{i_0} \widehat{I}_i(x_{i+1})\Bigr)+\widehat{I}_{i_0}(t) }

where \eqn{i_0 = \ }min\eqn{\{ m-1 \,  , \ \{i \ : \ x_i \le t \} \} } and

\deqn{I_j(x) = \int_{x_j}^x \widehat{F}(r) d r = (x-x_j)\widehat{F}(x_j)+\Delta x_{j+1}\Bigl(\frac{\Delta x_{j+1}}{\Delta \widehat\varphi_{j+1}}J\Bigl(\widehat\varphi_j,\widehat\varphi_{j+1}, \frac{x-x_j}{\Delta x_{j+1}}\Bigr)-\frac{\widehat f(x_j)(x-x_j)}{\Delta \widehat \varphi_{j+1}}\Bigr)}

for \eqn{x \in [x_j, x_{j+1}], \ j = 1,\ldots,m-1}, \eqn{\Delta v_{i+1} = v_{i+1}-v_i} for any vector \eqn{\bold{v}}
and the function \eqn{J} introduced in \code{\link{Jfunctions}}.

}
\usage{intF(s, x, phi, Fhat)}
\arguments{
  \item{s}{Vector of real numbers where the functions should be evaluated at.}
  \item{x}{Vector \eqn{{\bold{x}} = (x_1, \ldots, x_m)} of original observations (sorted).}
  \item{phi}{Vector \eqn{(\widehat \varphi_m(x_i))_{i=1}^m}, as computed by \code{\link{activeSetLogCon}}.}
  \item{Fhat}{Vector \eqn{(\widehat F_{m,i})_{i=1}^m} with entries \deqn{\widehat F_{m,i} = \int_{x_1}^{x_i} \exp(\widehat \varphi_m(t)) dt,} as computed by \code{\link{activeSetLogCon}}.}
}
\value{Vector of the same length as \eqn{\bold{s}}, containing the values of \eqn{\widehat I} at the elements of \eqn{\bold{s}}.}

\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}
\seealso{
This function uses the output of \code{\link{activeSetLogCon}}. The function \code{\link{intECDF}} is similar, but based
on the empirical distribution function.
}

\references{
Duembgen, L. and Rufibach, K. (2008)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency.
Bernoulli, to appear.

Rufibach K. (2006) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
\cr Available at \url{http://www.stub.unibe.ch/download/eldiss/06rufibach_k.pdf}.
}

\examples{
## estimate gamma density
set.seed(1977)
x <- sort(rgamma(200, 2, 1))
res <- activeSetLogCon(x, w = NA, print = FALSE)

## compute and plot the process D(t) in Duembgen and Rufibach (2008)
s <- seq(min(x), max(x), by = 10 ^ -3)
D1 <- intF(s, x, res$phi, res$F)
D2 <- intECDF(s, x)
par(mfrow = c(2, 1))
plot(x, res$phi, type = 'l'); rug(x)
plot(s, D1 - D2, type = 'l'); abline(h = 0, lty = 2)
}
\keyword{htest}
\keyword{nonparametric}
