\name{intF}
\alias{intF}
\title{Computes the Integral of the estimated CDF at Arbitrary Real Numbers in s}
\description{Based on an object of class \code{dlc} as output by the function \code{\link{logConDens}}, 
this function gives values of 

\deqn{\widehat I(t) = \int_{x_1}^t \widehat{F}(r) d r}{\hat I(t) = \int_{x_1}^t \hat{F}(r) d r}

at all numbers in \eqn{\bold{s}}{s}. Note that \eqn{t} (so all elements in \eqn{\bold{s}}{s}) must lie in 
\eqn{[x_1,x_m]}. The exact formula for \eqn{\widehat I(t)} is

\deqn{\widehat I(t)  = \Bigl(\sum_{i=1}^{i_0} \widehat{I}_i(x_{i+1})\Bigr)+\widehat{I}_{i_0}(t)}{\hat I(t) = (\sum_{i=1}^{i_0} \hat{I}_i(x_{i+1}))+\hat{I}_{i_0}(t)}

where \eqn{i_0 = }min\eqn{\{m-1 \, , \ \{i \ : \ x_i \le t \}\}}{{m-1, {i : x_i \le t}}} and

\deqn{I_j(x) = \int_{x_j}^x \widehat{F}(r) d r = (x-x_j)\widehat{F}(x_j)+\Delta x_{j+1}\Bigl(\frac{\Delta x_{j+1}}{\Delta \widehat\phi_{j+1}}J\Bigl(\widehat\phi_j,\widehat\phi_{j+1}, \frac{x-x_j}{\Delta x_{j+1}}\Bigr)-\frac{\widehat f(x_j)(x-x_j)}{\Delta \widehat \phi_{j+1}}\Bigr)}{I_j(x) = int_{x_j}^x \hat{F}(r) d r = (x-x_j)\hat{F}(x_j)+\Delta x_{j+1}((\Delta x_{j+1})/(\Delta \hat\phi_{j+1})J(\hat\phi_j, \hat \phi_{j+1}, (x-x_j)/(\Delta x_{j+1}))-(\hat f(x_j)(x-x_j))/(\Delta \hat \phi_{j+1}))}

for \eqn{x \in [x_j, x_{j+1}], \ j = 1,\ldots, m-1}{x \in [x_j, x_{j+1}], j = 1, \ldots, m-1}, \eqn{\Delta v_{i+1} = v_{i+1} - v_i} for any vector \eqn{\bold{v}}{v}
and the function \eqn{J} introduced in \code{\link{Jfunctions}}.

}
\usage{intF(s, res)}
\arguments{
  \item{s}{Vector of real numbers where the functions should be evaluated at.}
  \item{res}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens}.}
}
\value{Vector of the same length as \eqn{\bold{s}}, containing the values of \eqn{\widehat I} at the elements of 
\eqn{\bold{s}}{s}.}

\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}}
\seealso{
This function uses the output of \code{\link{activeSetLogCon}}. The function \code{\link{intECDF}} is similar, but based
on the empirical distribution function.
}

\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \doi{https://doi.org/10.18637/jss.v039.i06}

Rufibach K. (2006) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
\cr Available at \url{https://slsp-ube.primo.exlibrisgroup.com/permalink/41SLSP_UBE/17e6d97/alma99116730175505511}.
}

\examples{
## estimate gamma density
set.seed(1977)
x <- rgamma(200, 2, 1)
res <- logConDens(x, smoothed = FALSE, print = FALSE)

## compute and plot the process D(t) in Duembgen and Rufibach (2009)
s <- seq(min(res$x), max(res$x), by = 10 ^ -3)
D1 <- intF(s, res)
D2 <- intECDF(s, res$xn)
par(mfrow = c(2, 1))
plot(res$x, res$phi, type = 'l'); rug(res$x)
plot(s, D1 - D2, type = 'l'); abline(h = 0, lty = 2)
}
\keyword{htest}
\keyword{nonparametric}
