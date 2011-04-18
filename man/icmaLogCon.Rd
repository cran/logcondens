\name{icmaLogCon}
\alias{icmaLogCon}
\title{Computes a Log-Concave Probability Density Estimate via an Iterative Convex Minorant Algorithm}
\description{Given a vector of observations \eqn{{\bold{x}_n} = (x_1, \ldots, x_n)}{x_n = (x_1, \ldots, x_n)} 
with not necessarily equal entries, 
\code{\link{activeSetLogCon}} first computes vectors \eqn{{\bold{x}_m} = (x_1, \ldots, x_m)}{x_m = (x_1, \ldots, x_m)} and
\eqn{{\bold{w}} = (w_1, \ldots, w_m)}{w = (w_1, \ldots, w_m)} where \eqn{w_i} is the weight of each \eqn{x_i} s.t. 
\eqn{sum_{i=1}^m w_i = 1}.
Then, \code{\link{activeSetLogCon}} computes a concave, piecewise 
linear function \eqn{\widehat \phi_m} on \eqn{[x_1, x_m]} with knots only in \eqn{\{x_1, \ldots, x_m\}}{{x_1, \ldots, x_m}} such that 
  
  \deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{-\infty}^\infty \exp(\phi(t)) dt}{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \int_{-\infty}^\infty exp(\phi(t)) dt}
  
  is maximal. In order to be able to apply the pool - adjacent - violaters algorithm, computations are performed
  in the parametrization
  
  \deqn{{\bold{\eta}}({\bold{\phi}}) = \Bigl(\phi_1, \Bigl(\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m \Bigr).}{\eta(\phi) = (\phi_1, (\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i)_{i=2}^m ).}
  
  To find the maximum of \eqn{L}, a variant of the iterative convex minorant using the pool - adjacent - violaters 
  algorithm is used.}
\usage{icmaLogCon(x, xgrid = NULL, eps = 10^-8, T1 = 2000, 
    robustif = TRUE, print = FALSE)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, not necessarily equal.}
  \item{xgrid}{Governs the generation of weights for observations. See \code{\link{preProcess}} for details.}
  \item{eps}{An arbitrary real number, typically small. Iterations are halted if the directional derivative of \eqn{{\bold{\eta}} \to L({\bold{\eta}})} in the direction of the new candidate is \eqn{\le \varepsilon}.}
  \item{T1}{Maximal number of iterations to perform.}
  \item{robustif}{\code{robustif = TRUE} performs the robustification and Hermite interpolation procedure detailed in 
                Rufibach (2006, 2007), \code{robustif = FALSE} does not. In the latter case, convergence of the algorithm
                is no longer guaranteed.}
  \item{print}{\code{print = TRUE} outputs log-likelihood in every loop, \code{print = FALSE} does not. Make sure to 
  tell \code{R} to output (press CTRL+W).}
}
\value{
  \item{x}{Vector of observations \eqn{x_1, \ldots, x_m} that was used to estimate the density. }
  \item{w}{The vector of weights that had been used. Depends on the chosen setting for \code{xgrid}.}
  \item{f}{Vector with entries \eqn{\widehat f_m(x_i).}}
  \item{xn}{Vector with initial observations \eqn{x_1, \ldots, x_n}.}  
  \item{Loglik}{The value \eqn{L(\widehat \phi_m)} of the log-likelihood-function \eqn{L} at the maximum \eqn{\widehat \phi_m.}}
  \item{Iterations}{Number of iterations performed.}
  \item{sig}{The standard deviation of the initial sample \eqn{x_1, \ldots, x_n}.}
}
\references{
Rufibach K. (2006) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
\cr Available at \url{http://www.stub.unibe.ch/download/eldiss/06rufibach_k.pdf}.

Rufibach, K. (2007)
Computing maximum likelihood estimators of a log-concave density function.
\emph{J. Stat. Comput. Simul.} \bold{77}, 561--574.
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}
\seealso{
\code{\link{icmaLogCon}} can be used to estimate a log-concave density. However, to generate an object of
class \code{dlc} that allows application of \code{\link{summary}} and \code{\link{plot}} one has to 
use \code{\link{logConDens}}.

The following functions are used by \code{\link{icmaLogCon}}: 

\code{\link{phieta}}, \code{\link{etaphi}}, \code{\link{Lhat_eta}}, \code{\link{quadDeriv}},
\code{\link{robust}}, \code{\link{isoMean}}.
}
\examples{
set.seed(1977)
x <- rgamma(200, 2, 1)
\dontrun{
res <- icmaLogCon(x, T1 = 2000, robustif = TRUE, print = TRUE)

## plot resulting functions
par(mfrow = c(2, 1), mar = c(3, 2, 1, 2))
plot(x, exp(res$phi), type = 'l'); rug(x)
plot(x, res$phi, type = 'l'); rug(x)
}
}
\keyword{htest}
\keyword{nonparametric}

