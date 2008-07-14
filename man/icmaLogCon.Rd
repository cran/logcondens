\name{icmaLogCon}
\alias{icmaLogCon}
\title{Computes a Log-Concave Probability Density Estimate via an Iterative Convex Minorant Algorithm}
\description{Given a vector of observations \eqn{{\bold{x}} = (x_1, \ldots, x_m)} with pairwise distinct entries and
a vector of weights \eqn{{\bold{w}}=(w_1, \ldots, w_m)} s.t. \eqn{\sum_{i=1}^m w_i = 1}, \code{\link{icmaLogCon}}
computes iteratively a concave, piecewise linear function \eqn{\widehat \varphi_m} on \eqn{[x_1, x_m]} with knots only in
\eqn{\{x_1, \ldots, x_m\}} such that 
  
  \deqn{L(\varphi) = \sum_{i=1}^m w_i \varphi(x_i) - \int_{-\infty}^\infty \exp(\varphi(t)) dt}
  
  is maximal. In order to be able to apply the pool - adjacent - violaters algorithm, computations are performed
  in the parametrization
  
  \deqn{{\bold{\eta}}({\bold{\varphi}}) = \Bigl(\varphi_1, \Bigl(\eta_1 + \sum_{j=2}^i (x_i-x_{i-1})\eta_i\Bigr)_{i=2}^m \Bigr).} 
  
  To find the maximum of \eqn{L}, a variant of the iterative convex minorant using the pool - adjacent - violaters algorithm is used.}
\usage{icmaLogCon(x, w = NA, eps = 10^-8, T1 = 2000, robustif = TRUE, print = FALSE)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}}}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1 / m}.}
  \item{eps}{An arbitrary real number, typically small. Iterations are halted if the directional derivative of \eqn{{\bold{\eta}} \to L({\bold{\eta}})} in the direction of the new candidate is \eqn{\le \varepsilon}.}
  \item{T1}{Maximal number of iterations to perform.}
  \item{robustif}{robustif = TRUE performs the robustification and Hermite interpolation procedure detailed in 
                Rufibach (2006, 2007), robustif = FALSE does not. In the latter case, convergence of the algorithm
                is no longer guaranteed.}
  \item{print}{print = TRUE outputs log-likelihood in every loop, print = FALSE does not. Make sure to tell R to output (press CTRL+W).}
}
\value{
  \item{x}{Column vector with observations \eqn{x_1, \ldots, x_m.}}
  \item{f}{Column vector with entries \eqn{\widehat f_m(x_i).}}
  \item{Loglik}{The value \eqn{L(\widehat \varphi_m)} of the log-likelihood-function \eqn{L} at the maximum \eqn{\widehat \varphi_m.}}
  \item{Iterations}{Number of iterations performed.}
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
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}
\seealso{The following functions are used by \code{\link{icmaLogCon}}: 

\code{\link{phieta}}, \code{\link{etaphi}}, \code{\link{Lhat_eta}}, \code{\link{quadDeriv}},
\code{\link{robust}}, \code{\link{isoMean}}

Log concave density estimation via an Active Set Algorithm can be performed using
\code{\link{activeSetLogCon}}.
}
\examples{
set.seed(1977)
x <- sort(rgamma(200, 2, 1))
res <- icmaLogCon(x, w = NA, T1 = 2000, robustif = TRUE, print = TRUE)

## plot resulting functions
par(mfrow = c(2, 1), mar = c(3, 2, 1, 2))
plot(x, res$f, type = 'l'); rug(x)
plot(x, log(res$f), type = 'l'); rug(x)
}
\keyword{htest}
\keyword{nonparametric}

