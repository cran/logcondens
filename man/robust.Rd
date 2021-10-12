\name{robust}
\alias{robust}
\title{Robustification and Hermite Interpolation for ICMA}
\description{Performs robustification and Hermite interpolation in the iterative convex minorant algorithm as described in Rufibach (2006, 2007).}
\usage{robust(x, w, eta, etanew, grad)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}_m}}{x_m}.}
  \item{eta}{Current candidate vector.}
  \item{etanew}{New candidate vector.}
  \item{grad}{Gradient of L at current candidate vector \eqn{\eta.}}
}
\value{Returns a (possibly) new vector \eqn{\eta} on the segment 

\deqn{(1 - t_0)  \eta + t_0  \eta_{new} }

such that the log-likelihood of this new \eqn{\eta} is strictly greater than that of the initial \eqn{\eta} and \eqn{t_0} is chosen 
according to the Hermite interpolation procedure described in Rufibach (2006, 2007).}

\references{
Rufibach K. (2006) \emph{Log-concave Density Estimation and Bump Hunting for i.i.d. Observations.}
PhD Thesis, University of Bern, Switzerland and Georg-August University of Goettingen, Germany, 2006.
\cr Available at \url{https://slsp-ube.primo.exlibrisgroup.com/permalink/41SLSP_UBE/17e6d97/alma99116730175505511}.

Rufibach, K. (2007)
Computing maximum likelihood estimators of a log-concave density function.
\emph{J. Stat. Comput. Simul.} \bold{77}, 561--574.
}

\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}}

\keyword{htest}
\keyword{nonparametric}
