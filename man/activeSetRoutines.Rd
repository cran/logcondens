\name{activeSetRoutines}
\alias{activeSetRoutines}
\alias{LocalCoarsen}
\alias{LocalConvexity}
\alias{LocalExtend}
\alias{LocalF}
\alias{LocalNormalize}
\alias{LocalMLE}
\alias{LocalVariance}
\title{Auxiliary Numerical Routines for the Function activeSetLogCon}
\description{Functions that are used by activeSetLogCon.}
\usage{
LocalCoarsen(x, w, IsKnot)
LocalConvexity(x, phi)
LocalExtend(x, IsKnot, x2, phi2) 
LocalF(x, phi)
LocalNormalize(x, phi)
LocalMLE(x, w, IsKnot, phi_o, prec)
LocalVariance(x, w = NULL, phi)
}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}}}{x}.}
  \item{IsKnot}{Vector with entries IsKnot\eqn{_i = 1\{\phi} has a kink at \eqn{x_i\}.}}
  \item{phi}{Vector with entries \eqn{\phi(x_i).}}
  \item{x2}{Vector of same type as \eqn{{\bold{x}}}{x}.}
  \item{phi2}{Vector of same type as \eqn{{\bold{\phi}}}{\phi}.}
  \item{phi_o}{Optional starting vector.}
  \item{prec}{Threshold for the directional derivative during Newton-Raphson procedure.}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.imsv.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}}

\seealso{All the above functions are used by \code{\link{activeSetLogCon}} to estimate a log-concave
probability density.

Log concave density estimation via an iterative convex minorant algorithm can be performed using
\code{\link{icmaLogCon}}.
}
\keyword{htest}
\keyword{nonparametric}

