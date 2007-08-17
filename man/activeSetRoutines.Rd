\name{activeSetRoutines}
\alias{activeSetRoutines}
\alias{LocalCoarsen}
\alias{LocalConvexity}
\alias{LocalExtend}
\alias{LocalF}
\alias{LocalNormalize}
\alias{LocalMLE}
\title{Auxiliary Numerical Routines for the Function activeSetLogCon}
\description{Functions that are used by activeSetLogCon.}
\usage{
LocalCoarsen(x, w, IsKnot)
LocalConvexity(x, phi)
LocalExtend(x, IsKnot, x2, phi2) 
LocalF(x, phi)
LocalNormalize(x, phi)
LocalMLE(x, w, IsKnot, phi_o, prec)
}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}}}, where \eqn{w_1 > 0} and \eqn{w_m > 0}. These raw weights are 
            normalized in order to sum to one. Default: \eqn{w_i = 1 / m}.}
  \item{IsKnot}{Column vector with entries \cr
  IsKnot\eqn{_i = 1\{\varphi} has a kink at \eqn{x_i\}.}}
  \item{phi}{Column vector with entries \eqn{\varphi(x_i).}}
  \item{x2}{Vector of same type as \eqn{{\bold{x}}}.}
  \item{phi2}{Vector of same type as \eqn{{\bold{\varphi}}}.}
  \item{phi_o}{Optional starting vector.}
  \item{prec}{Threshold for the directional derivative during Newton-Raphson procedure.}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@stanford.edu}, \cr \url{http://www.stanford.edu/~kasparr} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/~duembgen}}

\seealso{All the above functions are used by \code{\link{activeSetLogCon}} to estimate a log-concave
probability density.

Log concave density estimation via an iterative convex minorant algorithm can be performed using
\code{\link{icmaLogCon}}.
}
\keyword{htest}
\keyword{nonparametric}

