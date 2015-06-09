\name{MLE}
\alias{MLE}
\title{Unconstrained piecewise linear MLE}
\description{
Given a vector of observations \eqn{{\bold{x}} = (x_1, \ldots, x_m)}{x = (x_1, \ldots, x_m)} with pairwise distinct entries and
a vector of weights \eqn{{\bold{w}}=(w_1, \ldots, w_m)}{w =(w_1, \ldots, w_m)} s.t. \eqn{\sum_{i=1}^m w_i = 1}, this function computes a function \eqn{\widehat \phi_{MLE}}{\hat \phi_{MLE}} (represented by the vector \eqn{(\widehat \phi_{MLE}(x_i))_{i=1}^m}{(\hat \phi_{MLE}(x_i))_{i=1}^m}) supported by \eqn{[x_1, x_m]} such that 

\deqn{L(\phi) = \sum_{i=1}^m w_i \phi(x_i) - \sum_{j=1}^{m-1} (x_{j+1} - x_j) J(\phi_j, \phi_{j+1})}

is maximal over all continuous, piecewise linear functions with knots in \eqn{\{x_1, \ldots, x_m\}}{{x_1, \ldots, x_m}.}
}
\usage{MLE(x, w = NA, phi_o = NA, prec = 1e-7, print = FALSE)}
\arguments{
  \item{x}{Vector of independent and identically distributed numbers, with strictly increasing entries.}
  \item{w}{Optional vector of nonnegative weights corresponding to \eqn{{\bold{x}_m}}{x_m}.}
  \item{phi_o}{Optional starting vector.}
  \item{prec}{Threshold for the directional derivative during the Newton-Raphson procedure.}
  \item{print}{print = TRUE outputs log-likelihood in every loop, print = FALSE does not. Make sure to tell R to output (press CTRL+W).}
}
\value{
  \item{phi}{Resulting column vector \eqn{(\widehat \phi_{MLE}(x_i))_{i=1}^m.}{(\hat \phi_{MLE}(x_i))_{i=1}^m.}}
  \item{L}{Value \eqn{L(\widehat \phi_{MLE})}{L(\hat \phi_{MLE})} of the log-likelihood at \eqn{\widehat \phi_{MLE}.}{\hat \phi_{MLE}.}}
  \item{Fhat}{Vector of the same length as \eqn{{\bold{x}}}{x} with entries \eqn{\widehat F_{MLE,1} = 0}{\hat F_{MLE,1} = 0} and 

\deqn{\widehat F_{MLE,k} = \sum_{j=1}^{k-1} (x_{j+1} - x_j) J(\phi_j, \phi_{j+1}) }{\hat F_{MLE,k} = \sum_{j=1}^{k-1} (x_{j+1} - x_j) J(\phi_j, \phi_{j+1})}

for \eqn{k \ge 2.}{k >= 2.}}
}
\note{This function is not intended to be invoked by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}}
\keyword{htest}
\keyword{nonparametric}


