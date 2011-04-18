\name{qloglin}
\alias{qloglin}
\title{Quantile Function In a Simple Log-Linear model}
\description{Suppose the random variable \eqn{X} has density function 

  \deqn{g_\theta(x) = \frac{\theta \exp(\theta x)}{\exp(\theta) - 1}}{g_\theta(x) = (\theta exp(\theta x))/(exp(\theta) - 1)}
  
  for an arbitrary real number \eqn{\theta} and \eqn{x \in [0,1]}. The function \code{\link{qloglin}} is simply the 
  quantile function 
  
  \deqn{G^{-1}_\theta(u) = \theta^{-1} \log \Big( 1 + (e^\theta - 1)u \Big)}{G^{-1}_\theta(u) = \theta^{-1} log (1 + (e^\theta - 1)u)}    
  
  in this model, for \eqn{u \in [0,1]}. This quantile function is used for the computation of quantiles of \eqn{\widehat F_m} in \code{\link{quantilesLogConDens}}.} 
\usage{qloglin(u, t)}
\arguments{
  \item{u}{Vector in \eqn{[0,1]^d} where quantiles are to be computed at.}
  \item{t}{Parameter \eqn{\theta}.}
}
\value{\item{z}{Vector containing the quantiles \eqn{G_n^{-1}(u_i)} for \eqn{i = 1, \ldots, d}.}}
\note{Taylor approximation is used if \eqn{\theta} is small.

This function is not intended to be called by the end user.}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}}

\keyword{htest}
\keyword{nonparametric}
