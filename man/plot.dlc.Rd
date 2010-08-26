\name{plot.dlc}
\alias{plot.dlc}
\title{Standard plots for a dlc object}
\description{\code{plot} method for class \code{"dlc"}.
Three plots (selectable by \code{which}) are currently available: a plot of the estimated density,
the estimated log-density, or the distribution function corresponding to the estimated log-concave density.
By default, a plot of the density estimate is provided. If \code{smoothed = TRUE}, the smoothed version of 
the log-concave density estimate (saved in \code{x}) is added to the density and log-density plot.
For the CDF, the smoothed version is not contained by default in a \code{dlc} object and needs to be computed
when asked to be plotted.}
\usage{\method{plot}{dlc}(x, which = c("density", "log-density", "CDF"), 
    add.title = TRUE, legend.pos = "topright", \dots)}
\arguments{
  \item{x}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens}.}
  \item{which}{One of \code{"density"}, \code{"log-density"}, or \code{"CDF"}.}
  \item{add.title}{Logical, if \code{TRUE} adds a standard title to the plot.}
  \item{legend.pos}{Placement of the legend. One of \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, \code{"left"}, 
  \code{"topleft"}, \code{"top"}, \code{"topright"}, \code{"right"}, \code{"center"}; or \code{"none"} for not displaying a legend. See \code{\link{legend}} for details.}
  \item{\dots}{Further arguments.}  
}
\details{See \code{\link{activeSetLogCon}} and \code{\link{evaluateLogConDens}} for details on the computations.}
\value{Chosen plot is generated.}
\references{
Duembgen, L, Huesler, A. and Rufibach, K. (2007).
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Technical report 61, IMSV, Univ. of Bern, available at \url{http://arxiv.org/abs/0707.4643}.

Duembgen, L. and Rufibach, K. (2009).
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.staff.unibe.ch/duembgen}
}
\examples{
## See help file of function "logConDens".
}
\keyword{htest}
\keyword{nonparametric}






