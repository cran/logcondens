\name{summary.dlc}
\Rdversion{1.1}
\alias{summary.dlc}
\title{Summarizing log-concave density estimation}
\description{\code{summary} method for class \code{"dlc"}.}
\usage{\method{summary}{dlc}(object, \dots)}
\arguments{
  \item{object}{An object of class \code{"dlc"}, usually a result of a call to \code{logConDens}.}
  \item{\dots}{Further arguments.}
}
\details{See \code{\link{activeSetLogCon}} and \code{\link{evaluateLogConDens}} for details on the computations.}
\value{The function \code{summary.dlc} returns a list of summary statistics of the estimated
log-concave density as well as of its smoothed version (depending on the value of \code{smoothed} when calling
\code{\link{logConDens}}).}
\section{Warning}{Note that the numbering of knots in the output relies on the vector of \emph{unique} observations.}
\references{
Duembgen, L, Huesler, A. and Rufibach, K. (2010).
Active set and EM algorithms for log-concave densities based on complete and censored data. 
Technical report 61, IMSV, Univ. of Bern, available at \url{https://arxiv.org/abs/0707.4643}.

Duembgen, L. and Rufibach, K. (2009).
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \doi{https://doi.org/10.18637/jss.v039.i06}
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{https://www.imsv.unibe.ch/about_us/staff/prof_dr_duembgen_lutz/index_eng.html}
}
\examples{
## See help file of function "logConDens".
}
\keyword{htest}
\keyword{nonparametric}
