\name{logconTwoSample}
\alias{logconTwoSample}
\title{Compute p-values for two-sample test based on log-concave CDF estimates}
\description{Compute \eqn{p}-values for a test for the null hypothesis of equal CDFs of two samples. The test
statistic is reminiscient of Kolmogorv-Smirnov's, but instead of computing it for the empirical CDFs, this function
computes it based on log-concave estimates for the CDFs.}
\usage{logconTwoSample(x, y, which = c("MLE", "smooth"), M = 999, 
    n.grid = 500, display = TRUE, seed0 = 1977)}
\arguments{
  \item{x}{First data sample.}
  \item{y}{Second data sample.}
  \item{which}{Indicate for which type of estimate the test statistic should be computed.}
  \item{M}{Number of permutations.}
  \item{n.grid}{Number of grid points in computation of maximal difference between smoothed log-concave CDFs. See \code{\link{maxDiffCDF}} for details.}
  \item{display}{If \code{TRUE} progress of computations is shown.}
  \item{seed0}{Set seed to reproduce results.}
}
\details{Given two i.i.d. samples \eqn{x_1, \ldots, x_{n_1}} and \eqn{y_1, \ldots, y_{n_2}} this function computes a permutation
test \eqn{p}-value that provides evidence against the null hypothesis

\deqn{H_0 : F_1 = F_2}

where \eqn{F_1, F_2} are the CDFs of the samples, respectively. A test either based on the log-concave MLE or on its
smoothed version (see Duembgen and Rufibach, 2009, Section 3) are provided. Note that computation of the smoothed
version takes considerably more time.
}
\value{
\item{p.value}{A two dimensional vector containing the \eqn{p}-values.}
\item{test.stat.orig}{The test statistics for the original samples.}
\item{test.stats}{A \eqn{M \times 2} matrix containing the test statistics for all the permutations.}
}
\section{Warning}{Note that the algorithm that finds the maximal difference for the smoothed estimate is of approximative
nature only. It may fail for very large sample sizes.}
\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.  

Duembgen, L. and Rufibach, K. (2011)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, \bold{39(6)}, 1--28. \url{http://www.jstatsoft.org/v39/i06}  
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@gmail.com}, \cr \url{http://www.kasparrufibach.ch} 

Lutz Duembgen, \email{duembgen@stat.unibe.ch}, \cr \url{http://www.stat.unibe.ch/content/staff/personalhomepages/duembgen/index_eng.html}
}
\examples{
\dontrun{
n1 <- 30
n2 <- 25
x <- rgamma(n1, 2, 1)
y <- rgamma(n2, 2, 1) + 1
twosample <- logconTwoSample(x, y, which = c("MLE", "smooth")[1], M = 999)
}
}
\keyword{htest}
\keyword{nonparametric}
