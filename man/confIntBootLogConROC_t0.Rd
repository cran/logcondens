\name{confIntBootLogConROC_t0}
\alias{confIntBootLogConROC_t0}
\title{Function to compute a bootstrap confidence interval for the ROC curve at a given t, based on the log-concave ROC curve}
\description{This function computes a bootstrap confidence interval for the ROC curve at a given value false negative fraction (1 - specificity) \eqn{t}. The ROC curve estimate is based on log-concave densities, as discussed in Rufibach (2011).}
\usage{confIntBootLogConROC_t0(controls, cases, grid = c(0.2, 0.8), conf.level = 0.95, 
M = 1000, smooth = TRUE, output = TRUE)}
\arguments{
  \item{cases}{Values of the continuous variable for the cases.}
  \item{controls}{Values of the continuous variable for the controls.}
  \item{grid}{Values of 1 - specificity where confidence intervals should be computed at (may be a vector).}
  \item{conf.level}{Confidence level of confidence interval.}
  \item{M}{Number of bootstrap replicates.}
  \item{smooth}{\code{Logical}. Compute confidence interval also for ROC curve estimate based on smoothed log-concave densities.}
  \item{output}{\code{Logical}. Show progress of computations?}
}
\value{A list containing the following elements:
  \item{qs}{\code{data.frame} with the columns \code{t} (false positive fractions where confidence interval is computed at) and the confidence intervals for the ROC curve at \code{grid}, based on the log-concave density estimate.}
  \item{boot.mat}{Bootstrap samples for the ROC curve based on the log-concave density estimate.}
  \item{qs.smooth}{If \code{smooth = TRUE}, same as \code{qs} but for the ROC curve based on the smooth log-concave density estimate.}
  \item{boot.mat.smooth}{If \code{smooth = TRUE}, bootstrap samples for the ROC curve based on the smoothed log-concave density estimate.}
}

\references{The reference for computation of these bootstrap confidence intervals is:

Rufibach, K. (2012).
A smooth ROC curve estimator based on log-concave density estimates.
\emph{Int. J. Biostat.}, \bold{8}(1), 1--29.

The bootstrap competitor based on the empirical ROC curve is described in:

Zhou, X.H. and Qin, G. (2005).
Improved confidence intervals for the sensitivity at a fixed level of specificity of a continuous-scale diagnostic test.
\emph{Statist. Med.}, \bold{24}, 465--477.
}

\author{Kaspar Rufibach (maintainer) \cr \email{kaspar.rufibach@gmail.com} \cr \url{http://www.kasparrufibach.ch}.}

\note{The confidence intervals are only valid if observations are \emph{independent}, i.e. eacht patient only contributes one measurement, e.g.}

\seealso{The ROC curve based on log-concave density estimates can be computed using \code{\link{logConROC}}. In the example below we analyze the \code{\link{pancreas}} data.}

\examples{
\dontrun{
## ROC curve for pancreas data 
data(pancreas)
status <- factor(pancreas[, "status"], levels = 0:1, labels = c("healthy", "diseased"))
var <- log(pancreas[, "ca199"])
cases <- var[status == "diseased"]
controls <- var[status == "healthy"]

## compute confidence intervals
res <- confIntBootLogConROC_t0(controls, cases, grid = c(0.2, 0.8), conf.level = 0.95, 
    M = 1000, smooth = TRUE, output = TRUE)
res
}
}

\keyword{htest}

