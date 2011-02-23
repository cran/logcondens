\name{logConROC}
\alias{logConROC}
\title{Compute ROC curve based on log-concave estimates for the constituent distributions}
\description{The receiver operating characteristic (ROC) curve for two constituent distributions \eqn{F} and \eqn{G} is defined as 

\deqn{R(t; F, G)  =  1 - G(F^{-1}(1 - t))}

for \eqn{t \in [0, 1]}. It is typically used to assess the performance of a diagnostic test used to discriminate between healthy and diseased
individuals based on a continuous variable.
}
\usage{logConROC(cases, controls, grid)}
\arguments{
  \item{cases}{A vector of measurements for the cases.}
  \item{controls}{A vector of measurements for the controls.}
  \item{grid}{A vector specifying the grid where the ROC curve is computed on.}
}
\details{In Rufibach (2011) it was shown that the ROC curve based on log-concave density estimates exhibit nice properties for finite sample sizes
as well as asymptotically. Its performance is typically much better than that of the empirical ROC curve and only, if at all, sligthly worse compared
to the binormal model when in fact the underlying densities are normal. However, log-concavity encompasses many parametric densities, so this
new model is much more flexible than the binormal one, at little efficiency sacrifice.}
\value{
\item{m}{Number of control measurements.}
\item{n}{Number of case measurements.}
\item{fROC}{Estimated ROC curve based on the log-concave density estimate.}
\item{fROC.smooth}{Estimated ROC curve based on the smoothed log-concave density estimate.}
\item{res0}{\code{dlc} object as a result of a call to \code{\link{logConDens}} for the data of the controls.}
\item{res1}{\code{dlc} object as a result of a call to \code{\link{logConDens}} for the data of the cases.}
}
\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68. 

Duembgen, L. and Rufibach, K. (2010)
logcondens: Computations Related to Univariate Log-Concave Density Estimation. 
\emph{Journal of Statistical Software}, to appear.  

Rufibach, K. (2011).
A smooth ROC curve estimator based on log-concave density estimates.
Preprint.
}
\author{
Kaspar Rufibach, \email{kaspar.rufibach@ifspm.uzh.ch}, \cr \url{http://www.biostat.uzh.ch/aboutus/people/rufibach.html}}
\note{Note.
}
\seealso{For the computation of the AUC the function \code{\link{ROCx}} is used. In the example below we analyze the 
\code{\link{pancreas}} data.}
\examples{
## ROC curve for pancreas data 
data(pancreas)
status <- factor(pancreas[, "status"], levels = 0:1, 
    labels = c("healthy", "diseased"))
var <- log(pancreas[, "ca199"])
cases <- var[status == "diseased"]
controls <- var[status == "healthy"]

## compute and plot empirical ROC curve
## code modified from https://stat.ethz.ch/pipermail/r-help/2008-October/178531.html
xx <- c(-Inf, sort(unique(c(cases, controls))), Inf)
sens <- sapply(xx, function(x){mean(cases >= x)})  
spec <- sapply(xx, function(x){mean(controls < x)})

## compute log-concave ROC curve
grid <- seq(0, 1, by = 1 / 500)
roc.logcon <- logConROC(cases, controls, grid)        

## plot
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'l', 
    main = "ROC curves for pancreas data", xlab = "1 - specificity", 
    ylab = "sensitivity", pty = 's')
legend("bottomright", c("empirical ROC", "log-concave ROC", "smooth log-concave ROC"), 
    lty = c(1, 1, 2), lwd = 2, col = 2:4, bty = "n")
segments(0, 0, 1, 1, col = 1)
lines(1 - spec, sens, type = 'l', col = 2, lwd = 2)
lines(grid, roc.logcon$fROC, col = 3, lwd = 2)
lines(grid, roc.logcon$fROC.smooth, col = 4, lwd = 2, lty = 2)
}
\keyword{nonparametric}