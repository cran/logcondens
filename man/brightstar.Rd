\name{brightstar}
\alias{brightstar}
\docType{data}
\title{Bright star dataset used to illustrate log-concave density estimation}
\description{Dataset that contains the data analyzed in Mizera and Koenker (2009, Section 5).
The sample consists of measurements of radial and rotational velocities for the stars from the Bright Star Catalog,
see Hoffleit and Warren (1991).}
\usage{data(brightstar)}
\format{
  A data frame with 9092 rows on the following 2 variables.
  \describe{
    \item{\code{nr}}{Location of measurements.}
    \item{\code{rad}}{Measurements of radial velocities.}
    \item{\code{rot}}{Measurements of rotational velocities.}
  }
}
\source{The data was downloaded from \url{ftp://adc.astro.umd.edu/pub/adc/archives/catalogs/5/5050/} on September 3, 2009.}
\references{
Duembgen, L. and Rufibach, K. (2009)
Maximum likelihood estimation of a log--concave density and its distribution function: basic properties and uniform consistency. 
\emph{Bernoulli}, \bold{15(1)}, 40--68.

Hoffleit, D., Warren, W.H. (1991).
\emph{The Bright Star Catalog.}
Yale University Observatory, New Heaven.  

Mizera, I., Koenker, R. (2010).
Quasi-concave density estimation.
\emph{Ann. Statist.}, \bold{38(5)}, 2998--3027. 
}
\examples{
# ---- load rotational velocity data ----
data(brightstar)

# ---- compute and plot log-concave estimate ----
# See also Figure 3 in Koenker & Mizera (2009)
x0 <- sort(brightstar[, 3])
res <- logConDens(x0, print = FALSE, smoothed = FALSE)
plot(res, which = "density")
}
\keyword{datasets}









