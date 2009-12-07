evaluateLogConDens <-
function (x0, x, phi, Fhat, IsKnot) 
{
    n <- length(x)
    if (x0 < x[1]) {
        res <- c(-Inf, 0, 0)
    }
    if (x0 == x[1]) {
        res <- c(phi[1], exp(phi[1]), 0)
    }
    if (x0 > x[1] && x0 <= x[n]) {
        res <- NA
        x.knot <- x[IsKnot > 0]
        phi.knot <- phi[IsKnot > 0]
        k <- length(x.knot[x.knot < x0])
        phi.x0 <- (1 - (x0 - x.knot[k])/(x.knot[k + 1] - x.knot[k])) * 
            phi.knot[k] + (x0 - x.knot[k])/(x.knot[k + 1] - x.knot[k]) * 
            phi.knot[k + 1]
        f.x0 <- exp(phi.x0)
        j <- length(x[x < x0])
        Fhat.x0 <- Fhat[j] + (x[j + 1] - x[j]) * J00(phi[j], 
            phi[j + 1], (x0 - x[j])/(x[j + 1] - x[j]))
        res <- c(phi.x0, f.x0, Fhat.x0)
    }
    if (x0 > x[n]) {
        res <- c(-Inf, 0, 1)
    }
    return(res)
}
