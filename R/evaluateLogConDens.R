"evaluateLogConDens" <-
function (x0, x, f, F, IsKnot) 
{
    n <- length(x)
    f <- as.vector(f)
    phi <- log(f)
    res <- NA
    if (x0 > x[1] && x <= x[n]) {
        knot <- as.vector(IsKnot)
        x.knot <- x[knot > 0]
        phi.knot <- phi[knot > 0]
        f.knot <- f[knot > 0]
        F.knot <- F[knot > 0]
        dx.knot <- c(0, diff(x.knot))
        dp.knot <- c(0, diff(phi.knot))
        i0 <- length(x.knot[x.knot < x0]) + 1
        phix0 <- dp.knot[i0]/dx.knot[i0] * (x0 - x.knot[i0 - 
            1]) + phi.knot[i0 - 1]
        fx0 <- exp(phix0)
        Fx0 <- F.knot[i0 - 1] + dx.knot[i0]/dp.knot[i0] * f.knot[i0 - 
            1] * (exp(dp.knot[i0]/dx.knot[i0] * (x0 - x.knot[i0 - 
            1])) - 1)
        res <- c(phix0, fx0, Fx0)
    }
    if (x0 < min(x)) {
        res <- c(-Inf, 0, 0)
    }
    if (x0 > max(x)) {
        res <- c(-Inf, 0, 1)
    }
    if (x0 == x[1]) {
        res <- c(phi[1], exp(phi[1]), 0)
    }
    return(res)
}
