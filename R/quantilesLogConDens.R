"quantilesLogConDens" <-
function (p0, x, f, F, IsKnot) 
{
    n <- length(x)
    if (p0 <= 0) {
        q <- -Inf
    }
    if (p0 >= 1) {
        q <- x[n]
    }
    if (p0 > 0 && p0 < 1) {
        n <- length(x)
        phi <- log(f)
        knot <- as.vector(IsKnot)
        x.knot <- x[knot > 0]
        phi.knot <- phi[knot > 0]
        f.knot <- f[knot > 0]
        F.knot <- F[knot > 0]
        dx.knot <- c(0, diff(x.knot))
        dp.knot <- c(0, diff(phi.knot))
        x0 <- max(x.knot[F.knot <= p0])
        i0 <- min(length(x.knot[x.knot <= x0]) + 1, n)
        q <- x.knot[i0 - 1] + (dx.knot[i0]/dp.knot[i0]) * (log(f.knot[i0] + 
            (p0 - F.knot[i0]) * dp.knot[i0]/dx.knot[i0]) - phi.knot[i0 - 
            1])
    }
    return(q)
}
