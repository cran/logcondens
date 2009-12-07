logconSmoothed <-
function (x, w, phi, gam = NULL, xs = NULL, sig = NA) 
{
    if ((is.null(gam) == TRUE) & (is.na(sig) == TRUE)) {
        cat("Please provide sig!\n")
    }
    n <- length(x)
    f <- exp(phi)
    if (is.null(gam) == TRUE) {
        k <- 1:(n - 1)
        m.x0 <- weighted.mean(x, w)
        VarFn <- sum(w * (x^2)) - sum(diff(x)[k]^3 * J11(x = phi[k], 
            y = phi[k + 1])) - 2 * m.x0 * sum(w * x) + m.x0^2
        gam <- sqrt(sig^2 - VarFn)
    }
    meanX <- weighted.mean(x, w)
    s <- c(NA, diff(phi)/diff(x))
    if (identical(xs, NULL)) {
        r <- diff(range(x))
        xs <- seq(min(x) - 0.05 * r, max(x) + 0.05 * r, length = 500)
    }
    js <- 1:(n - 1)
    f.smoothed <- rep(NA, length(xs))
    for (i in 1:length(xs)) {
        xi <- xs[i]
        f.smoothed[i] <- sum(f[js] * exp(s[js + 1] * (gam^2/2 * 
            s[js + 1] + (xi - x[js]))) * (pnorm(x[js + 1], mean = xi + 
            gam^2 * s[js + 1], sd = gam) - pnorm(x[js], mean = xi + 
            gam^2 * s[js + 1], sd = gam)))
    }
    res <- list(f.smoothed = f.smoothed, gam = gam, xs = xs)
    return(res)
}
