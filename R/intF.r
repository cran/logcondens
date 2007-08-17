"intF" <-
function (s, x, phi) 
{
    if (min(s) < min(x) | max(s) > max(x)) {
        cat("All elements of s must be in [x_1, x_n]!")
    }
    else {
        n <- length(x)
        dx <- c(NA, diff(x))
        dphi <- c(NA, diff(phi))
        f <- exp(phi)
        F <- res$F
        intF.xi <- c(0, rep(NA, n - 1))
        for (i in 2:n) {
            intF.xi[i] <- dx[i] * (F[i - 1] + dx[i]/dphi[i] * 
                (J00(phi[i - 1], phi[i], 1) - f[i - 1]))
        }
        intF.xi <- cumsum(intF.xi)
        intF.s <- rep(NA, length(s))
        for (k in 1:length(s)) {
            j <- max((1:n)[x <= s[k]])
            j <- min(j, n - 1)
            xj <- x[j]
            Fj <- F[j]
            intF.s[k] <- intF.xi[j] + (s[k] - xj) * F[j] + 
                dx[j + 1] * (dx[j + 1]/dphi[j + 1] * J00(phi[j], 
                  phi[j + 1], (s[k] - xj)/(dx[j + 1])) - (s[k] - 
                  xj)/(dphi[j + 1]) * f[j])
        }
        return(intF.s)
    }
}

