rlogcon <-
function (n, x0) 
{
    n0 <- length(x0)
    x <- unique(x0)
    w <- as.vector(table(x0))
    w <- w/n0
    res <- activeSetLogCon(x, w, print = FALSE)
    phi <- as.vector(res$phi[, 1])
    f.smoothed <- logconSmoothed(x, w, phi, gam = NULL, sig = sd(x0))
    Fj <- res$Fhat[res$IsKnot == 1]
    m <- length(Fj[-1])
    J <- 1:m
    dist_J <- diff(Fj)
    x_s <- res$x[res$IsKnot == 1][-1]
    phi_s <- phi[res$IsKnot == 1][-1]
    X <- rep(NA, length = n)
    J <- sample(x = 1:m, size = n, prob = dist_J, replace = TRUE)
    U <- runif(n)
    Z <- rnorm(n)
    for (i in 1:n) {
        Ji <- J[i]
        xj <- x_s[Ji]
        if (Ji == 1) {
            xj1 <- x[1]
            phij1 <- phi[1]
        }
        else {
            xj1 <- x_s[Ji - 1]
            phij1 <- phi_s[Ji - 1]
        }
        theta <- phi[Ji] - phij1
        X[i] <- xj1 + (x_s[Ji] - xj1) * ((theta == 0) * U[i] + 
            (theta != 0) * qloglin(U[i], theta))
    }
    X_star <- X + f.smoothed$gam * Z
    res <- list(X = X, X_star = X_star, U = U, Z = Z, f = exp(phi), 
        f.smoothed = f.smoothed, x = x, w = w)
    return(res)
}
