"icmaLogCon" <-
function (x, w = NA, eps = 10^-8, T1 = 2000, robust = TRUE, print = FALSE) 
{
    n <- length(x)
    if (sum(x[2:n] <= x[1:n - 1]) > 0) {
        cat("We need strictly increasing numbers x(i)!\n")
    }
    if (max(is.na(w)) == 1) {
        w <- rep(1/n, n)
    }
    if (sum(w <= 0) > 0) {
        cat("We need strictly positive weights w(i)!\n")
    }
    w <- w/sum(w)
    dx <- c(0, diff(x))
    iter1 <- 0
    dirder <- 2 * eps
    phi <- LocalNormalize(x, 1:n * 0)
    phi <- LocalMLE(x, w, c(1, rep(0, n - 2), 1), phi, eps)$phi
    eta <- phieta(x, phi)
    loglik <- Lhat_eta(x, w, eta)$ll
    etanew <- 1:n * 0
    while (abs(dirder) > eps && iter1 < T1) {
        iter1 <- iter1 + 1
        derivs <- quadDeriv(dx, w, eta)
        grad <- derivs[, 1]
        hess <- -derivs[, 2]
        y <- eta + grad/hess
        etanew[1] <- y[1]
        etanew[2:n] <- -isoMean(-y[2:n], hess[2:n])
        if (robust == TRUE) {
            etanew <- robust(x, w, eta, etanew, grad)
        }
        dirder <- as.numeric(t(grad) %*% (etanew - eta))
        eta <- etanew
        loglik <- Lhat_eta(x, w, eta)$ll
        out <- data.frame(Iteration = iter1, LogLikelihood = loglik, 
            dirder = dirder)
        if (print == TRUE) {
            print(out)
        }
    }
    return(list(x = x, f = exp(etaphi(x, eta)), Loglik = loglik, 
        Iterations = iter1))
}
