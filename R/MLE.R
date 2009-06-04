`MLE` <-
function (x, w = NA, phi_o = NA, prec = 10^(-7), print = FALSE) 
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
    ww <- w/sum(w)
    if (max(is.na(phi_o)) == 1) {
        m <- sum(ww * x)
        s2 <- sum(ww * (x - m)^2)
        phi <- LocalNormalize(x, -(x - m)^2/(2 * s2))
    }
    else {
        phi <- LocalNormalize(x, phi_o)
    }
    iter0 <- 0
    res <- Local_LL_all(x, ww, phi)
    L <- res$ll
    phi_new <- res$phi_new
    dirderiv <- res$dirderiv
    while ((dirderiv >= prec) & (iter0 < 100)) {
        iter0 <- iter0 + 1
        L_new <- Local_LL(x, ww, phi_new)
        iter1 <- 0
        while ((L_new < L) & (iter1 < 20)) {
            iter1 <- iter1 + 1
            phi_new <- 0.5 * (phi + phi_new)
            L_new <- Local_LL(x, ww, phi_new)
            dirderiv <- 0.5 * dirderiv
            if (print == TRUE) {
                print(paste("iter0=", iter0, " / iter1=", iter1, 
                  " / L=", round(L_new, 4), " / dirderiv=", round(dirderiv, 
                    4), sep = ""))
            }
        }
        if (L_new >= L) {
            tstar <- max((L_new - L)/dirderiv)
            if (tstar >= 0.5) {
                phi <- LocalNormalize(x, phi_new)
            }
            else {
                tstar <- max(0.5/(1 - tstar))
                phi <- LocalNormalize(x, (1 - tstar) * phi + 
                  tstar * phi_new)
            }
            res <- Local_LL_all(x, ww, phi)
            L <- res$ll
            phi_new <- res$phi_new
            dirderiv <- res$dirderiv
        }
        else {
            dirderiv <- 0
        }
        if (print == TRUE) {
            print(paste("iter0=", iter0, " / iter1=", iter1, 
                " / L=", round(L, 4), " / dirderiv=", round(dirderiv, 
                  4), sep = ""))
        }
    }
    res <- list(phi = matrix(phi, ncol = 1), L = L, Fhat = matrix(LocalF(x, 
        phi), ncol = 1))
    return(res)
}
