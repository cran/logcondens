Q00 <- function(x, a, u, v, gamma, QFhat = FALSE){

    ## x, gamma in R
    ## a, u, v in R^m
    q <- exp(a * (x - u) + a ^ 2 * gamma ^ 2 / 2) * (pnorm(v - x - a * gamma ^ 2, mean = 0, sd = gamma) - pnorm(u - x - a * gamma ^ 2, mean = 0, sd = gamma))

    Q <- NA
    if (QFhat == TRUE){
        m <- length(a)
        Q <- rep(NA, m)
        
        II <- (1:m)[abs(a) > 10^-6]
        Q[II] <- q[II] / a[II] + (exp(a[II] * (v[II] - u[II])) * pnorm(x - v[II], mean = 0, sd = gamma) - pnorm(x - u[II], mean = 0, sd = gamma)) / a[II]
    
        II <- (1:m)[abs(a) <= 10^-6]
        Q[II] <- (x - u[II]) * pnorm(x - u[II], mean = 0, sd = gamma) - (x - v[II]) * pnorm(x - v[II], mean = 0, sd = gamma) + gamma ^ 2 * (dnorm(x - u[II], mean = 0, sd = gamma) - dnorm(x - v[II], mean = 0, sd = gamma))
    }

    res <- list("q" = q, "Q" = Q)
    return(res)
}
