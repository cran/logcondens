"qloglin" <-
function (u, t) 
{
    m <- length(u)
    z <- 1:m * NA
    II <- (1:m)[abs(t) > 10^(-4)]
    z[II] <- log( 1 + ( (exp(t) - 1 ) * u[II])) / t
    II <- (1:m)[abs(t) <= 10^(-4)]
    z[II] <- u[II] + t * u[II]*(1-u[II])/2
    return(matrix(z, ncol = 1))
}
