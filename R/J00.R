`J00` <-
function (x, y, v = 1) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 10^(-6)]
    z[II] <- z[II] * (exp(v * d[II]) - 1)/d[II]
    II <- (1:m)[abs(d) <= 10^(-6)]
    z[II] <- z[II] * (v + v^2 * d[II]/2 + v^3 * d[II]^2/6)
    return(matrix(z, ncol = 1))
}
