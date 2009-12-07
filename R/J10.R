J10 <-
function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 10^(-4)]
    z[II] <- z[II] * (exp(d[II]) - 1 - d[II])/(d[II]^2)
    II <- (1:m)[abs(d) <= 10^(-4)]
    z[II] <- z[II] * (1/2 + d[II]/6 + d[II]^2/24)
    return(matrix(z, ncol = 1))
}
