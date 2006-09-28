"J00" <-
function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 10^(-6)]
    z[II] <- z[II] * (exp(d[II]) - 1)/d[II]
    II <- (1:m)[abs(d) <= 10^(-6)]
    z[II] <- z[II] * (1 + d[II]/2 + d[II]^2/6)
    return(matrix(z, ncol = 1))
}
