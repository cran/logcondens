`J11` <-
function (x, y) 
{
    m <- length(x)
    z <- exp(x)
    d <- y - x
    II <- (1:m)[abs(d) > 10^(-3)]
    z[II] <- z[II] * (exp(d[II]) * (d[II] - 2) + 2 + d[II])/(d[II]^3)
    II <- (1:m)[abs(d) <= 10^(-3)]
    z[II] <- z[II] * (1/6 + d[II]/12 + d[II]^2/60 + d[II]^3/180)
    return(matrix(z, ncol = 1))
}
