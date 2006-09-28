"LocalF" <-
function (x, phi) 
{
    n <- length(x)
    F <- 1:n * 0
    dx <- diff(x)
    F[2:n] <- cumsum(dx * J00(phi[1:(n - 1)], phi[2:n]))
    return(matrix(F, ncol = 1))
}
