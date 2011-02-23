logConROC <- function(cases, controls, grid){
        
    m <- length(controls) 
    n <- length(cases)        
        
    res0 <- logConDens(controls, smoothed = TRUE, print = FALSE)
    res1 <- logConDens(cases, smoothed = TRUE, print = FALSE)

    qs0 <- quantilesLogConDens(1 - grid, res0)
    F1 <- evaluateLogConDens(qs0[, "quantile"], res1, which = c(3, 5), gam = NULL, print = FALSE)

    rocLogCon <- 1 - F1[, "CDF"]
    rocLogConSmooth <- 1 - F1[, "smooth.CDF"]

    ## manually adjust
    rocLogCon[grid == 0] <- 0
    rocLogCon[grid == 1] <- 1
    rocLogConSmooth[grid == 0] <- 0
    rocLogConSmooth[grid == 1] <- 1
    
    ## compute AUC
    auc <- integrate(f = ROCx, lower = 0, upper = 1, res0 = res0, res1 = res1, smooth = FALSE)$value
    aucSmooth <- integrate(f = ROCx, lower = 0, upper = 1, res0 = res0, res1 = res1, smooth = TRUE)$value
    
    return(list("m" = m, "n" = n, "fROC" = rocLogCon, "fROC.smooth" = rocLogConSmooth, "auc" = auc, "aucSmooth" = aucSmooth, "res0" = res0, "res1" = res1))
}







