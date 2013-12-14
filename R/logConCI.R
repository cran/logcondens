logConCI <- function(res, xx0, conf.level = c(0.8, 0.9, 0.95, 0.99)[3], type = c("DR", "ks", "nrd")[2], htype = c("hscv", "hlscv", "hpi", "hns")[4]){

                ## re-translate to Hanna's arguments:
                data        <-  res$xn
                mle         <-  res
                alpha       <-  1 - conf.level

                ## compute the CIs
				fmle		<-	evaluateLogConDens(xx0, mle, which = 2)
				fhat		<-	fmle[, 3]
				n		<-	length(data)					
				c2fit		<-	c2hat(res = mle, xx0 = xx0, type = type, htype = htype)

				if(conf.level == 0.99){
				zup			<-	3.6881	
				zlo			<-	-3.0905
				}
				if(conf.level == 0.95){
				zup			<-	2.7536	
				zlo			<-	-2.4157
				}
				if(conf.level == 0.9){
				zup			<-	2.2653	
				zlo			<-	-2.0574
				}
				if(conf.level == 0.80){
				zup			<-	1.7421	
				zlo			<-	-1.6440
				}
	
				up_DR		<-	fhat - zlo * (n ^ (-0.4)) * c2fit$DR
				up_ks_hscv	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hscv
				up_ks_hlscv	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hlscv
				up_ks_hpi	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hpi
				up_ks_hns	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$ks_hns
				up_nrd	<-	fhat - zlo * (n ^ (-0.4)) * c2fit$nrd
	
				lo_DR       <-    fhat - zup * (n ^ (-0.4)) * c2fit$DR
				lo_DR		<-	if (length(lo_DR) > 0){lo_DR <- pmax(0, lo_DR)}

				lo_ks_hscv	<-	fhat - zup * (n ^ (-0.4)) * c2fit$ks_hscv
				lo_ks_hscv	<-	if (length(lo_ks_hscv) > 0){lo_ks_hscv <- pmax(0, lo_ks_hscv)}

				lo_ks_hlscv	<-	fhat - zup * (n ^ (-0.4)) * c2fit$ks_hlscv
				lo_ks_hlscv	<-	if (length(lo_ks_hlscv) > 0){lo_ks_hlscv <- pmax(0, lo_ks_hlscv)}

				lo_ks_hpi	<-	fhat - zup * (n ^ (-0.4)) * c2fit$ks_hpi
				lo_ks_hpi	<-	if (length(lo_ks_hpi) > 0){lo_ks_hpi <- pmax(0, lo_ks_hpi)}

				lo_ks_hns	<-	fhat - zup * (n ^ (-0.4)) * c2fit$ks_hns
				lo_ks_hns	<-	if (length(lo_ks_hns) > 0){lo_ks_hns <- pmax(0, lo_ks_hns)}

				lo_nrd	<-	fhat - zup * (n ^ (-0.4)) * c2fit$nrd
				lo_nrd	<-	if (length(lo_nrd) > 0){lo_nrd <- pmax(0, lo_nrd)}
	
	            	res <- list(fhat = fhat, up_DR = up_DR, lo_DR = lo_DR, up_ks_hscv = up_ks_hscv, lo_ks_hscv = lo_ks_hscv, up_ks_hlscv = up_ks_hlscv, lo_ks_hlscv = lo_ks_hlscv, up_ks_hpi = up_ks_hpi, lo_ks_hpi = lo_ks_hpi, up_ks_hns = up_ks_hns, lo_ks_hns = lo_ks_hns, up_nrd = up_nrd, lo_nrd = lo_nrd)
				return(res)
}
