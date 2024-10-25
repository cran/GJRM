k.tau <- function(x, prob.lev = 0.05){

if(x$Model %in% c("T","TSS","TESS","BPO0","ROY")) stop("Function not extended yet to the chosen model.")

kt <- theta2tau(x$VC$BivD, x$VC$nCa, x$theta, tau.res = TRUE)$tau

if(length(x$theta) == 1) names(kt) <- "tau"
if(length(x$theta) > 1) {kt <- as.data.frame(kt); dimnames(kt)[[2]] <- "tau"} 

ktS <- theta2tau(x$VC$BivD, x$VC$nCa, summary(x)$est.RHOb, tau.res = TRUE)$tau

#if(length(x$theta) > 1) ktS <- t(as.matrix(ktS))

CIktS <- rowQuantiles(ktS, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

#if(length(x$theta) > 1) CIktS <- t(CIktS)


if(length(x$theta) == 1) res <- c(kt, CIktS)
if(length(x$theta) > 1)  res <- cbind(kt, CIktS)

#class(res) <- "tau"

res

}




