ggmtrust.path <- function(s, n, data = NULL, nlambda = 10, lambda.min.ratio = 0.1, pen = "lasso", method = "BHHH", w.alasso = NULL, gamma = 1, a = 3.7){


if( !(pen %in% c("ridge", "lasso", "alasso", "scad")) ) stop("pen should be one of: ridge, lasso, alasso, scad")

if( !(method %in% c("H", "BHHH", "TB")) ) stop("method should be one of: H, BHHH, TB")


  lambda.max <- 2*max(abs(s[upper.tri(s, diag = FALSE)]))  # max(abs(s[upper.tri(s, diag = FALSE)]))
  lambda.min <- lambda.min.ratio*lambda.max
  lambda.seq <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))


paramsO <- NULL


modO <- Vb <- Ve <- edf <- params <- omega <- list()

iterations <- countPD <- loglik <- aic <- bic <- t.edf <- NA


for(i in 1:length(lambda.seq)){

ggmtrustOB <- ggmtrust(s, n, data = data, lambda = lambda.seq[i], pen = pen, params = paramsO, method = method, w.alasso = w.alasso, gamma = gamma, a = a)


paramsO <- params[[i]] <- ggmtrustOB$coefficients
paramsO[ggmtrustOB$idx] <- exp(paramsO[ggmtrustOB$idx])

omega[[i]]    <- ggmtrustOB$omega
edf[[i]]      <- ggmtrustOB$edf
t.edf[i]      <- ggmtrustOB$t.edf
aic[i]        <- AIC(ggmtrustOB)
bic[i]        <- BIC(ggmtrustOB)
loglik[i]     <- logLik(ggmtrustOB)
countPD[i]    <- ggmtrustOB$countPD
iterations[i] <- ggmtrustOB$iter.if
Vb[[i]]       <- ggmtrustOB$Vb
modO[[i]]     <- ggmtrustOB

print(i)

}


# external function
#
# 1 estimates original scale, insert alpha prob level, sim.rep
# 2 2.5%
# 3 97.5%
# 4 empirical distribution, hist


L <- list(modO = modO, lambda.seq = lambda.seq, params = params, omega = omega, edf = edf, t.edf = t.edf, 
          aic = aic, bic = bic, loglik = loglik, countPD = countPD, iterations = iterations, Vb = Vb, idx = ggmtrustOB$idx)   


L

}










