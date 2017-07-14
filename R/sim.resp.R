sim.resp <- function(margin, rsim, eta, sigma2, nu, setseed = TRUE){

if(margin == "ZTP") rZTP <- function(n, mu) qpois(runif(n, dpois(0, mu), 1), mu)

if(setseed == TRUE) set.seed(1)

if(margin == "N")     y <- rNO(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "N2")    y <- rNO(   rsim,    mu =     eta,     sigma =      sigma2) 
if(margin == "GU")    y <- rGU(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "rGU")   y <- rRG(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "LO")    y <- rLO(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "LN")    y <- rLOGNO(rsim,    mu =     eta,     sigma = sqrt(sigma2)) 
if(margin == "WEI")   y <- rWEI(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
#if(margin == "GO")    y <- rgompertz(rsim, shape = exp(eta), rate = sqrt(sigma2))
if(margin == "iG")    y <- rIG(   rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "GA")    y <- rGA(   rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "GAi")   y <- rGA(   rsim,    mu =     eta,     sigma = sqrt(sigma2)) 

#if(margin == "GA2")   y <- rgamma(rsim, shape = sqrt(sigma2), rate = exp(eta) ) 
#if(margin == "GGA")   y <- exp( log(rgamma(rsim, shape = nu))/sqrt(sigma2) + log(exp(eta)) )

if(margin == "DAGUM") y <- rGB2(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2), nu = nu, tau = 1) 
if(margin == "SM")    y <- rGB2(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2), nu = 1 , tau = nu) 
if(margin == "BE")    y <- rBE(   rsim,    mu = plogis(eta), sigma = sqrt(sigma2)) 
if(margin == "FISK")  y <- rGB2(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2), nu = 1 , tau = 1)
if(margin == "NBI")   y <- rNBI(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "NBII")  y <- rNBII( rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "PIG")   y <- rPIG(  rsim,    mu = exp(eta),    sigma = sqrt(sigma2)) 
if(margin == "PO")    y <- rPO(   rsim,    mu = exp(eta)) 
if(margin == "ZTP")   y <- rZTP(  rsim,    mu = exp(eta)) 

if(setseed == TRUE) rm(list = ".Random.seed", envir = globalenv()) 

y

}    