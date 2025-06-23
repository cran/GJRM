sim.resp <- function(margin, rsim, eta, sigma2, nu, setseed = TRUE, left.trunc = 0){



if(margin %in% c("DGP","DGPII","DGP0")){

       rDGP <- function(n, mu, sigma, margin){ # mu is eta here

                minp <- distrHsATDiscr(0, mu, sigma, NULL, margin, NULL, robust = TRUE, min.dn = 1e-323, min.pr = 1e-16, max.pr = 0.999999)$pdf2
                
                p    <- runif(rsim)

                indv <- p <= minp
                

                if(margin == "DGP")   q <- ifelse(indv == TRUE, 0, ceiling( sigma/mu*(     (1 - p)^(-mu)   - 1   )) - 1 )
                if(margin == "DGPII") q <- ifelse(indv == TRUE, 0, ceiling( sigma/(exp(mu))*( (1 - p)^(-exp(mu)) - 1   )) - 1 )
                
                if(margin == "DGP0")  q <- ifelse(indv == TRUE, 0, ceiling(-exp(mu)*log(1-p)) - 1 )                 
                
                
                
                
                q
                
                                      }
}


if(setseed == TRUE) set.seed(1)

if(margin == "N")     y <- rNO(   rsim,    mu =     eta,     sigma = sigma2) 
if(margin == "tN")    y <- rNtr(  rsim,    mu =     eta,     sigma = sigma2, left.trunc = left.trunc) 
if(margin == "GU")    y <- rGU(   rsim,    mu =     eta,     sigma = sigma2) 
if(margin == "rGU")   y <- rRG(   rsim,    mu =     eta,     sigma = sigma2) 
if(margin == "LO")    y <- rLO(   rsim,    mu =     eta,     sigma = sigma2) 
if(margin == "LN")    y <- rLOGNO(rsim,    mu =     eta,     sigma = sigma2) 
if(margin == "WEI")   y <- rWEI(  rsim,    mu = exp(eta),    sigma = sigma2) 
if(margin == "IG")    y <- rIG(   rsim,    mu = exp(eta),    sigma = sigma2) 
if(margin == "GA")    y <- rGA(   rsim,    mu = exp(eta),    sigma = sigma2) 
if(margin == "GAi")   y <- rGA(   rsim,    mu =     eta,     sigma = sigma2)






if(margin == "TW"){    


     if(rsim == 1) y <- rTweedie(mu = exp(eta), phi = sigma2, p = nu)
     
     if(rsim > 1){
           y <- NA
           if(length(sigma2) == 1) sigma2 <- rep(sigma2, rsim) 
           if(length(nu)     == 1) nu     <- rep(nu, rsim) 
           for(i in 1:rsim) y[i] <- rTweedie(mu = exp(eta[i]), phi = sigma2[i], p = nu[i])                
                 }     
     
}


if(margin %in% c("GP","GPII","GPo")){

             if(rsim == 1){ 
             
               if(margin == "GP")    y <- rgpd(  rsim, loc = 0, shape = eta,            scale = sigma2)
               if(margin == "GPII")  y <- rgpd(  rsim, loc = 0, shape = exp(eta) - 0.5, scale = sigma2)
               if(margin == "GPo")   y <- rgpd(  rsim, loc = 0, shape = exp(eta) - 0.5, scale = sigma2/(1+(exp(eta) - 0.5)))
               
             
             
             }
             
             if(rsim > 1){
                y <- NA
                if(length(sigma2) == 1) sigma2 <- rep(sigma2, rsim) 
                for(i in 1:rsim){
                  if(margin == "GP")     y[i] <- rgpd(  1, loc = 0, shape = eta[i],            scale = sigma2[i])
                  if(margin == "GPII")   y[i] <- rgpd(  1, loc = 0, shape = exp(eta[i]) - 0.5, scale = sigma2[i])
                  if(margin == "GPo")    y[i] <- rgpd(  1, loc = 0, shape = exp(eta[i]) - 0.5, scale = sigma2[i]/( 1 + (exp(eta[i]) - 0.5) )  )
                }
             
             
             }
}


if(margin == "DAGUM") y <- rGB2(  rsim,    mu = exp(eta),    sigma = sigma2, nu = nu, tau = 1) 
if(margin == "SM")    y <- rGB2(  rsim,    mu = exp(eta),    sigma = sigma2, nu = 1 , tau = nu) 
if(margin == "BE")    y <- rBE(   rsim,    mu = plogis(eta), sigma = sigma2) 
if(margin == "FISK")  y <- rGB2(  rsim,    mu = exp(eta),    sigma = sigma2, nu = 1 , tau = 1)
if(margin == "NBI")   y <- rNBI(  rsim,    mu = exp(eta),    sigma = sigma2) 
if(margin == "NBII")  y <- rNBII( rsim,    mu = exp(eta),    sigma = sigma2) 
if(margin == "PIG")   y <- rPIG(  rsim,    mu = exp(eta),    sigma = sigma2) 
if(margin == "P")     y <- rPO(   rsim,    mu = exp(eta)) 
if(margin == "tP")     y <- rPtr(  rsim,    mu = exp(eta), left.trunc = left.trunc)
if(margin == "tNBI")   y <- rNBItr(  rsim,    mu = exp(eta),    sigma = sigma2, left.trunc = left.trunc) 
if(margin == "tNBII")  y <- rNBIItr( rsim,    mu = exp(eta),    sigma = sigma2, left.trunc = left.trunc) 
if(margin == "tPIG")   y <- rPIGtr(  rsim,    mu = exp(eta),    sigma = sigma2, left.trunc = left.trunc) 

if(margin %in% c("DGP", "DGPII", "DGP0")){

    if(margin %in% c("DGP0")) sigma2 <- NULL

y <- rDGP(rsim, mu = eta, sigma = sigma2, margin = margin) # 


}

if(margin %in% c("probit", "logit", "cloglog"))  y <- rbinom(rsim, 1, prob = probm(eta, margin, min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999)$pr )
                                                

if(setseed == TRUE) rm(list = ".Random.seed", envir = globalenv()) 

y

}    