r.resp <- function(margin, rsim, eta1, eta2, eta3){



if(margin == "ZTP") rZTP <- function(n, mu) qpois(runif(n, dpois(0, mu), 1), mu)

if(margin %in% c("DGP","DGPII")){

       rDGP <- function(n, eta1, eta2, margin){  

                minp <- distrHsATDiscr(0, eta1, esp.tr(eta2, margin)$vrb, NULL, margin, NULL, robust = TRUE)$pdf2 
                
                p    <- runif(rsim)

                indv <- p <= minp
                
                if(margin == "DGP")   q <- ifelse(indv == TRUE, 0, ceiling( esp.tr(eta2, margin)$vrb/eta1*( (1 - p)^(-eta1) - 1   ) - 1 ))
                if(margin == "DGPII") q <- ifelse(indv == TRUE, 0, ceiling( esp.tr(eta2, margin)$vrb/(eta1^2)*( (1 - p)^(-eta1^2) - 1   ) - 1 )) 
                
                q
                
                                      }
}


if(margin == "N")     y <- rNO(   rsim,    mu =     eta1,     sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "GU")    y <- rGU(   rsim,    mu =     eta1,     sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "rGU")   y <- rRG(   rsim,    mu =     eta1,     sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "LO")    y <- rLO(   rsim,    mu =     eta1,     sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "LN")    y <- rLOGNO(rsim,    mu =     eta1,     sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "WEI")   y <- rWEI(  rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "iG")    y <- rIG(   rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "GA")    y <- rGA(   rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "GAi")   y <- rGA(   rsim,    mu =     eta1,     sigma = esp.tr(eta2, margin)$vrb)


if(margin == "TW"){

    if(rsim == 1) y <- rTweedie(mu = exp(eta1), p = esp.tr(eta2, margin)$vrb, phi = exp(eta3))
    if(rsim > 1){ # this needs to be generalised

      y <- NA
      for(i in 1:rsim){
        y[i] <- rTweedie(mu = exp(eta1[i]), p = esp.tr(eta2[i], margin)$vrb, phi = exp(eta3[i]))
                      }
                   }

}



if(margin %in% c("GP","GPII","GPo")){

             if(rsim == 1){ 
             
               if(margin == "GP")    y <- rgpd(  rsim, loc = 0, shape = eta1,            scale = esp.tr(eta2, margin)$vrb)
               if(margin == "GPII")  y <- rgpd(  rsim, loc = 0, shape = exp(eta1) - 0.5, scale = esp.tr(eta2, margin)$vrb)
               if(margin == "GPo")   y <- rgpd(  rsim, loc = 0, shape = exp(eta1) - 0.5, scale = esp.tr(eta2, margin)$vrb/(1+(exp(eta1) - 0.5)))
               
             
             
             }
             
             if(rsim > 1){
                y <- NA
                if(length(eta2) == 1) sigma <- rep(exp(eta2), rsim) 
                for(i in 1:rsim){
                  if(margin == "GP")     y[i] <- rgpd(  1, loc = 0, shape = eta1[i],            scale = esp.tr(eta2[i], margin)$vrb)
                  if(margin == "GPII")   y[i] <- rgpd(  1, loc = 0, shape = exp(eta1[i]) - 0.5, scale = esp.tr(eta2[i], margin)$vrb)
                  if(margin == "GPo")    y[i] <- rgpd(  1, loc = 0, shape = exp(eta1[i]) - 0.5, scale = esp.tr(eta2[i], margin)$vrb/( 1 + (exp(eta1[i]) - 0.5) )  )
                }
             
             
             }
}


if(margin == "DAGUM") y <- rGB2(  rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb, nu = exp(eta3), tau = 1) 
if(margin == "SM")    y <- rGB2(  rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb, nu = 1 , tau = exp(eta3)) 
if(margin == "BE")    y <- rBE(   rsim,    mu = plogis(eta1), sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "FISK")  y <- rGB2(  rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb, nu = 1 , tau = 1)
if(margin == "NBI")   y <- rNBI(  rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "NBII")  y <- rNBII( rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "PIG")   y <- rPIG(  rsim,    mu = exp(eta1),    sigma = esp.tr(eta2, margin)$vrb) 
if(margin == "PO")    y <- rPO(   rsim,    mu = exp(eta1)) 
if(margin == "ZTP")   y <- rZTP(  rsim,    mu = exp(eta1)) 

if(margin %in% c("DGP", "DGPII"))   y <- rDGP(rsim, eta1, eta2, margin) 

if(margin %in% c("probit", "logit", "cloglog"))  y <- rbinom(rsim, 1, prob = probm(eta1, margin)$pr )
                                                

y

}    