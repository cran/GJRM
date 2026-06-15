res.check <- function(x, joint = FALSE, intervals = FALSE, n.sim = 1000, prob.lev = 0.05){

qr <- qr1 <- qr2 <- qr3 <- NULL
u1 <- u2 <- NA

if(!inherits(x, "gjrm.pcc")){

if(x$triv == TRUE) stop("Diagnostic plots not available for trivariate binary models.")
if(x$Model %in% c("BPO", "BPO0") ) stop("Diagnostic plots not available for partial observability models.")

}

if(joint == FALSE & inherits(x, "gjrm.pcc") ) stop("For PCC models, only joint = TRUE is supported (overall model assessment).")

if(joint == TRUE & inherits(x, "gjrm.pcc") ){


ros.resids <- x$ros.resids

qr1 <- qnorm(ros.resids[, 1])
qr2 <- qnorm(ros.resids[, 2])
qr3 <- qnorm(ros.resids[, 3])

qr <- qr1^2 + qr2^2 + qr3^2

 if(intervals == FALSE){
  qqplot(qchisq((1:length(qr) - 0.5)/length(qr), df = 3),
         sort(qr),
         xlab = "Theoretical Quantiles",
         ylab = "Sample Quantiles")
  abline(0, 1)
                       }

if(intervals == TRUE) .qqenv.chisq2(qr, n.sim, prob.lev, df = 3) 


}



if(joint == TRUE & !inherits(x, "gjrm.pcc") ){##**##


margins <- x$VC$margins 

if(margins[1] %in% x$VC$bl & margins[2] == "TW") stop("Joint residual diagnostic plots are not currently available for the fitted model. Get in touch to check progress.")
if(x$Model %in% c("BSS") ) stop("Joint residual diagnostic plots not available for sample selection models.")
if(x$Model == "ROY") stop("Joint residual diagnostic plots are not available for switching models.")
if(margins[1] %in% x$VC$bl & margins[2] %in% x$VC$bl & x$surv.flex == TRUE & x$end.surv == TRUE) stop("Joint residual diagnostic plots are not currently available for the fitted model. Get in touch to check progress.")


y1 <- x$y1
y2 <- x$y2
n  <- x$n

min.pr <- x$VC$min.pr
max.pr <- x$VC$max.pr
min.dn <- x$VC$min.dn

p1  <- x$fit$p1
p2  <- x$fit$p2
pf1 <- x$fit$pdf1
pf2 <- x$fit$pdf2

n1 <- sum(x$y1) 
n0 <- n - n1        


if(margins[1] %in% x$VC$bl & margins[2] %in% x$VC$bl & is.null(x$VC$K1) & is.null(x$VC$K2) & x$surv.flex == FALSE){ ####**** BINARY - BINARY ****#### 

u1[y1 == 0] <- runif(n0, 0              , 1 - p1[y1 == 0])
u1[y1 == 1] <- runif(n1, 1 - p1[y1 == 1], 1              )

## zero part

P00 <- x$fit$p00[y1 == 0]

u2[y1 == 0 & y2 == 0] <- runif(sum(y1 == 0 & y2 == 0), 0,                                               P00[y2[y1 == 0] == 0]/(1-p1[y1 == 0 & y2 == 0]))     
u2[y1 == 0 & y2 == 1] <- runif(sum(y1 == 0 & y2 == 1), P00[y2[y1 == 0] == 1]/(1-p1[y1 == 0 & y2 == 1]), 1                                              )
                            
## one part

P10 <- x$fit$p10[y1 == 1] 

u2[y1 == 1 & y2 == 0] <- runif(sum(y1 == 1 & y2 == 0), 0,                                           P10[y2[y1 == 1] == 0]/p1[y1 == 1 & y2 == 0])     
u2[y1 == 1 & y2 == 1] <- runif(sum(y1 == 1 & y2 == 1), P10[y2[y1 == 1] == 1]/p1[y1 == 1 & y2 == 1], 1                                          )
                              
							} ####**** END: BINARY - BINARY ****#### 





if(margins[1] %in% c(x$VC$m2, x$VC$m3) &&  margins[2] %in% c(x$VC$m2, x$VC$m3) & x$surv.flex == FALSE){ ####**** CONT - CONT ****#### & is.null(x$VC$K1) & is.null(x$VC$K2)

u1 <- p1 ### ** TW not dealt with at the moment ** ##
u2 <- x$fit$c.copula.be1   
										}  ####**** END: CONT - CONT ****#### 




if(margins[1] %in% x$VC$bl && margins[2] %in% c(x$VC$m1d, x$VC$m2d) & x$surv.flex == FALSE){ ####**** BINARY - DISCR ****#### & is.null(x$VC$K1) & is.null(x$VC$K2)


u1[y1 == 0] <- runif(n0, 0              , 1 - p1[y1 == 0])
u1[y1 == 1] <- runif(n1, 1 - p1[y1 == 1], 1              )


A1 <- x$fit$C1[y1 == 0]
A2 <- x$fit$C2[y1 == 0] 
             
A1.0 <- A1/(1-p1[y1 == 0])
A2.0 <- A2/(1-p1[y1 == 0])

u2[y1 == 0] <- runif(n0, A2.0, A1.0)
    
  
  
A1 <- mm( p2[y1 == 1] - x$fit$C1[y1 == 1], min.pr = min.pr, max.pr = max.pr )
A2 <- mm( (p2[y1 == 1] - pf2[y1 == 1]) - x$fit$C2[y1 == 1] , min.pr = min.pr, max.pr = max.pr )
                                                                                              
A1.1 <- A1/p1[y1 == 1]
A2.1 <- A2/p1[y1 == 1]

u2[y1 == 1] <- runif(n1, A2.1, A1.1)



} ####**** END: BINARY - DISCR ****#### 









if(margins[1] %in% c(x$VC$m1d, x$VC$m2d) &&  margins[2] %in% c(x$VC$m1d, x$VC$m2d) & x$surv.flex == FALSE){ ####**** DISCR - DISCR ****#### & is.null(x$VC$K1) & is.null(x$VC$K2)

u1 <- runif(n, p1 - pf1, p1)

upB <- (x$fit$C11 - x$fit$C01) / pf1
loB <- (x$fit$C10 - x$fit$C00) / pf1

tmp <- pmin(loB, upB)
upB <- pmax(loB, upB)
loB <- tmp

loB <- mm(loB, min.pr = min.pr, max.pr = max.pr)
upB <- mm(upB, min.pr = min.pr, max.pr = max.pr)

u2 <- runif(n, loB, upB)
                                                                                     } ####**** END: DISC - DISCR ****#### 







if(margins[1] %in% x$VC$bl & margins[2] %in% x$VC$bl & !is.null(x$VC$K1) & !is.null(x$VC$K2) & x$surv.flex == FALSE){ ####**** ORDINAL - ORDINAL ****#### 


u1 <- runif(n, x$fit$p1.m, x$fit$p1.p)

pf1 <- x$fit$p1.p - x$fit$p1.m

upB <- (x$fit$C11 - x$fit$C01) / pf1
loB <- (x$fit$C10 - x$fit$C00) / pf1

tmp <- pmin(loB, upB)
upB <- pmax(loB, upB)
loB <- tmp

loB <- mm(loB, min.pr = min.pr, max.pr = max.pr)
upB <- mm(upB, min.pr = min.pr, max.pr = max.pr)

u2 <- runif(n, loB, upB)

                                                                                               } ####**** ORDINAL - ORDINAL ****#### 




if(margins[1] %in% x$VC$bl & margins[2] %in% c(x$VC$m2, x$VC$m3) & !is.null(x$VC$K1) & is.null(x$VC$K2) & x$surv.flex == FALSE){ ####**** ORDINAL - CONT ****#### 

u1  <- runif(n, x$fit$p1.m, x$fit$p1.p)
u2  <- mm( (x$fit$C11 - x$fit$C01) / (x$fit$p1.p - x$fit$p1.m), min.pr = min.pr, max.pr = max.pr)

                                                                                               } ####**** ORDINAL - CONT ****#### 




if(margins[1] %in% x$VC$bl &&  margins[2] %in% c(x$VC$m2, x$VC$m3) & x$surv.flex == FALSE & is.null(x$VC$K1) & is.null(x$VC$K2)){ ####**** BINARY - CONT ****#### 


u1[y1 == 0] <- runif(n0, 0              , 1 - p1[y1 == 0])
u1[y1 == 1] <- runif(n1, 1 - p1[y1 == 1], 1              )


u2[y1 == 0] <- x$fit$C01[y1 == 0]/(1 - p1[y1 == 0])
u2[y1 == 1] <- (p2[y1 == 1] - x$fit$C01[y1 == 1])/p1[y1 == 1]

u2 <- mm(u2, min.pr = min.pr, max.pr = max.pr)


} ####**** END: BINARY - CONT ****#### 






if(margins[1] %in% c(x$VC$m1d, x$VC$m2d) & margins[2] %in% c(x$VC$m2, x$VC$m3) & x$surv.flex == FALSE){ ####**** DISCR - CONT ****#### & is.null(x$VC$K1) & is.null(x$VC$K2)


u1 <- runif(n, p1 - pf1, p1)

u2 <- (x$fit$C11 - x$fit$C01)/pf1

u2 <- mm(u2, min.pr = min.pr, max.pr = max.pr)



}




if(margins[1] %in% x$VC$bl & margins[2] %in% x$VC$bl & x$surv.flex == TRUE ){ ####**** SURV - SURV ****#### 

F1U <- F2U <- NULL

onlyUR <- is.numeric(x$cens1) && is.numeric(x$cens2)

                        F1L <- 1 - x$fit$p1
if(any(x$cens1 == "I")) F1U <- 1 - x$fit$p1.2

                        F2L <- 1 - x$fit$p2
if(any(x$cens2 == "I")) F2U <- 1 - x$fit$p2.2
  
  
u1 <- pnorm( rqres_nrsp_F(F1 = F1L, cens = x$cens1, F2 = F1U) )



if(onlyUR){

indUU <- x$cens1 == 1 & x$cens2 == 1
indUR <- x$cens1 == 1 & x$cens2 == 0
indRU <- x$cens1 == 0 & x$cens2 == 1
indRR <- x$cens1 == 0 & x$cens2 == 0


if(any(indUU)) u2[indUU] <- 1 - x$fit$c.copula.be1[indUU]

if(any(indUR)) u2[indUR] <- runif(sum(indUR), 1 - x$fit$c.copula.be1[indUR], 1)

if(any(indRU)){

        C12L <- F1L[indRU] + F2L[indRU] - 1 + x$fit$p00[indRU]
        u2[indRU] <- mm((F2L[indRU] - C12L) / x$fit$p1[indRU], min.pr, max.pr)
        
               }

if(any(indRR)){

        C12L <- F1L[indRR] + F2L[indRR] - 1 + x$fit$p00[indRR]
        U2L <- mm((F2L[indRR] - C12L) / x$fit$p1[indRR], min.pr, max.pr)
        u2[indRR] <- runif(sum(indRR), U2L, 1)
        
               }

}







if(!onlyUR){


indII <- x$cens1 == "I" & x$cens2 == "I"
indIL <- x$cens1 == "I" & x$cens2 == "L"
indIR <- x$cens1 == "I" & x$cens2 == "R"
indIU <- x$cens1 == "I" & x$cens2 == "U"

indLI <- x$cens1 == "L" & x$cens2 == "I"
indLL <- x$cens1 == "L" & x$cens2 == "L"
indLR <- x$cens1 == "L" & x$cens2 == "R"
indLU <- x$cens1 == "L" & x$cens2 == "U"

indRI <- x$cens1 == "R" & x$cens2 == "I"
indRL <- x$cens1 == "R" & x$cens2 == "L"
indRR <- x$cens1 == "R" & x$cens2 == "R"
indRU <- x$cens1 == "R" & x$cens2 == "U"

indUI <- x$cens1 == "U" & x$cens2 == "I"
indUL <- x$cens1 == "U" & x$cens2 == "L"
indUR <- x$cens1 == "U" & x$cens2 == "R"
indUU <- x$cens1 == "U" & x$cens2 == "U"


if(any(indII)){

  den <- x$fit$p1[indII] - x$fit$p1.2[indII]

  U2L <- mm( 1 - (x$fit$p00[indII]      - x$fit$p00.mix2[indII]) / den , min.pr, max.pr)
  U2U <- mm( 1 - (x$fit$p00.mix1[indII] - x$fit$p00.2[indII])    / den , min.pr, max.pr)

  u2[indII] <- runif(sum(indII), U2L, U2U)

               }


if(any(indIL)){

  den <- x$fit$p1[indIL] - x$fit$p1.2[indIL]

  U2U <- mm( 1 - (x$fit$p00[indIL] - x$fit$p00.mix2[indIL]) / den, min.pr, max.pr)

  u2[indIL] <- runif(sum(indIL), 0, U2U)

}

if(any(indIR)){

  den <- x$fit$p1[indIR] - x$fit$p1.2[indIR]

  U2L <- mm(1 - (x$fit$p00[indIR] - x$fit$p00.mix2[indIR]) / den, min.pr, max.pr)

  u2[indIR] <- runif(sum(indIR), U2L, 1)

}

if(any(indIU)){

  den <- x$fit$p1[indIU] - x$fit$p1.2[indIU]

  u2[indIU] <- mm(1 - (x$fit$p00[indIU] - x$fit$p00.mix2[indIU]) / den, min.pr, max.pr)


}




if(any(indLI)){

  	C12L <- F1L[indLI] + F2L[indLI] - 1 + x$fit$p00[indLI]
  	C12U <- F1L[indLI] + F2U[indLI] - 1 + x$fit$p00.mix1[indLI]
  	U2L <- mm(C12L / F1L[indLI], min.pr, max.pr)
  	U2U <- mm(C12U / F1L[indLI], min.pr, max.pr)
  	u2[indLI] <- runif(sum(indLI), U2L, U2U)

	      }

if(any(indLL)){

  	C12L <- F1L[indLL] + F2L[indLL] - 1 + x$fit$p00[indLL]
  	U2U  <- mm(C12L / F1L[indLL], min.pr, max.pr)
  	u2[indLL] <- runif(sum(indLL), 0, U2U)

	      }

if(any(indLR)){

  	C12L <- F1L[indLR] + F2L[indLR] - 1 + x$fit$p00[indLR]
  	U2L <- mm(C12L / F1L[indLR], min.pr, max.pr)
  	u2[indLR] <- runif(sum(indLR), U2L, 1)

              }


if(any(indLU)){

  	C12L <- F1L[indLU] + F2L[indLU] - 1 + x$fit$p00[indLU]
  	U2 <- mm(C12L / F1L[indLU], min.pr, max.pr)
  	u2[indLU] <- U2

              }

if(any(indRI)){  

        C12L <- F1L[indRI] + F2L[indRI] - 1 + x$fit$p00[indRI]
        C12U <- F1L[indRI] + F2U[indRI] - 1 + x$fit$p00.mix1[indRI]
        U2L <- mm( (F2L[indRI] - C12L)/x$fit$p1[indRI], min.pr, max.pr)
        U2U <- mm( (F2U[indRI] - C12U)/x$fit$p1[indRI], min.pr, max.pr) 
        u2[indRI] <- runif(sum(indRI), U2L, U2U) 

              }
              
if(any(indRL)){  

        C12L <- F1L[indRL] + F2L[indRL] - 1 + x$fit$p00[indRL]
        U2L <- mm( (F2L[indRL] - C12L)/x$fit$p1[indRL], min.pr, max.pr)
        u2[indRL] <- runif(sum(indRL), 0, U2L) # should be called U2U but it is the same

              }            
 
if(any(indRR)){

        C12L <- F1L[indRR] + F2L[indRR] - 1 + x$fit$p00[indRR]
        U2L <- mm((F2L[indRR] - C12L) / x$fit$p1[indRR], min.pr, max.pr)
        u2[indRR] <- runif(sum(indRR), U2L, 1)
        
               } 
               
if(any(indRU)){

        C12L <- F1L[indRU] + F2L[indRU] - 1 + x$fit$p00[indRU]
        u2[indRU] <- mm((F2L[indRU] - C12L) / x$fit$p1[indRU], min.pr, max.pr)
        
               }               
 

if(any(indUI)) u2[indUI] <- runif(sum(indUI), 1 - x$fit$c.copula.be1[indUI], 1 - x$fit$c.copula.be1.mix1[indUI])   
if(any(indUL)) u2[indUL] <- runif(sum(indUL), 0, 1 - x$fit$c.copula.be1[indUL])
if(any(indUR)) u2[indUR] <- runif(sum(indUR), 1 - x$fit$c.copula.be1[indUR], 1)
if(any(indUU)) u2[indUU] <- 1 - x$fit$c.copula.be1[indUU]

} # end IF cond


}






qr1 <- qnorm(u1)
qr2 <- qnorm(u2)
qr  <- qr1^2 + qr2^2

lims <- max(abs(c(qr1, qr2)))
lims <- c(-lims, lims)

biv_std_norm <- function(x, y) { (1 / (2 * pi)) * exp(-0.5 * (x^2 + y^2)) }

# Grid for contours
grid_len <- 100
x_seq <- seq(lims[1], lims[2], length.out = grid_len)
y_seq <- seq(lims[1], lims[2], length.out = grid_len)
z <- outer(x_seq, y_seq, biv_std_norm)


radii <- seq(0.5, 3, length.out = 5)
levels_gauss <- (1 / (2 * pi)) * exp(-0.5 * radii^2)

# par(mar = c(5, 5, 2, 2)) #   oma = c(0, 0, 0, 0)) 

par(mfrow = c(1, 2), mar = c(5, 5, 2, 2))

plot(qr1, qr2, xlab = expression(e[1]), ylab = expression(e[2*"|"*1]), xlim = lims, ylim = lims, asp = 1, col = "grey50")
contour(x_seq, y_seq, z, levels = levels_gauss, add = TRUE, drawlabels = FALSE)

if(intervals == FALSE){ qqplot(qchisq((1:length(qr) - 0.5) / length(qr), df = 2), sort(qr), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1) }
if(intervals == TRUE) .qqenv.chisq2(qr, n.sim, prob.lev, df = 2)

print(cor.test(qr1, qr2, method = "kendall"))





} ##**## joint = TRUE







if(joint == FALSE){##**##


if(inherits(x, "gjrm.pcc")) stop("For PCC models, only joint = TRUE is supported (overall model assessment).")


y1m <- y2m <- y3m <- NA # no longer used...but left for the time being



if(x$univar.gamlss == FALSE){

if(x$surv.flex == TRUE){

if(x$end.surv == FALSE) par(mfrow = c(1, 2))

if(x$margins[1] %in% x$bl && x$end.surv == FALSE){


qr1 <- rqres_nrsp_F(F1 = 1 - x$fit$p1, cens = x$cens1, F2 = 1 - x$fit$p1.2)

if(intervals == FALSE){qqplot(qnorm((1:length(qr1) - 0.5) / length(qr1)), sort(qr1), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}         
if(intervals == TRUE){.qqenv(qr1, n.sim, prob.lev)}


}



if(x$margins[1] %in% c(x$VC$m2,x$VC$m3)){

p1  <- x$fit$p1

qr1 <- qnorm(p1)

if(intervals == FALSE){qqplot(qnorm((1:length(qr1) - 0.5) / length(qr1)), sort(qr1), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}         
if(intervals == TRUE){.qqenv(qr1, n.sim, prob.lev)}

}


qr2 <- rqres_nrsp_F(F1 = 1 - x$fit$p2, cens = x$cens2, F2 = 1 - x$fit$p2.2)


if(intervals == FALSE){qqplot(qnorm((1:length(qr2) - 0.5) / length(qr2)), sort(qr2), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}         
if(intervals == TRUE){.qqenv(qr2, n.sim, prob.lev)}



}






if(x$surv.flex == FALSE & !is.null(x$VC$K1)){

par(mfrow = c(1, 2))

is_ordcon <- !is.null(x$VC$K1) &  is.null(x$VC$K2)
is_ordord <- !is.null(x$VC$K1) & !is.null(x$VC$K2)



if(is_ordord == TRUE){

qr1 <- rqres_ordinal(x$fit$p1.m, x$fit$p1.p)   
qr2 <- rqres_ordinal(x$fit$p2.m, x$fit$p2.p)   

if(intervals == FALSE){qqplot(qnorm((1:length(qr1) - 0.5) / length(qr1)), sort(qr1), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)
                       qqplot(qnorm((1:length(qr2) - 0.5) / length(qr2)), sort(qr2), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)
                      }         
if(intervals == TRUE){.qqenv(qr1, n.sim, prob.lev); .qqenv(qr2, n.sim, prob.lev)}

                } # is_ordord


if(is_ordcon == TRUE){

qr1 <- rqres_ordinal(x$fit$p1.m, x$fit$p1.p)   
p2  <- x$fit$p2
qr2 <- qnorm(p2)

if(intervals == FALSE){qqplot(qnorm((1:length(qr1) - 0.5) / length(qr1)), sort(qr1), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)
                       qqplot(qnorm((1:length(qr2) - 0.5) / length(qr2)), sort(qr2), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)
                      }         
if(intervals == TRUE){.qqenv(qr1, n.sim, prob.lev); .qqenv(qr2, n.sim, prob.lev)}

                } # is_ordcon




}





if(x$surv.flex == FALSE & is.null(x$VC$K1)){##


if(x$VC$margins[1] %in% x$bl && x$VC$margins[2] %in% x$bl ) stop("Diagnostic plots not available for binary margins.")






if(x$Cont == "NO"){


y2     <- x$y2
eta2   <- x$eta2
sigma2 <- x$sigma2
nu     <- x$nu


if(x$VC$ccss == "yes"){ # for the time being, to be coherent with eta1, I use eta2 etc. calculated from X2s etc. (more relevant for splines)
                        # same reasoning for Roy models
   
   p1  <- x$p1[x$inde]                                            
   p0  <- mm( 1 - p1, min.pr = 1e-300, max.pr = 0.9999999999999999)     
   
   eta2   <- eta2[x$inde]
   theta  <- x$theta
   
   if(length(theta) > 1) theta <- theta[x$inde]

   if(!is.null(x$X3) && !(x$VC$margins[2] %in% x$VC$m1d) ){
   sigma2 <- sigma2[x$inde]
   nu     <- nu[x$inde]
                                                          }
                                                          
                      }
                      
                      




if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3) && x$VC$ccss != "yes" && x$Model != "ROY"){

                                            p <- x$fit$p2
                                              if(x$VC$margins[2] %in% c("TW")){ if( any(y2 == 0) == TRUE ) p[y2 == 0] <- runif(sum(y2 == 0), min = 1e-300, max = p[y2 == 0]) } 
                                            qr2 <- qnorm(p)
 
                                            } 
                                            
                                            
if(x$VC$margins[2] %in% c(x$VC$m2, x$VC$m3) && x$VC$ccss == "yes" && x$Model != "ROY"){

                                            p2  <- distrHsAT(y2, eta2, sigma2, nu, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0, p2, x$nC, theta, x$dof), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm( (p2 - p12)/p1, min.pr = 1e-300, max.pr = 0.9999999999999999)                                                  
                                              if(x$VC$margins[2] %in% c("TW")){ if( any(y2 == 0) == TRUE ) p[y2 == 0] <- runif(sum(y2 == 0), min = 1e-300, max = p[y2 == 0]) }
                                            qr2  <- qnorm(p)
                                                                                        
                                            } 
                                            
                                            
if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d) && x$VC$ccss != "yes" && x$Model != "ROY"){ 

                                            p <- x$fit$p2 
                                            d <- x$fit$pdf2 
                                            qr2 <- qnorm( runif(length(y2), p - d, p))                                           
                                            } 
                                            
                                            
if(x$VC$margins[2] %in% c(x$VC$m1d, x$VC$m2d, x$VC$m3d) && x$VC$ccss == "yes" && x$Model != "ROY"){ 

                                            p2  <- distrHsATDiscr(y2, eta2, sigma2, nu, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0, p2, x$nC, theta, x$dof), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm( (p2 - p12)/p1, min.pr = 1e-300, max.pr = 0.9999999999999999)  
                                            
                                            y2a <- y2; y2a <- y2 - 1; ind.y2a <- y2a < 0; y2a <- ifelse(y2a < 0, 0, y2a)
                                            p2  <- distrHsATDiscr(y2a, eta2, sigma2, nu, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2 
                                            p12 <- mm(BiCDF(p0, p2, x$nC, theta, x$dof), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            pm1 <- mm( (p2 - p12)/p1, min.pr = 1e-300, max.pr = 0.9999999999999999)   
                                            pm1[ind.y2a] <- 0  
                                            
                                            pm1 <- ifelse(pm1 < 0, 0, pm1) # these needed for very bad model fits
                                            p   <- ifelse(p   < 0, 0, p  )
                                            pm1 <- ifelse(pm1 > p, p, pm1)
                                            ru  <- runif(length(y2), pm1, p)
                                            ru  <- ifelse(ru == 0, 1e-16, ru)      
                                            ru  <- ifelse(ru > 0.999999, 0.999999, ru)                                     
                                            qr2 <- qnorm( ru)  
                                            
                                            }                                             
                                            
                                             
   
if(x$Model == "ROY"){                                           

   y2 <- x$y2
   y3 <- x$y3
   
   eta2 <- x$eta2[x$inde0]
   eta3 <- x$eta3[x$inde1]
   
   sigma1 <- x$sigma2[x$inde0]
   sigma2 <- x$sigma3[x$inde1]

   nu1 <- x$nu2[x$inde0]
   nu2 <- x$nu3[x$inde1] 
   
   theta1 <- x$theta12[x$inde0]
   theta2 <- x$theta13[x$inde1]
   
   p1eq2 <- x$p1[x$inde0]  
   p0eq2 <- 1 - p1eq2
   p1eq3 <- x$p1[x$inde1]   
   p0eq3 <- 1 - p1eq3 
   

if(x$VC$margins[2] %in% c(x$VC$m2, x$VC$m3)){

                                            p2  <- distrHsAT(y2, eta2, sigma1, nu1, x$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2
                                            p12 <- mm(BiCDF(p0eq2, p2, x$nC1, theta1, x$dof12), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)
                                            p   <- mm(p12/p0eq2, min.pr = 1e-300, max.pr = 0.9999999999999999)                                              
                                              if(x$VC$margins[2] %in% c("TW")){ if( any(y2 == 0) == TRUE ) p[y2 == 0] <- runif(sum(y2 == 0), min = 1e-300, max = p[y2 == 0]) }
                                            qr2  <- qnorm(p)
                                                                                        
                                            }

if(x$VC$margins[3] %in% c(x$VC$m2, x$VC$m3)){

                                            p2  <- distrHsAT(y3, eta3, sigma2, nu2, x$margins[3], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0eq3, p2, x$nC2, theta2, x$dof13), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm((p2 - p12)/p1eq3, min.pr = 1e-300, max.pr = 0.9999999999999999)                                             
                                              if(x$VC$margins[3] %in% c("TW")){ if( any(y3 == 0) == TRUE ) p[y3 == 0] <- runif(sum(y3 == 0), min = 1e-300, max = p[y3 == 0]) }
                                            qr3 <- qnorm(p)
                                                                                        
                                            }                                             
                                            
if(x$VC$margins[2] %in% c(x$VC$m1d, x$VC$m2d, x$VC$m3d)){ 

                                            p2  <- distrHsATDiscr(y2, eta2, sigma1, nu1, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2
                                            p12 <- mm(BiCDF(p0eq2, p2, x$nC1, theta1, x$dof12), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm(p12/p0eq2, min.pr = 1e-300, max.pr = 0.9999999999999999) 
                                            
                                            y2a <- y2; y2a <- y2 - 1; ind.y2a <- y2a < 0; y2a <- ifelse(y2a < 0, 0, y2a)
                                            p2  <- distrHsATDiscr(y2a, eta2, sigma1, nu1, x$margins[2], y2m = y2m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc1)$p2 
                                            p12 <- mm(BiCDF(p0eq2, p2, x$nC1, theta1, x$dof12), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            pm1 <- mm(p12/p0eq2, min.pr = 1e-300, max.pr = 0.9999999999999999)  
                                            pm1[ind.y2a] <- 0  
                                            
                                            pm1 <- ifelse(pm1 < 0, 0, pm1) 
                                            p   <- ifelse(p   < 0, 0, p  )
                                            pm1 <- ifelse(pm1 > p, p, pm1)
                                            ru  <- runif(length(y2), pm1, p)
                                            ru  <- ifelse(ru == 0, 1e-16, ru)  
                                            ru  <- ifelse(ru > 0.999999, 0.999999, ru)                                                                                
                                            qr2 <- qnorm(ru)                                              
                                            
                                            }

if(x$VC$margins[3] %in% c(x$VC$m1d, x$VC$m2d, x$VC$m3d)){ 

                                            p2  <- distrHsATDiscr(y3, eta3, sigma2, nu2, x$margins[3], y2m = y3m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2
                                            p12 <- mm(BiCDF(p0eq3, p2, x$nC2, theta2, x$dof13), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            p   <- mm((p2 - p12)/p1eq3, min.pr = 1e-300, max.pr = 0.9999999999999999) 
                                            
                                            y3a <- y3; y3a <- y3 - 1; ind.y3a <- y3a < 0; y3a <- ifelse(y3a < 0, 0, y3a)
                                            p2  <- distrHsATDiscr(y3a, eta3, sigma2, nu2, x$margins[3], y2m = y3m, min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, 
                                                                  max.pr = x$VC$max.pr, left.trunc = x$VC$left.trunc2)$p2 
                                            p12 <- mm(BiCDF(p0eq3, p2, x$nC2, theta2, x$dof13), min.pr = x$VC$min.pr, max.pr = x$VC$max.pr )
                                            pm1 <- mm((p2 - p12)/p1eq3, min.pr = 1e-300, max.pr = 0.9999999999999999)  
                                            pm1[ind.y3a] <- 0  
                                            
                                            pm1 <- ifelse(pm1 < 0, 0, pm1) 
                                            p   <- ifelse(p   < 0, 0, p  )  
                                            pm1 <- ifelse(pm1 > p, p, pm1)
                                            ru  <- runif(length(y3), pm1, p)
                                            ru  <- ifelse(ru == 0, 1e-16, ru)   
                                            ru  <- ifelse(ru > 0.999999, 0.999999, ru)                                            
                                            qr3 <- qnorm(ru)  
                                            
                                            }  


                                          
}                                            
                                            
                                            
                                                    


if(x$Model != "ROY"){

if(intervals == FALSE){qqplot(qnorm((1:length(qr2) - 0.5) / length(qr2)), sort(qr2), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}
if(intervals == TRUE){.qqenv(qr2, n.sim, prob.lev)}

}


if(x$Model == "ROY"){

par(mfrow = c(1,2))

if(intervals == FALSE){ qqplot(qnorm((1:length(qr2) - 0.5) / length(qr2)), sort(qr2), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1); 
                        qqplot(qnorm((1:length(qr3) - 0.5) / length(qr3)), sort(qr3), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}
if(intervals == TRUE){.qqenv(qr2, n.sim, prob.lev); .qqenv(qr3, n.sim, prob.lev)}

}

} # end Cont == NO




if(x$Cont == "YES"){

y1 <- x$y1
y2 <- x$y2

y1m <- y2m <- NA




if(x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3)) {p1 <- x$fit$p1  

                                            if(x$VC$margins[1] %in% c("TW")){ if( any(x$y1 == 0) == TRUE ) p1[x$y1 == 0] <- runif(sum(x$y1 == 0), min = 1e-300, max = p1[x$y1 == 0]) }


}


if(x$VC$margins[1] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){

p <- x$fit$p1
d <- x$fit$pdf1   
p1 <- runif(length(y1), p - d, p) 

}



if(x$VC$margins[2] %in% c(x$VC$m2,x$VC$m3)){   p2 <- x$fit$p2  

                                            if(x$VC$margins[2] %in% c("TW")){ if( any(x$y2 == 0) == TRUE ) p2[x$y2 == 0] <- runif(sum(x$y2 == 0), min = 1e-300, max = p2[x$y2 == 0])   }

}

if(x$VC$margins[2] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){


p <- x$fit$p2 
d <- x$fit$pdf2   

p2 <- runif(length(y2), p - d, p) 

}


par(mfrow = c(1,2))

qr1 <- qnorm(p1)
if(intervals == FALSE){qqplot(qnorm((1:length(qr1) - 0.5) / length(qr1)), sort(qr1), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}
if(intervals == TRUE){.qqenv(qr1, n.sim, prob.lev)}

qr2 <- qnorm(p2)
if(intervals == FALSE){qqplot(qnorm((1:length(qr2) - 0.5) / length(qr2)), sort(qr2), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}
if(intervals == TRUE){.qqenv(qr2, n.sim, prob.lev)}

}


}##

}###




if(x$univar.gamlss == TRUE){

if(x$surv.flex == FALSE){ ###

if(x$VC$margins[1] %in% c("GEVlink") ) stop("Diagnostic plot not available since residual checks are uninformative.")


y1 <- x$y1

                         


if(x$VC$margins[1] %in% c(x$VC$m2,x$VC$m3)) p1 <- x$fit$p 

        if(x$VC$margins[1] %in% c("TW")){ if( any(x$y1 == 0) == TRUE )  p1[x$y1 == 0] <- runif(sum(x$y1 == 0), min = 1e-300, max = p1[x$y1 == 0]) }

if(x$VC$margins[1] %in% c(x$VC$m1d,x$VC$m2d,x$VC$m3d)){

      p <- x$fit$p
      d <- x$fit$pdf2   
      p1 <- runif(length(y1), p - d, p) 
                                                      }

qr <- qnorm(p1)

if(intervals == FALSE){qqplot(qnorm((1:length(qr) - 0.5) / length(qr)), sort(qr), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}
if(intervals == TRUE){.qqenv(qr, n.sim, prob.lev)}

}###


if(x$surv.flex == TRUE){ ###

qr <- rqres_nrsp_F(F1 = 1 - x$fit$p1, cens = x$cens, F2 = 1 - x$fit$p2)


if(intervals == FALSE){qqplot(qnorm((1:length(qr) - 0.5) / length(qr)), sort(qr), xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); abline(0, 1)}
if(intervals == TRUE){.qqenv(qr, n.sim, prob.lev)}

}###


}




} ##**## joint == FALSE


L <- list(qr = as.numeric(qr), qr1 = as.numeric(qr1), qr2 = as.numeric(qr2), qr3 = as.numeric(qr3))

invisible(L)


}

