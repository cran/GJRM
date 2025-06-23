jc.probs1 <- function(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, cumul = "no", tau.res = FALSE){


######################################################################################################

nu1 <- nu2 <- nu <- sigma2 <- 1
CIp12 <- dof <- p12s <- dofs <- C1s <- C2s <- C11s <- C01s <- C10s <- C00s <- p1 <- NULL
CIkt <- tau <- theta <- CItheta <- NULL
p12s  <- matrix(0, 1, 2) 

######################################################################################################


if(type == "joint"){ 


######################################################################
######################################################################
# set up
######################################################################
######################################################################

if(!missing(newdata)){ #################

nu1 <- nu2 <- dof <- sigma21 <- sigma22 <- NA

eta1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata)
eta2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata)

if( !is.null(x$X3) ){ # X3


if(x$margins[1] %in% c(cont2par, cont3par))                                           sigma21 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata), x$margins[1])$vrb
if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par, cont3par)) sigma22 <- esp.tr(predict.SemiParBIV(x, eq = 4, newdata = newdata), x$margins[2])$vrb
if(x$margins[1] %in% cont1par && x$margins[2] %in% c(cont2par,cont3par))              sigma22 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata), x$margins[2])$vrb


if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) ){

if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6;  eq.dof <- 7; eq.th <- 8}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.dof <- 5; eq.th <- 6}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5;  eq.dof <- 6; eq.th <- 7}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA; eq.dof <- 6; eq.th <- 7}

dof <- dof.tr(predict.SemiParBIV(x, eq = eq.dof, newdata = newdata))$vao 


                   }else{
                   
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6;  eq.th <- 7}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 5}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5;  eq.th <- 6}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA; eq.th <- 6} 

if(x$margins[1] %in% cont1par && x$margins[2] %in% cont1par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 3}
if(x$margins[1] %in% cont1par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA; eq.th <- 4}
if(x$margins[1] %in% cont1par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 4 ; eq.th <- 5}
                   
}                   
                   


if(x$margins[1] %in% cont3par) nu1 <- enu.tr(predict.SemiParBIV(x, eq = eq.nu1, newdata = newdata), x$margins[1])$vrb
if(x$margins[2] %in% cont3par) nu2 <- enu.tr(predict.SemiParBIV(x, eq = eq.nu2, newdata = newdata), x$margins[2])$vrb

theta <- teta.tr(x$VC, predict.SemiParBIV(x, eq = eq.th, newdata = newdata))$teta

} # X3 ok


if( is.null(x$X3) ){

sigma21 <- x$sigma21
sigma22 <- x$sigma22

nu1 <- x$nu1 
nu2 <- x$nu2

theta <- x$theta 
dof <- x$dof

                   }


} ############ ok 



if(missing(newdata)){

eta1 <- x$eta1
eta2 <- x$eta2

sigma21 <- x$sigma21
sigma22 <- x$sigma22

nu1 <- x$nu1 
nu2 <- x$nu2

theta <- x$theta 
dof   <- x$dof

}



eta1 <- as.numeric(eta1)
eta2 <- as.numeric(eta2)

sigma21 <- as.numeric(sigma21)
sigma22 <- as.numeric(sigma22)

nu1 <- as.numeric(nu1) 
nu2 <- as.numeric(nu2)

theta <- as.numeric(theta) 
dof   <- as.numeric(dof)



if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p1 <- as.numeric(distrHsAT(y1, eta1, sigma21, nu1, x$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc1)$p2) 
if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)) p2 <- as.numeric(distrHsAT(y2, eta2, sigma22, nu2, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2) 

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)) {y1rep <- rep(y1, length(eta1))
                                             p1pdf1 <- distrHsATDiscr(y1rep, eta1, sigma21, nu = 1, x$margins[1], x$VC$y1m, min.dn = min.pr, min.pr = min.pr, 
                                                                      max.pr = max.pr, left.trunc = x$VC$left.trunc1)
                                             p1pdf1$p2 <- as.numeric(p1pdf1$p2)
                                             p1pdf1$pdf2 <- as.numeric(p1pdf1$pdf2)
                                             p1 <- p1pdf1$p2 
                                             }
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)) {y2rep <- rep(y2, length(eta2))
                                             p2pdf2 <- distrHsATDiscr(y2rep, eta2, sigma22, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, 
                                                                      max.pr = max.pr, left.trunc = x$VC$left.trunc2)
                                             p2pdf2$p2 <- as.numeric(p2pdf2$p2)
                                             p2pdf2$pdf2 <- as.numeric(p2pdf2$pdf2)
                                             p2 <- p2pdf2$p2 
                                             }



######################################################################
######################################################################














if(cond == 0 || ( cond == 1 && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3) ) || (cond == 1 && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)) || (cond == 2 && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d))  ){  ############# COND 0 #############


if(x$BivD %in% x$BivD2){ ########## BivD2 ##########

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- C11 <- C01 <- C10 <- C00 <- NA
 
 
if(x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ #### CONT - CONT ####

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- mm(BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  )
if(length(theta) == 1) p12[x$teta.ind1] <- mm(BiCDF(p1[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )
                          }                       
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- mm(BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  )
if(length(theta) == 1) p12[x$teta.ind2] <- mm(BiCDF(p1[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )
                          }                            
                                                                                   } #### CONT - CONT ####




if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ #### DISCR - CONT ####

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- mm(BiCDF(p1pdf1$p2[x$teta.ind1], p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  ) - mm(BiCDF(mm(p1pdf1$p2[x$teta.ind1] - p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), p2[x$teta.ind1], nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  ); p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < min.pr, min.pr, p12[x$teta.ind1]) 
if(length(theta) == 1) p12[x$teta.ind1] <- mm(BiCDF(p1pdf1$p2[x$teta.ind1], p2[x$teta.ind1], nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )              - mm(BiCDF(mm(p1pdf1$p2[x$teta.ind1] - p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), p2[x$teta.ind1], nC1, theta,              dof), min.pr = min.pr, max.pr = max.pr  ); p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < min.pr, min.pr, p12[x$teta.ind1]) 
                          }                      
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- mm(BiCDF(p1pdf1$p2[x$teta.ind2], p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  ) - mm(BiCDF(mm(p1pdf1$p2[x$teta.ind2] - p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), p2[x$teta.ind2], nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  ); p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < min.pr, min.pr, p12[x$teta.ind2]) 
if(length(theta) == 1) p12[x$teta.ind2] <- mm(BiCDF(p1pdf1$p2[x$teta.ind2], p2[x$teta.ind2], nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )              - mm(BiCDF(mm(p1pdf1$p2[x$teta.ind2] - p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), p2[x$teta.ind2], nC2, theta,              dof), min.pr = min.pr, max.pr = max.pr  ); p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < min.pr, min.pr, p12[x$teta.ind2]) 
                          } 
                          
if(cond == 1) p12 <- p12/p1pdf1$pdf2
                          
                          
                                                                                    } #### DISCR - CONT ####






if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){ #### DISCR - DISCR ####


if( length(x$teta1) != 0){

if(length(theta) > 1){  

  C11[x$teta.ind1] <- mm(BiCDF(p1pdf1$p2[x$teta.ind1],                              p2pdf2$p2[x$teta.ind1],                              nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  )
  C01[x$teta.ind1] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), p2pdf2$p2[x$teta.ind1],                              nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  )
  C10[x$teta.ind1] <- mm(BiCDF(p1pdf1$p2[x$teta.ind1],                              mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1], min.pr, max.pr), nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  )
  C00[x$teta.ind1] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1], min.pr, max.pr), nC1, theta[x$teta.ind1], dof), min.pr = min.pr, max.pr = max.pr  )

  p12[x$teta.ind1] <- C11[x$teta.ind1] - C01[x$teta.ind1] - C10[x$teta.ind1] + C00[x$teta.ind1]
  p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < min.pr, min.pr, p12[x$teta.ind1])
                     }


if(length(theta) == 1){  

  C11[x$teta.ind1] <- mm(BiCDF(p1pdf1$p2[x$teta.ind1],                              p2pdf2$p2[x$teta.ind1],                              nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C01[x$teta.ind1] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), p2pdf2$p2[x$teta.ind1],                              nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C10[x$teta.ind1] <- mm(BiCDF(p1pdf1$p2[x$teta.ind1],                              mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1], min.pr, max.pr), nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C00[x$teta.ind1] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind1]-p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), mm(p2pdf2$p2[x$teta.ind1]-p2pdf2$pdf2[x$teta.ind1], min.pr, max.pr), nC1, theta, dof), min.pr = min.pr, max.pr = max.pr  )

  p12[x$teta.ind1] <- C11[x$teta.ind1] - C01[x$teta.ind1] - C10[x$teta.ind1] + C00[x$teta.ind1]
  p12[x$teta.ind1] <- ifelse(p12[x$teta.ind1] < min.pr, min.pr, p12[x$teta.ind1])

                     } 
                          }
                                                                                     
if( length(x$teta2) != 0){

if(length(theta) > 1){  

  C11[x$teta.ind2] <- mm(BiCDF(p1pdf1$p2[x$teta.ind2],                              p2pdf2$p2[x$teta.ind2],                              nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  )
  C01[x$teta.ind2] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), p2pdf2$p2[x$teta.ind2],                              nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  )
  C10[x$teta.ind2] <- mm(BiCDF(p1pdf1$p2[x$teta.ind2],                              mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2], min.pr, max.pr), nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  )
  C00[x$teta.ind2] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2], min.pr, max.pr), nC2, theta[x$teta.ind2], dof), min.pr = min.pr, max.pr = max.pr  )

  p12[x$teta.ind2] <- C11[x$teta.ind2] - C01[x$teta.ind2] - C10[x$teta.ind2] + C00[x$teta.ind2]
  p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < min.pr, min.pr, p12[x$teta.ind2])

                     }

if(length(theta) == 1){  

  C11[x$teta.ind2] <- mm(BiCDF(p1pdf1$p2[x$teta.ind2],                              p2pdf2$p2[x$teta.ind2],                              nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C01[x$teta.ind2] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), p2pdf2$p2[x$teta.ind2],                              nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C10[x$teta.ind2] <- mm(BiCDF(p1pdf1$p2[x$teta.ind2],                              mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2], min.pr, max.pr), nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C00[x$teta.ind2] <- mm(BiCDF(mm(p1pdf1$p2[x$teta.ind2]-p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), mm(p2pdf2$p2[x$teta.ind2]-p2pdf2$pdf2[x$teta.ind2], min.pr, max.pr), nC2, theta, dof), min.pr = min.pr, max.pr = max.pr  )

  p12[x$teta.ind2] <- C11[x$teta.ind2] - C01[x$teta.ind2] - C10[x$teta.ind2] + C00[x$teta.ind2]
  p12[x$teta.ind2] <- ifelse(p12[x$teta.ind2] < min.pr, min.pr, p12[x$teta.ind2])

                     }                     
                          } 


if(cond == 1) p12 <- p12/p1pdf1$pdf2
if(cond == 2) p12 <- p12/p2pdf2$pdf2


                                                                                       } #### DISCR - DISCR ####                          

} ########## BivD2 ########## ok 



if(!(x$BivD %in% x$BivD2)){ ########## !BivD2 ##########


if(x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3))  p12 <- mm(BiCDF(p1, p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  ) #### CONT - CONT ####


if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ #### DISCR - CONT ####

C1  <- mm(BiCDF(p1pdf1$p2, p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  ); C2  <- 0 
if(cumul == "no") C2  <- mm(BiCDF(mm(p1pdf1$p2 - p1pdf1$pdf2, min.pr, max.pr), p2, x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )  
p12 <- mm(C1 - C2, min.pr = min.pr, max.pr = max.pr  ) 

if(cond == 1) p12 <- p12/p1pdf1$pdf2

                                                                                    } #### DISCR - CONT ####



if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){ #### DISCR - DISCR ####

  C11 <- mm(BiCDF(p1pdf1$p2,                 p2pdf2$p2,                 x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C01 <- mm(BiCDF(mm(p1pdf1$p2-p1pdf1$pdf2, min.pr, max.pr), p2pdf2$p2,                 x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C10 <- mm(BiCDF(p1pdf1$p2,                 mm(p2pdf2$p2-p2pdf2$pdf2, min.pr, max.pr), x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )
  C00 <- mm(BiCDF(mm(p1pdf1$p2-p1pdf1$pdf2, min.pr, max.pr), mm(p2pdf2$p2-p2pdf2$pdf2, min.pr, max.pr), x$nC, theta, dof), min.pr = min.pr, max.pr = max.pr  )

  p12 <- mm(C11 - C01 - C10 + C00, min.pr, max.pr)
  
  if(cond == 1) p12 <- p12/p1pdf1$pdf2
  if(cond == 2) p12 <- p12/p2pdf2$pdf2
  
  } #### DISCR - DISCR #### 


                         } ########## !BivD2 ########## ok 





} ############# COND 0 ############# ok 










if(cond == 1){ ########### COND 1 ########### 

if(x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ ##### CONT - CONT ####

if(!(x$BivD %in% x$BivD2)) p12 <- copgHsCond(p1, p2, theta, dof = dof, x$BivD , min.pr = min.pr, max.pr = max.pr)$c.copula.be1


if(x$BivD %in% x$BivD2){

p12 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta[x$teta.ind1], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be1
if(length(theta) == 1) p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta, dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be1
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta[x$teta.ind2], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be1
if(length(theta) == 1) p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta, dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be1
                         }                            
                                                                           
                       }

                                                                                  } ##### CONT - CONT #### 

} ########### COND 1 ########### ok










if(cond == 2){#############


if(x$margins[1] %in% c(x$VC$m2, x$VC$m3) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ # CONT - CONT


if(!(x$BivD %in% x$BivD2)) p12 <- copgHsCond(p1, p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2


if(x$BivD %in% x$BivD2){

p12 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta[x$teta.ind1], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12[x$teta.ind1] <- copgHsCond(p1[x$teta.ind1], p2[x$teta.ind1], theta, dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta[x$teta.ind2], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12[x$teta.ind2] <- copgHsCond(p1[x$teta.ind2], p2[x$teta.ind2], theta, dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }                            
                                                                           
                       }

} # CONT - CONT ok





if(x$BivD %in% x$BivD2){ ########## BivD2 ##########

nC1 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop1),2] 
nC2 <- x$VC$ct[which(x$VC$ct[,1] == x$Cop2),2]

p12 <- C11 <- C01 <- C10 <- C00 <- NA
 
 
if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ #### DISCR - CONT ####


p12h1 <- p12h2 <- NA
 
if( length(x$teta1) != 0){
if(length(theta) > 1)  p12h1[x$teta.ind1] <- copgHsCond(p1pdf1$p2[x$teta.ind1], p2[x$teta.ind1], theta[x$teta.ind1], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12h1[x$teta.ind1] <- copgHsCond(p1pdf1$p2[x$teta.ind1], p2[x$teta.ind1], theta,              dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12h1[x$teta.ind2] <- copgHsCond(p1pdf1$p2[x$teta.ind2], p2[x$teta.ind2], theta[x$teta.ind2], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12h1[x$teta.ind2] <- copgHsCond(p1pdf1$p2[x$teta.ind2], p2[x$teta.ind2], theta,              dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }                            
                                                                           
                       #} 
 

if( length(x$teta1) != 0){
if(length(theta) > 1)  p12h2[x$teta.ind1] <- copgHsCond(mm(p1pdf1$p2[x$teta.ind1] - p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), p2[x$teta.ind1], theta[x$teta.ind1], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12h2[x$teta.ind1] <- copgHsCond(mm(p1pdf1$p2[x$teta.ind1] - p1pdf1$pdf2[x$teta.ind1], min.pr, max.pr), p2[x$teta.ind1], theta,              dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }  
                          
if( length(x$teta2) != 0){
if(length(theta) > 1)  p12h2[x$teta.ind2] <- copgHsCond(mm(p1pdf1$p2[x$teta.ind2] - p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), p2[x$teta.ind2], theta[x$teta.ind2], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
if(length(theta) == 1) p12h2[x$teta.ind2] <- copgHsCond(mm(p1pdf1$p2[x$teta.ind2] - p1pdf1$pdf2[x$teta.ind2], min.pr, max.pr), p2[x$teta.ind2], theta,              dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                         }                            
                                                                           
                       #} 
  
diffh1.h2 <- p12h1 - p12h2 
p12 <- ifelse(diffh1.h2 < min.pr, min.pr, diffh1.h2)  

  

                          
                                                                                    } #### DISCR - CONT ####

                    
} ########## BivD2 ##########



if(!(x$BivD %in% x$BivD2)){ ########## !BivD2 ##########


if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)){ #### DISCR - CONT ####

p12h1 <- copgHsCond(p1pdf1$p2,                   p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
p12h2 <- copgHsCond(mm(p1pdf1$p2 - p1pdf1$pdf2, min.pr, max.pr), p2, theta, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2

diffh1.h2 <- p12h1 - p12h2 
p12 <- ifelse(diffh1.h2 < min.pr, min.pr, diffh1.h2)                                  


                                                                                     } #### DISCR - CONT ####


                         } ########## !BivD2 ##########








}################### ok







# kendalls' tau

if(x$BivD %in% x$BivD2 && tau.res == TRUE)    {x$SemiParFit <- x; tau <- Reg2Copost(x$SemiParFit, x$VC, theta, tau.res = tau.res)$tau} 
if(!(x$BivD %in% x$BivD2) && tau.res == TRUE) tau <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = tau.res)$tau

if(x$BivD %in% x$BivD2)    {x$SemiParFit <- x; theta <- Reg2Copost(x$SemiParFit, x$VC, theta, tau.res = FALSE)$theta} 
if(!(x$BivD %in% x$BivD2)) theta <- theta2tau(x$VC$BivD, x$VC$nCa, theta, tau.res = FALSE)$theta





if(intervals == TRUE){



bs <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)  
lf <- length(x$coefficients)


#############  
# etas
#############  

# try with 1 number

if(!missing(newdata)){ X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix") 
                       X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix") }
                       
if( missing(newdata)){ X1 <- x$X1 
                       X2 <- x$X2 }                       


eta1s <- eta.tr( X1%*%t(bs[,1:x$X1.d2])                     , x$VC$margins[1]) 
eta2s <- eta.tr( X2%*%t(bs[,(x$X1.d2+1):(x$X1.d2+x$X2.d2)]) , x$VC$margins[2])

#############  
# thetas
#############  



if(  is.null(x$X3) ){ epds <- bs[, lf]
                      if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) ) dofs <- bs[, lf - 1] else dofs <- dof
                    } 
   
if( !is.null(x$X3) ){ 
  	
 if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){ 
  
  if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)){
                            X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix") 
                            X6 <- predict.SemiParBIV(x, eq = 6, newdata = newdata, type = "lpmatrix")               
                            }  
       if( missing(newdata)){X5 <- x$X5; X6 <- x$X6} 
 
       dofs <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
       epds <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
 
                                                                     }
 
  if((x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) || (x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) ){
  
       if(!missing(newdata)){
                             X6 <- predict.SemiParBIV(x, eq = 6, newdata = newdata, type = "lpmatrix")  
                             X7 <- predict.SemiParBIV(x, eq = 7, newdata = newdata, type = "lpmatrix")               
                            }
       
       if( missing(newdata)){X6 <- x$X6; X7 <- x$X7}   
  
       dofs <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
       epds <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
  
                                            }
  
  if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par){
  
       if(!missing(newdata)){
                             X7 <- predict.SemiParBIV(x, eq = 7, newdata = newdata, type = "lpmatrix")
                             X8 <- predict.SemiParBIV(x, eq = 8, newdata = newdata, type = "lpmatrix")               
                             }
       
       if( missing(newdata)){X7 <- x$X7; X8 <- x$X8}    
  
       dofs <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
       epds <- X8%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2+x$X8.d2)])
  
     }  
  
  
  
                }else{
                
  dofs <- dof   
  
  if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont1par){
  
  
       if(!missing(newdata)) X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X3 <- x$X3 
 
       epds <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)])
 
                                                                     }  
  
   if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)) X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X4 <- x$X4 
 
       epds <- X4%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)])
 
                                                                     }   
  
  
   if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont3par){
  
  
       if(!missing(newdata)) X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X5 <- x$X5 
 
       epds <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
 
                                                                     }   
  
  
       
       
  if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par){
  
  
       if(!missing(newdata)) X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X5 <- x$X5 
 
       epds <- X5%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2)])
 
                                                                     }
 
  if((x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) || (x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) ){
  
       if(!missing(newdata)) X6 <- predict.SemiParBIV(x, eq = 6, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X6 <- x$X6   
  
       epds <- X6%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2)])
  
                                                                                                                                            }
  
  if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par){
  
       if(!missing(newdata)) X7 <- predict.SemiParBIV(x, eq = 7, newdata = newdata, type = "lpmatrix")               
       if( missing(newdata)) X7 <- x$X7    
  
       epds <- X7%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2+x$X5.d2+x$X6.d2+x$X7.d2)])
  
                                                                    }                
                              
}  # is.null              
                
                
                
 
}
  	                        


est.RHOb <- teta.tr(x$VC, epds)$teta
if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3) ) dofs <- dof.tr(dofs)$vao else dofs <- dof 

   
######   
   
#if(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)) tc <- 1 else tc <- 0  
      
#############  
# sigmas
#############  

      if( is.null(x$X3) ) { # could remove the conditions below and make it more general but need to define
                            # XX.d2 and change the multiplication operator. For now we will leave it like this   
  
if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% c(cont2par,cont3par) ){ ps2 <- x$X1.d2 + x$X2.d2 + 1 }  
if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont2par ){ ps1 <- x$X1.d2 + x$X2.d2 + 1; ps2 <- x$X1.d2 + x$X2.d2 + 1 + 1 }
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){ ps1 <- x$X1.d2 + x$X2.d2 + 1; ps2 <- x$X1.d2 + x$X2.d2 + 1 + 1 }
if((x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par) || (x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par) ){ ps1 <- x$X1.d2 + x$X2.d2 + 1; ps2 <- x$X1.d2 + x$X2.d2 + 1 + 1}
      
        if(!(x$VC$margins[1] %in% cont1par) ) sigma2.1.star <- bs[, ps1] 
        if(!(x$VC$margins[2] %in% cont1par) ) sigma2.2.star <- bs[, ps2] 
                                
                                } ### is.null
  
  
      if( !is.null(x$X3) ) {


if(x$VC$margins[1] %in% c(cont2par,cont3par) && x$VC$margins[2] %in% c(cont2par,cont3par) )  {    
      
       if(!missing(newdata)){ X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")
                              X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X3 <- x$X3; X4 <- x$X4 }        
      
       sigma2.1.star <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 
       sigma2.2.star <- X4%*%t(bs[,(x$X1.d2+x$X2.d2+x$X3.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2+x$X4.d2)]) 

                                                                                              }


if(x$VC$margins[1] %in% c(cont1par) && x$VC$margins[2] %in% c(cont2par,cont3par) ){    
      
       if(!missing(newdata)){ X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X3 <- x$X3}        
      
       sigma2.2.star <- X3%*%t(bs[,(x$X1.d2+x$X2.d2+1):(x$X1.d2+x$X2.d2+x$X3.d2)]) 

                                                                                   }

}  ### is.null


    if(!(x$VC$margins[1] %in% cont1par) ) sigma21 <- esp.tr(sigma2.1.star, x$VC$margins[1])$vrb else sigma21 <- 1  
    if(!(x$VC$margins[2] %in% cont1par) ) sigma22 <- esp.tr(sigma2.2.star, x$VC$margins[2])$vrb else sigma22 <- 1   
 
 

 
#############  
# NUs
#############    
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3) )  {    pn1 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1 
                            pn2 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1 + 1
                            nu1.st <- bs[, pn1]    
                            nu2.st <- bs[, pn2]  } 
  
  if( !is.null(x$X3) ) {  
  
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")
                              X6 <- predict.SemiParBIV(x, eq = 6, newdata = newdata, type = "lpmatrix") }
       if( missing(newdata)){ X5 <- x$X5; X6 <- x$X6 }   
 
       nu1.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
       nu2.st <- X6%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2)])
  
                       }   
   
nu1 <- enu.tr(nu1.st, x$VC$margins[1])$vrb   
nu2 <- enu.tr(nu2.st, x$VC$margins[2])$vrb   
  
} 


if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  {  pn2 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1; nu2.st <- bs[, pn2]  } 
  

  if( !is.null(x$X3) ) { 
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X5 <- x$X5}     
  
       nu2.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
  
                        }
  
  nu2 <- enu.tr(nu2.st, x$VC$margins[2])$vrb   

} 



if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  {  pn2 <- x$X1.d2 + x$X2.d2 + 1 + 1; nu2.st <- bs[, pn2]  } 
  
  
  
    if( !is.null(x$X3) ) { 
  
       if(!missing(newdata)){ X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X4 <- x$X4}     
  
       nu2.st <- X4%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
                          }
     
  
  nu2 <- enu.tr(nu2.st, x$VC$margins[2])$vrb   

}




  
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par ){  
    
  if( is.null(x$X3) )  {  pn1 <- x$X1.d2 + x$X2.d2 + 1 + 1 + 1; nu1.st <- bs[, pn1]  } 
  
  
  
  if( !is.null(x$X3) ) {
  
       if(!missing(newdata)){ X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")}
       if( missing(newdata)){ X5 <- x$X5}    
  
       nu1.st <- X5%*%t(bs[,(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2)]) 
  
                       }
  
  nu1 <- enu.tr(nu1.st, x$VC$margins[1])$vrb  
   
}  


####################################################


if( is.null(x$X3) ){


if(!(x$BivD == "T" && x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3))) dofs <- dof

if(is.null(sigma21)) sigma21 <- 1
if(is.null(sigma22)) sigma22 <- 1
if(is.null(nu1))     nu1 <- 1
if(is.null(nu2))     nu2 <- 1


                   }



if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){

p1s  <- distrHsAT(y1, eta1s, sigma21, nu1, x$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc1)$p2
p2s  <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2

}


if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){

y1rep <- rep(y1, length(eta1s))
ppdf1 <- distrHsATDiscr(y1rep, eta1s, sigma21, nu = 1, x$margins[1], x$VC$y1m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc1) 
p1s   <- ppdf1$p2
pdf1s <- ppdf1$pdf2

p2s  <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2 

}


if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d,x$VC$m2d)){

y1rep <- rep(y1, length(eta1s))
ppdf1 <- distrHsATDiscr(y1rep, eta1s, sigma21, nu = 1, x$margins[1], x$VC$y1m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc1) 
p1s   <- ppdf1$p2
pdf1s <- ppdf1$pdf2

y2rep <- rep(y2, length(eta2s))
ppdf2 <- distrHsATDiscr(y2rep, eta2s, sigma22, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2) 
p2s   <- ppdf2$p2 
pdf2s <- ppdf2$pdf2


}









if(cond == 0 || (cond == 1 && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m2, x$VC$m3)) || (cond == 1 && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)) || (cond == 2 && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d))  ){##################################




if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){ # CONT - CONT




if(!(x$BivD %in% x$BivD2)) p12s <- mm(BiCDF(p1s, p2s, x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )

if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(eta1s)[1])

if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- mm(BiCDF(p1s[x$teta.ind1,], p2s[x$teta.ind1,], nC1,  est.RHOb[x$teta.ind1,]), min.pr = min.pr, max.pr = max.pr  )                  
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- mm(BiCDF(p1s[x$teta.ind2,], p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,]), min.pr = min.pr, max.pr = max.pr  )
                      
                        }

                                                                                                                          


                                                                                 } # CONT - CONT


if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){ # DISCR - CONT


if(x$VC$BivD %in% c("N","T")){ C1s <- mm(BiCDF(p1s, p2s,       x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  ); C2s <- 0 
                               if(cumul == "no") C2s <- mm(BiCDF(mm(p1s-pdf1s, min.pr, max.pr), p2s, x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
                               p12s <- ifelse(C1s - C2s < min.pr, min.pr, C1s - C2s)

                             }else{


if(x$BivD %in% x$BivD2){

p12s <- C1s <- C2s <- matrix(NA, ncol = n.sim, nrow = dim(eta1s)[1])

if( length(x$teta1) != 0){

C1s[x$teta.ind1,]  <- mm(BiCDF(p1s[x$teta.ind1,],                     p2s[x$teta.ind1,], nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
C2s[x$teta.ind1,]  <- mm(BiCDF(mm(p1s[x$teta.ind1,]-pdf1s[x$teta.ind1,], min.pr = min.pr, max.pr = max.pr  ), p2s[x$teta.ind1,], nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
p12s[x$teta.ind1,] <- ifelse(C1s[x$teta.ind1,] - C2s[x$teta.ind1,] < min.pr, min.pr, C1s[x$teta.ind1,] - C2s[x$teta.ind1,])

                         }                            
 
if( length(x$teta2) != 0){

C1s[x$teta.ind2,]  <- mm(BiCDF(p1s[x$teta.ind2,],                     p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
C2s[x$teta.ind2,]  <- mm(BiCDF(mm(p1s[x$teta.ind2,]-pdf1s[x$teta.ind2,], min.pr = min.pr, max.pr = max.pr  ), p2s[x$teta.ind2,], nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
p12s[x$teta.ind2,] <- ifelse(C1s[x$teta.ind2,] - C2s[x$teta.ind2,] < min.pr, min.pr, C1s[x$teta.ind2,] - C2s[x$teta.ind2,])

                         }  
                            
                            
                      
                        } # Biv2



if(!(x$BivD %in% x$BivD2)){ C1s <- mm(BiCDF(p1s,       p2s, x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  ); C2s <- 0
                            if(cumul == "no") C2s <- mm(BiCDF(mm(p1s-pdf1s,min.pr = min.pr, max.pr = max.pr  ), p2s, x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
                            p12s <- ifelse(C1s - C2s < min.pr, min.pr, C1s - C2s)
                          }

                                   } # else


if(cond == 1) p12s <- p12s/pdf1s

                                                                                  }# DISCR - CONT











if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d,x$VC$m2d)){ # DISCR - DISCR


if(x$VC$BivD %in% c("N","T")){ 

  C11s <- mm(BiCDF(p1s,       p2s,       x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C01s <- mm(BiCDF(mm(p1s-pdf1s, min.pr = min.pr, max.pr = max.pr), p2s,       x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C10s <- mm(BiCDF(p1s, mm(p2s-pdf2s, min.pr = min.pr, max.pr = max.pr),       x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C00s <- mm(BiCDF(mm(p1s-pdf1s, min.pr = min.pr, max.pr = max.pr), mm(p2s-pdf2s, min.pr = min.pr, max.pr = max.pr), x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )

 p12s <- C11s - C01s - C10s + C00s

                              }else{


if(x$BivD %in% x$BivD2){

p12s <- C11s <- C01s <- C10s <- C00s <- matrix(NA, ncol = n.sim, nrow = dim(eta1s)[1])

if( length(x$teta1) != 0){

  C11s[x$teta.ind1,] <- mm(BiCDF(p1s[x$teta.ind1,],                     p2s[x$teta.ind1,],                     nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C01s[x$teta.ind1,] <- mm(BiCDF(mm(p1s[x$teta.ind1,]-pdf1s[x$teta.ind1,],min.pr = min.pr, max.pr = max.pr  ), p2s[x$teta.ind1,],                     nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C10s[x$teta.ind1,] <- mm(BiCDF(p1s[x$teta.ind1,],                     mm(p2s[x$teta.ind1,]-pdf2s[x$teta.ind1,],min.pr = min.pr, max.pr = max.pr  ), nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C00s[x$teta.ind1,] <- mm(BiCDF(mm(p1s[x$teta.ind1,]-pdf1s[x$teta.ind1,],min.pr = min.pr, max.pr = max.pr  ), mm(p2s[x$teta.ind1,]-pdf2s[x$teta.ind1,],min.pr = min.pr, max.pr = max.pr  ), nC1, est.RHOb[x$teta.ind1,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )

 p12s[x$teta.ind1,] <- C11s[x$teta.ind1,] - C01s[x$teta.ind1,] - C10s[x$teta.ind1,] + C00s[x$teta.ind1,]  
 p12s[x$teta.ind1,] <- ifelse(p12s[x$teta.ind1,] < min.pr, min.pr, p12s[x$teta.ind1,]) 

                        }


if( length(x$teta2) != 0){

  C11s[x$teta.ind2,] <- mm(BiCDF(p1s[x$teta.ind2,],                     p2s[x$teta.ind2,],                     nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C01s[x$teta.ind2,] <- mm(BiCDF(mm(p1s[x$teta.ind2,]-pdf1s[x$teta.ind2,],min.pr = min.pr, max.pr = max.pr  ), p2s[x$teta.ind2,],                     nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C10s[x$teta.ind2,] <- mm(BiCDF(p1s[x$teta.ind2,],                     mm(p2s[x$teta.ind2,]-pdf2s[x$teta.ind2,],min.pr = min.pr, max.pr = max.pr  ), nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C00s[x$teta.ind2,] <- mm(BiCDF(mm(p1s[x$teta.ind2,]-pdf1s[x$teta.ind2,],min.pr = min.pr, max.pr = max.pr  ), mm(p2s[x$teta.ind2,]-pdf2s[x$teta.ind2,],min.pr = min.pr, max.pr = max.pr  ), nC2, -est.RHOb[x$teta.ind2,], dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )

 p12s[x$teta.ind2,] <- C11s[x$teta.ind2,] - C01s[x$teta.ind2,] - C10s[x$teta.ind2,] + C00s[x$teta.ind2,]  
 p12s[x$teta.ind2,] <- ifelse(p12s[x$teta.ind2,] < min.pr, min.pr, p12s[x$teta.ind2,]) 

                         }



                      }



if(!(x$BivD %in% x$BivD2)){

                           
  C11s <- mm(BiCDF(p1s, p2s,             x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C01s <- mm(BiCDF(mm(p1s-pdf1s,min.pr = min.pr, max.pr = max.pr  ), p2s,       x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C10s <- mm(BiCDF(p1s, mm(p2s-pdf2s,min.pr = min.pr, max.pr = max.pr  ),       x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )
  C00s <- mm(BiCDF(mm(p1s-pdf1s,min.pr = min.pr, max.pr = max.pr  ), mm(p2s-pdf2s,min.pr = min.pr, max.pr = max.pr  ), x$nC, est.RHOb, dofs, test = FALSE), min.pr = min.pr, max.pr = max.pr  )

 p12s <- C11s - C01s - C10s + C00s  
 p12s <- ifelse(p12s < min.pr, min.pr, p12s)  


                          }


                                     } # else

if(cond == 1) p12s <- p12s/pdf1s 
if(cond == 2) p12s <- p12s/pdf2s 



} # DISCR - DISCR





}###################





      
if(cond == 1){ # only for CONT - CONT

if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){

if(!(x$BivD %in% x$BivD2)) p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = dofs, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be1


if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])
 
if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- copgHsCond(p1s[x$teta.ind1,], p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = dofs, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be1                                               
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- copgHsCond(p1s[x$teta.ind2,], p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = dofs, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be1
                                                                                                     
                       }

}


}




if(cond == 2){


if(x$margins[1] %in% c(x$VC$m2,x$VC$m3) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){ # CONT - CONT

if(!(x$BivD %in% x$BivD2)) p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = dofs, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2



if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])
 
if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- copgHsCond(p1s[x$teta.ind1,], p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = dofs, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2                                               
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- copgHsCond(p1s[x$teta.ind2,], p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = dofs, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2
                                                                                                     
                       }

                                                                                 } # CONT - CONT
                                                                                 
                                                                                 
                                                                                 

if(x$margins[1] %in% c(x$VC$m1d,x$VC$m2d) && x$margins[2] %in% c(x$VC$m2,x$VC$m3)){ # DISCR - CONT
                                                                                                                                                        
if(!(x$BivD %in% x$BivD2)) p12s <- copgHsCond(p1s, p2s, est.RHOb, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2 - copgHsCond(mm(p1s - pdf1s, min.pr, max.pr), p2s, est.RHOb, dof = dof, x$BivD, min.pr = min.pr, max.pr = max.pr)$c.copula.be2



if(x$BivD %in% x$BivD2){

p12s <- matrix(NA, ncol = n.sim, nrow = dim(p1s)[1])
                                                                                                                                                                                                                                                                             
if( length(x$teta1) != 0) p12s[x$teta.ind1,] <- copgHsCond(p1s[x$teta.ind1,], p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2 - copgHsCond(mm((p1s - pdf1s)[x$teta.ind1,], min.pr, max.pr) , p2s[x$teta.ind1,],  est.RHOb[x$teta.ind1,], dof = dof, x$Cop1, min.pr = min.pr, max.pr = max.pr)$c.copula.be2                                                
if( length(x$teta2) != 0) p12s[x$teta.ind2,] <- copgHsCond(p1s[x$teta.ind2,], p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2 - copgHsCond(mm((p1s - pdf1s)[x$teta.ind2,], min.pr, max.pr) , p2s[x$teta.ind2,], -est.RHOb[x$teta.ind2,], dof = dof, x$Cop2, min.pr = min.pr, max.pr = max.pr)$c.copula.be2 
                                                                                                     
                       }
                       
p12s <- ifelse( p12s < min.pr, min.pr, p12s ) 
                       

                                                                                 } # DISCR - CONT

                                                                                                                                  
} # cond 2








nCa   <- x$VC$nCa
BivDt <- x$VC$BivD

  if(x$BivD %in% x$BivD2){
  
  if(x$BivD %in% x$BivD2[c(1:4,13:16)]) { BivDt <- "C0"; nCa <- 3} 
  if(x$BivD %in% x$BivD2[5:8]) { BivDt <- "J0"; nCa <- 6}
  if(x$BivD %in% x$BivD2[9:12]){ BivDt <- "G0"; nCa <- 4}
  
                         }
  
ass.msR <- theta2tau(BivDt, nCa, est.RHOb, tau.res = tau.res)
if(tau.res == TRUE) taus <- ass.msR$tau
thetas <- ass.msR$theta

                                
                                                  
if(tau.res == TRUE){ taus <- matrix(taus, 1, n.sim); CIkt <- rowQuantiles(taus, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)}

thetas <- matrix(thetas, 1, n.sim)
CItheta <- rowQuantiles(thetas, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

#if( is.null(x$X3) ) CIkt <- t(CIkt) 


if(tau.res == TRUE){
 if(x$BivD %in% x$BivD2){ 
 
   if(length(x$theta) > 1){
 
     if( length(x$teta2) != 0) CIkt[x$teta.ind2, ] <- -CIkt[x$teta.ind2, ]; CIkt[x$teta.ind2, c(1,2)] <- CIkt[x$teta.ind2, c(2,1)] 
                                 
                          }else{
 
     if( length(x$teta2) != 0) CIkt <- -CIkt; CIkt[, c(1,2)] <- CIkt[, c(2,1)]
                                 
                                }
 }
}

 if(x$BivD %in% x$BivD2){ 
 
   if(length(x$theta) > 1){
 
     if( length(x$teta2) != 0) CItheta[x$teta.ind2, ] <- -CItheta[x$teta.ind2, ]; CItheta[x$teta.ind2, c(1,2)] <- CItheta[x$teta.ind2, c(2,1)] 
                                 
                          }else{
 
     if( length(x$teta2) != 0) CItheta <- -CItheta; CItheta[, c(1,2)] <- CItheta[, c(2,1)]
                                 
                                }
 }







} # interv






                        }## biv
                        
######################################################################################################
######################################################################################################





######################################################################################################
######################################################################################################

      if(type == "independence"){



if(!missing(newdata)){

nu1 <- nu2 <- sigma21 <- sigma22 <- NA

eta1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$coefficients[1:x$X1.d2]
eta2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[1:x$X2.d2]

if( !is.null(x$X3) ){ 


if(x$margins[1] %in% cont1par && x$margins[2] %in% c(cont2par,cont3par))             sigma22 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[(x$X2.d2+1):(x$X2.d2+x$X3.d2)], x$margins[2])$vrb
if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)) sigma22 <- esp.tr(predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[(x$X2.d2+1):(x$X2.d2+x$X4.d2)], x$margins[2])$vrb
if(!(x$margins[1] %in% cont1par))                                                    sigma21 <- esp.tr(predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$coefficients[(x$X1.d2+1):(x$X1.d2+x$X3.d2)], x$margins[1])$vrb


if(x$margins[1] %in% cont1par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 4  ; indn1 <- NA; indn2 <- (x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont3par){ eq.nu1 <- 5;  eq.nu2 <- 6  ; indn1 <- (x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2); indn2 <- (x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X6.d2)}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont2par){ eq.nu1 <- NA; eq.nu2 <- NA ; indn1 <- NA; indn2 <- NA}
if(x$margins[1] %in% cont2par && x$margins[2] %in% cont3par){ eq.nu1 <- NA; eq.nu2 <- 5  ; indn1 <- NA; indn2 <- (x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X5.d2)}
if(x$margins[1] %in% cont3par && x$margins[2] %in% cont2par){ eq.nu1 <- 5;  eq.nu2 <- NA ; indn1 <- (x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2); indn2 <- NA}

if(x$margins[1] %in% cont3par) nu1 <- enu.tr(predict.SemiParBIV(x, eq = eq.nu1, newdata = newdata, type = "lpmatrix")%*%x$gamlss1$coefficients[indn1], x$margins[1])$vrb
if(x$margins[2] %in% cont3par) nu2 <- enu.tr(predict.SemiParBIV(x, eq = eq.nu2, newdata = newdata, type = "lpmatrix")%*%x$gamlss2$coefficients[indn2], x$margins[2])$vrb

}

if( is.null(x$X3) ){ 

sigma21 <- x$gamlss1$sigma2
sigma22 <- x$gamlss2$sigma2

nu1 <- x$gamlss1$nu
nu2 <- x$gamlss2$nu

}

}



if(missing(newdata)){


eta1 <- x$gamlss1$eta1
eta2 <- x$gamlss2$eta1

sigma21 <- x$gamlss1$sigma2
sigma22 <- x$gamlss2$sigma2

nu1 <- x$gamlss1$nu
nu2 <- x$gamlss2$nu

}






if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p1 <- p1a <- as.numeric(distrHsAT(y1, eta1, sigma21, nu1, x$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc1)$p2) 
if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)) p2 <- p2a <- as.numeric(distrHsAT(y2, eta2, sigma22, nu2, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2) 

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           y1rep <- rep(y1, length(eta1))
                                           p1pdf1 <- distrHsATDiscr(y1rep, eta1, sigma21, nu = 1, x$margins[1], x$VC$y1m, min.dn = min.pr, min.pr = min.pr, 
                                                                    max.pr = max.pr, left.trunc = x$VC$left.trunc1)
                                           p1     <- as.numeric(p1pdf1$p2) # as.numeric(p1pdf1$p2 - p1pdf1$pdf2)
                                           p1a    <- as.numeric(p1pdf1$pdf2)
                                        
                                           }
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           y2rep <- rep(y2, length(eta2))
                                           p2pdf2 <- distrHsATDiscr(y2rep, eta2, sigma22, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, 
                                                                    max.pr = max.pr, left.trunc = x$VC$left.trunc2)
                                           p2     <- as.numeric(p2pdf2$p2) # as.numeric(p2pdf2$p2 - p2pdf2$pdf2)
                                           p2a    <- as.numeric(p2pdf2$pdf2)
                                           }  


p12 <- p1a*p2a

if(cond == 1) p12 <- p2a
if(cond == 2) p12 <- p1a



if(intervals == TRUE){

bs1 <- rMVN(n.sim, mean = x$gamlss1$coefficients, sigma=x$gamlss1$Vb)
bs2 <- rMVN(n.sim, mean = x$gamlss2$coefficients, sigma=x$gamlss2$Vb)


#############  
# etas
############# 

if(!missing(newdata)){ X1 <- predict.SemiParBIV(x, eq = 1, newdata = newdata, type = "lpmatrix")
                       X2 <- predict.SemiParBIV(x, eq = 2, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X1 <- x$X1; X2 <- x$X2}  

eta1s <- eta.tr( X1%*%t(bs1[,1:x$X1.d2]), x$VC$margins[1]) 
eta2s <- eta.tr( X2%*%t(bs2[,1:x$X2.d2]), x$VC$margins[2]) 

#############  
# sigmas
#############  

      if( is.null(x$X3) ) {
      
      
      if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)){
      
       sigma2.1.star <- bs1[, x$X1.d2 + 1] 
       sigma2.2.star <- bs2[, x$X2.d2 + 1] 
      
      }
      
      if(x$margins[1] %in% c(cont1par) && x$margins[2] %in% c(cont2par,cont3par)){
       
       sigma2.2.star <- bs2[, x$X2.d2 + 1] 
      
      }      
      
      
                          }
  
  
  
  
      if( !is.null(x$X3) ) {

if(x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)){

if(!missing(newdata)){ X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")
                       X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X3 <- x$X3; X4 <- x$X4}  

       sigma2.1.star <- X3%*%t(bs1[,(x$X1.d2+1):(x$X1.d2+x$X3.d2)]) 
       sigma2.2.star <- X4%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X4.d2)]) 
                                                                                     }
                                                                                     
if(x$margins[1] %in% c(cont1par) && x$margins[2] %in% c(cont2par,cont3par)){

if(!missing(newdata)) X3 <- predict.SemiParBIV(x, eq = 3, newdata = newdata, type = "lpmatrix")  
if( missing(newdata)) X3 <- x$X3 

       sigma2.2.star <- X3%*%t(bs2[,(x$X2.d2+1):(x$X2.d2+x$X3.d2)]) 
                                                                            }                                                                                     
                           }  

if(!(x$margins[1] %in% c(cont1par))) sigma21 <- esp.tr(sigma2.1.star, x$VC$margins[1])$vrb   
                                     sigma22 <- esp.tr(sigma2.2.star, x$VC$margins[2])$vrb  

#############  
# NUs
#############    
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont3par ){  
    
  if( is.null(x$X3) )  {     
       nu1.st <- bs1[, x$X1.d2 + 1 + 1]  
       nu2.st <- bs2[, x$X2.d2 + 1 + 1]     
                        } 
  
  if( !is.null(x$X3) ) {  
  
  
if(!missing(newdata)){ X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")
                       X6 <- predict.SemiParBIV(x, eq = 6, newdata = newdata, type = "lpmatrix")}
if( missing(newdata)){ X5 <- x$X5; X6 <- x$X6}    
  
       nu1.st <- X5%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2)]) 
       nu2.st <- X6%*%t(bs2[,(x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X6.d2)])
                       }   
   
nu1 <- enu.tr(nu1.st, x$VC$margins[1])$vrb   
nu2 <- enu.tr(nu2.st, x$VC$margins[2])$vrb   
  
} 


if(x$VC$margins[1] %in% cont2par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  nu2.st <- bs2[, x$X2.d2 + 1 + 1]  
  
  
  if( !is.null(x$X3) ){
  
  
    if(!missing(newdata)){ X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X5 <- x$X5}  
  
    nu2.st <- X5%*%t(bs2[,(x$X2.d2 + x$X4.d2 + 1):(x$X2.d2 + x$X4.d2 + x$X5.d2)]) 
  
  
                      }
  
  
                       nu2    <- enu.tr(nu2.st, x$VC$margins[2])$vrb   
} 



if(x$VC$margins[1] %in% cont1par && x$VC$margins[2] %in% cont3par ){  
  
  if( is.null(x$X3) )  nu2.st <- bs2[, x$X2.d2 + 1 + 1]   
  
  
  if( !is.null(x$X3) ){
  
  
    if(!missing(newdata)){ X4 <- predict.SemiParBIV(x, eq = 4, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X4 <- x$X4}  
  
    nu2.st <- X4%*%t(bs2[,(x$X2.d2 + x$X3.d2 + 1):(x$X2.d2 + x$X3.d2 + x$X4.d2)]) 
  
  
                      }
  
   nu2    <- enu.tr(nu2.st, x$VC$margins[2])$vrb   
} 


  
  
  
if(x$VC$margins[1] %in% cont3par && x$VC$margins[2] %in% cont2par ){  
    
  if( is.null(x$X3) )  nu1.st <- bs1[, x$X1.d2 + 1 + 1]   


  if( !is.null(x$X3) ){
    
    if(!missing(newdata)){ X5 <- predict.SemiParBIV(x, eq = 5, newdata = newdata, type = "lpmatrix")}
    if( missing(newdata)){ X5 <- x$X5}    
  
    nu1.st <- X5%*%t(bs1[,(x$X1.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X3.d2 + x$X5.d2)])                  
  
                      }
  
  nu1    <- enu.tr(nu1.st, x$VC$margins[1])$vrb    
}  


#if( is.null(x$X3) ){
#
#if(x$VC$margins[1] %in% cont2par) sigma21  <- matrix(rep(sigma21, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
#if(x$VC$margins[2] %in% cont2par) sigma22  <- matrix(rep(sigma22, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
#if(x$VC$margins[1] %in% cont3par) nu1      <- matrix(rep(nu1, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
#if(x$VC$margins[2] %in% cont3par) nu2      <- matrix(rep(nu2, each = dim(eta1s)[1]), ncol = n.sim, byrow=FALSE)
#
#                   }




if(x$margins[1] %in% c(x$VC$m2, x$VC$m3)) p1s <- p1sa <- distrHsAT(y1, eta1s, sigma21, nu1, x$margins[1], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc1)$p2
if(x$margins[2] %in% c(x$VC$m2, x$VC$m3)) p2s <- p2sa <- distrHsAT(y2, eta2s, sigma22, nu2, x$margins[2], min.dn = min.pr, min.pr = min.pr, max.pr = max.pr, left.trunc = x$VC$left.trunc2)$p2

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           y1rep <- rep(y1, length(eta1s))
                                           p1pdf1s <- distrHsATDiscr(y1rep, eta1s, sigma21, nu = 1, x$margins[1], x$VC$y1m, min.dn = min.pr, min.pr = min.pr, 
                                                                     max.pr = max.pr, left.trunc = x$VC$left.trunc1)
                                           p1s     <- p1pdf1s$p2 # p1pdf1s$p2 - p1pdf1s$pdf2, dim(eta1s)[1], n.sim)
                                           p1sa    <- p1pdf1s$pdf2
                                           }
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)){ 
                                           y2rep <- rep(y2, length(eta2s))
                                           p2pdf2s <- distrHsATDiscr(y2rep, eta2s, sigma22, nu = 1, x$margins[2], x$VC$y2m, min.dn = min.pr, min.pr = min.pr, 
                                                                     max.pr = max.pr, left.trunc = x$VC$left.trunc2)
                                           p2s     <- p2pdf2s$p2 # p2pdf2s$p2 - p2pdf2s$pdf2, dim(eta1s)[1], n.sim)
                                           p2sa    <- p2pdf2s$pdf2
                                           }  

p12s <- p1sa*p2sa

if(cond == 1) p12s <- p2sa
if(cond == 2) p12s <- p1sa



} # intervals


} # independence



list(p12 = p12, p12s = matrix(p12s, 1, length(p12s)), p1 = p1, p2 = p2, p3 = NULL, CItau = CIkt, tau = tau, theta = theta, CItheta = CItheta)


}



