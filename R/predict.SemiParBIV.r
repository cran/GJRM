predict.SemiParBIV<- function(object, eq, ...){

if(missing(eq)) stop("You must provide the equation number.")

if(eq > object$l.flist) stop("The fitted model has a smaller number of equations.") 

                         
 if(eq==1){ ss.pred <- object$gam1
            ind <- 1:object$X1.d2 
            } 
 if(eq==2){ ss.pred <- object$gam2
            
            #if(object$VC$informative == "yes") ind <- object$X1.d2 + (1:dim(object$VC$fgam$X)[2])
            #else 
            
            ind <- (object$X1.d2+1):(object$X1.d2+object$X2.d2) 
            
            }
 if(eq==3){ ss.pred <- object$gam3
            ind <- (object$X1.d2+object$X2.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2) }   

 if(eq==4){ ss.pred <- object$gam4
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2) }             

 if(eq==5){ ss.pred <- object$gam5
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2) }             
  
 if(eq==6){ ss.pred <- object$gam6
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2) }             

 if(eq==7){ ss.pred <- object$gam7
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2) }             

 if(eq==8){ ss.pred <- object$gam8
            ind <- (object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2+1):(object$X1.d2+object$X2.d2+object$X3.d2+object$X4.d2+object$X5.d2+object$X6.d2+object$X7.d2+object$X8.d2) }             
                                             
                                   
           ss.pred$coefficients   <- object$coefficients[ind]
           ss.pred$coefficients.t <- object$coef.t[ind]
           ss.pred$Vp   <- object$Vb[ind,ind]
           ss.pred$Vp.t <- object$Vb.t[ind,ind]
           ss.pred$sig2 <- 1
           ss.pred$scale.estimated <- FALSE 

  #predict.gam(ss.pred, ...)
  predict(ss.pred, ...)


}
