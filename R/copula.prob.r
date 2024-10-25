copula.prob <- function(x, y1, y2, y3 = NULL, newdata, joint = TRUE, cond = 0, intervals = FALSE, 
                        n.sim = 100, prob.lev = 0.05, theta = FALSE, tau = FALSE, 
                        min.pr = 1e-323, max.pr = 1){

######################################################################################################
# preliminary checks
######################################################################################################

if(missing(newdata))                               stop("You must provide a data frame.")
if(!missing(newdata)){ 
                       if(!is.data.frame(newdata)) stop("You must provide a data frame.")
                       if(dim(newdata)[1] > 1)     stop("This calculation is supported for single row data frames.")
                     }

#if(!(type %in% c("surv", "hazard", "cumhaz"))) stop("The type argument can either be surv, hazard or cumhaz")

if(x$surv == TRUE && x$end.surv == TRUE) stop("Function not ready yet for the fitted model. Get in touch to check progress.")

if(n.sim < 2) stop("n.sim has to be a value greater than 2.")

sf <- "surv" 

if(joint == TRUE) type <- "joint" else type <- "independence" 

if(x$VC$Model == "ROY") stop("This function is not designed for the type of model chosen for modelling. Get in touch for more info.")

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

cont1par <- x$VC$m1d   
cont2par <- c(x$VC$m2,x$VC$m2d) 
cont3par <- x$VC$m3 
bin.link <- x$VC$bl

if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% c(cont2par,cont3par)) stop("This function is currently not suitable for survival models with\nparametric margins. Consider using semiparametric margins instead. ")


if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% bin.link && x$margins[2] %in% bin.link) y1 <- y2 <- 1
if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% bin.link ) y2 <- 1


if(x$univar.gamlss == TRUE) stop("This function is not suitable for univariate models.")
#if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)) stop("This function is currently not suitable for models involving discrete margins.")
#if(x$VC$Cont == "YES" && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)) stop("This function is currently not suitable for models involving a discrete margin.")

if(x$margins[2] == "TW") stop("Function not extended yet to the case of Tweedie margin. \nGet in touch for more info.")

if(missing(y1)) stop("You must provide a value for y1.")
if(missing(y2)) stop("You must provide a value for y2.")
if(x$triv == TRUE && is.null(y3)) stop("You must provide a binary value for y3.")
if(x$triv == FALSE && !is.null(y3)) stop("You can not provide a value for y3.")


if(missing(newdata)){

if(length(y1) != 1) stop("You can only provide one value for y1.")
if(length(y2) != 1) stop("You can only provide one value for y2.")

if( length(y1) == 1 && !is.null(x$ordinal) ) y1 <- rep(y1, length(x$eta1))         
if( length(y2) == 1 && !is.null(x$ordinal) ) y2 <- rep(y2, length(x$eta1))         


if( !is.null(y3) ){  

        if( length(y3) != 1 ) stop("You can only provide one value for y3.")    
        if( length(y3) == 1 && !is.null(x$ordinal) ) y3 <- rep(y3, length(x$eta1))                 
                  }

}



if(!missing(newdata)){

if(length(y1) != 1) stop("You can only provide one value for y1.")
if(length(y2) != 1) stop("You can only provide one value for y2.")

if( length(y1) == 1 && !is.null(x$ordinal)) y1 <- rep(y1, dim(newdata)[1])         
if( length(y2) == 1 && !is.null(x$ordinal)) y2 <- rep(y2, dim(newdata)[1])         


if( !is.null(y3) ){  

        if( length(y3) != 1 ) stop("You can only provide one value for y3.")    
        if( length(y3) == 1 && !is.null(x$ordinal)) y3 <- rep(y3, dim(newdata)[1])                 
                  }



}





if(!(type %in% c("joint","independence"))) stop("Error in parameter type value. It should be one of: joint, independence.")
if(!(cond %in% c(0,1,2, 3))) stop("Error in parameter cond value. It should be one of: 0, 1, 2 or 3 (for the trivariate model).")

#if( type %in% c("independence") && x$VC$gamlssfit == FALSE && is.null(x$VC$K1)) stop("You need to re-fit the model and set uni.fit = TRUE to obtain probabilities under independence.") # CopulaCLM
if( type %in% c("independence") && x$VC$gamlssfit == FALSE) stop("You need to re-fit the model and set uni.fit = TRUE to obtain probabilities under independence.") # This is now possible in CopulaCLM

if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && (!is.wholenumber(y1) || y1 < 0)) stop("The value for y1 must be discrete and positive.")
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d) && (!is.wholenumber(y2) || y2 < 0)) stop("The value for y2 must be discrete and positive.")


if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !all(y1 %in% c(0,1)) && is.null(x$VC$K1) ) stop("The value for y1 must be either 0 or 1.")

if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !is.null(x$VC$K1) && !all(y1 %in% seq.int(1, x$VC$K1))) stop(paste(paste("The value for y1 must be an integer between 1 and", x$VC$K1, "")), ".") # CopulaCLM

if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !is.null(x$VC$K2) && !all(y2 %in% seq.int(1, x$VC$K2))) stop(paste(paste("The value for y2 must be an integer between 1 and", x$VC$K2, "")), ".") # CopulaCLM

if( x$VC$Cont == "NO" && x$margins[2] %in% bin.link && is.null(x$VC$K1) ){ if( !all(y1 %in% c(0,1)) || !all(y2 %in% c(0,1))   ) stop("The value for y1 and/or y2 must be either 0 or 1.") } # CopulaCLM

if(x$triv == FALSE) {if(!missing(newdata) && x$BivD %in% x$BivD2) stop("Prediction for models based on mixed copulae and a new dataset is not feasible.")}



if(x$triv == TRUE){

if( !(y1 %in% c(0,1)) ) stop("The value for y1 must be either 0 or 1.")
if( !(y2 %in% c(0,1)) ) stop("The value for y2 must be either 0 or 1.")
if( !(y3 %in% c(0,1)) ) stop("The value for y3 must be either 0 or 1.")

}



lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "")


######################################################################################################
######################################################################################################
if(x$triv == FALSE){

if(x$VC$Cont == "YES" && x$surv == FALSE )              rr <- jc.probs1(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, cumul = "no", tau.res = tau)


if(x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) ){ # CopulaCLM: This selects the ordinal-continuous model

                                                                                                                                                                                 
 if( is.null(x$VC$K1)) rr <- jc.probs2(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)                       
                                                                                                                                                                                 
 if(!is.null(x$VC$K1)) {                                                                                                                                       
	#if( type %in% c("joint")        && x$VC$ind.ord == "TRUE") stop("You need to provide the fitted joint model as input to obtain predictive probabilities.")               
	
	if( type %in% c("independence") && x$VC$is_ordcon == "TRUE" && x$VC$gamlssfit == FALSE ) stop("You have to set uni.fit = TRUE when fitting the model to obtain probabilities under independence.")

	rr <- jc.probs7(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau) # CopulaCLM
 }

}

if(x$VC$Cont == "NO" && x$margins[2] %in% bin.link) { # CopulaCLM: This selects the ordinal-ordinal model
	if ( is.null(x$VC$K1)) rr <- jc.probs3(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)    
	if (!is.null(x$VC$K1)) rr <- jc.probs8(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)  
}

if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% bin.link && x$margins[2] %in% bin.link )             rr <- jc.probs4(x, y1, y2, newdata, type, cond, sf, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)
if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% bin.link ) rr <- jc.probs5(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)
# if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% bin.link ) rr <- jc.probs5(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)
}

if(x$triv == TRUE)  rr <- jc.probs6(x, y1, y2, y3, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, min.pr, max.pr, tau.res = tau)

# copulaReg
# SemiParBIV and copulaSampleSel
# SemiParBIV

######################################################################################################
###################################################################################################### 


  p12  <- rr$p12
  p1   <- rr$p1
  p2   <- rr$p2



if(x$triv == TRUE){ 

   p3 <- rr$p3 

   theta12  <- rr$theta12
   theta13  <- rr$theta13
   theta23  <- rr$theta23 
   
   p123 <- p12

}


quant.names <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))) 



   if(x$triv == FALSE && !is.null(rr$theta) && !is.null(rr$CItheta)){

	if(is.null(dim(rr$CItheta))){ 
	                             CIthetaLB <- rr$CItheta[1]; CIthetaUB <- rr$CItheta[2]
	                             if(tau == TRUE){ CItauLB <- rr$CItau[1]; CItauUB <- rr$CItau[2]}	                             
	                             }
	
	if(!is.null(dim(rr$CItheta))){
	                              CIthetaLB <- rr$CItheta[,1]; CIthetaUB <- rr$CItheta[,2]
	                              if(tau == TRUE){ CItauLB <- rr$CItau[,1]; CItauUB <- rr$CItau[,2]}
	                              }
   }


######################################################################################################
###################################################################################################### 

if(intervals == TRUE){ 


CIp12 <- rowQuantiles(rr$p12s, probs = c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE)
# above is for both TRIV and BIV



if(x$triv == TRUE){ # triv

    	temp.names1 <- c( rep( quant.names, 4) ) 
    	temp.names2 <- c("", "", "theta12", "theta12", "theta13", "theta13", "theta23", "theta23") # c("p123", "p123", "theta12", "theta12", "theta13", "theta13", "theta23", "theta23")
    	p.tempnames <- paste(temp.names2, temp.names1)
    
    
   if(type == "joint"){ # joint
   
     CItheta12 <- rowQuantiles(rr$theta12s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     CItheta13 <- rowQuantiles(rr$theta13s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     CItheta23 <- rowQuantiles(rr$theta23s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
     

  if(length(p123) > 1){ 
                       res <- data.frame(p123, CIp12, theta12, CItheta12, theta13, CItheta13, theta23, CItheta23)#, p1, p2, p3)
                       names(res)[c(2,3,5,6,8,9,11,12)] <- p.tempnames                      
                     }
                      
                      
  if(length(p123) == 1) {
                        res <- data.frame(t(c(p123, CIp12, theta12, CItheta12, theta13, CItheta13, theta23, CItheta23)))#, p1, p2, p3)))
                        names(res)[c(1,4,7,10)] <- c("p123", "theta12", "theta13", "theta23")#,"p1","p2","p3")
                        names(res)[c(2,3,5,6,8,9,11,12)] <- c(p.tempnames); dimnames(res)[[1]] <- ""
                        }
                        
                        
                      }  # joint
                        
               
   if(type != "joint"){ # indep.
      

  if(length(p123) > 1){ 
                       res <- data.frame(p123, CIp12)#, p1, p2, p3)
                       names(res)[c(2, 3)] <- c(paste(lbn),paste(ubn)) # c(paste("p123",lbn),paste("p123",ubn)) # c("p123 2.5%", "p123 97.5%")  #                     
                     }
                      
                      
  if(length(p123) == 1) {
                        res <- data.frame(t(c(p123, CIp12)))          #, p1, p2, p3)))
                        names(res) <- c("p123", paste(lbn), paste(ubn)) # c("p123", paste("p123",lbn), paste("p123",ubn)) # c("p123", "p123 2.5%", "p123 97.5%")       #, "p1","p2","p3")
                        dimnames(res)[[1]] <- ""
                        }
                        
                        
                      }  # !joint               
               
               
  rm(p12, CIp12)                      
                      
                  } # triv 






if(x$triv == FALSE){ # biv - BPO is a bit of repetition but clear for now


  temp.names1 <- quant.names
  temp.names2 <- c("p12", "p12")
  p.tempnames <- temp.names1 # paste(temp.names2, temp.names1) 
  long.names  <- c("p12", paste(lbn), paste(ubn), "theta", paste("theta",lbn), paste("theta",ubn)) # c("p12", paste("p12",lbn), paste("p12",ubn), "theta", paste("theta",lbn), paste("theta",ubn)) # c("p12", "p12 2.5%", "p12 97.5%", "theta", "theta 2.5%", "theta 97.5%")   #, "p1", "p2")
  long.names2 <- c("p12", paste(lbn), paste(ubn), "theta", paste("theta",lbn), paste("theta",ubn),"tau", paste("tau",lbn), paste("tau",ubn)) # c("p12", paste("p12",lbn), paste("p12",ubn), "theta", paste("theta",lbn), paste("theta",ubn),"tau", paste("tau",lbn), paste("tau",ubn)) # c("p12", "p12 2.5%", "p12 97.5%", "theta", "theta 2.5%", "theta 97.5%", "tau", "tau 2.5%", "tau 97.5%")   #, "p1", "p2")


  if(x$Model != "BPO0"){

  if(is.null(rr$theta) && is.null(rr$CItheta)){  

  if(length(p12) > 1)  {res <- data.frame(p12, CIp12);       names(res)[c(2,3)] <- p.tempnames }
  if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12))); names(res) <- c("p12", p.tempnames); dimnames(res)[[1]] <- "" }
 
                                              }
  
  if(!is.null(rr$theta) && !is.null(rr$CItheta) && tau == FALSE){  

  if(length(p12) > 1)  {res <-     data.frame(p12, CIp12, theta = rr$theta, CIthetaLB = CIthetaLB, CIthetaUB = CIthetaUB)
                        names(res) <- long.names}
  if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, theta = rr$theta, CIthetaLB = CIthetaLB, CIthetaUB = CIthetaUB)))#, p1, p2))) 
                        names(res) <- long.names; dimnames(res)[[1]] <- ""}  
                      
                                                                 }
  

  if(!is.null(rr$theta) && !is.null(rr$CItheta) && tau == TRUE){  

  if(length(p12) > 1)  {res <-     data.frame(p12, CIp12, theta = rr$theta, CIthetaLB = CIthetaLB, CIthetaUB = CIthetaUB, tau = rr$tau, CItauLB = CItauLB, CItauUB = CItauUB)#, p1, p2)    
                        names(res) <- long.names2}
  if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, theta = rr$theta, CIthetaLB = CIthetaLB, CIthetaUB = CIthetaUB, tau = rr$tau, CItauLB = CItauLB, CItauUB = CItauUB)))#, p1, p2)))
                        names(res) <- long.names2; dimnames(res)[[1]] <- ""}
                        
                                                                }
                        } # BPO0
                        
  if(x$Model == "BPO0"){ # repeats lies above but avoids tau conditions etc
  
       #if(length(p12) > 1)  {res <- data.frame(p12, CIp12, p1, p2);       names(res)[c(2,3)] <- p.tempnames }
       #if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, p1, p2))); names(res) <- c("p12", p.tempnames, "p1", "p2"); dimnames(res)[[1]] <- "" }
       
       if(length(p12) > 1)  {res <- data.frame(p12, CIp12);       names(res)[c(2,3)] <- p.tempnames }
       if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12))); names(res) <- c("p12", p.tempnames); dimnames(res)[[1]] <- "" }  
                        }
                        
    
                   } # biv 






}else{ # intervals = TRUE - bit above

if(x$Model != "BPO0"){
if(x$triv == FALSE &&  is.null(rr$theta) &&  is.null(rr$CItheta) && (tau == FALSE || tau == TRUE)) res <- data.frame(p12)#, p1, p2)
if(x$triv == FALSE && !is.null(rr$theta) &&  is.null(rr$CItheta) && tau == FALSE) res <- data.frame(p12, theta = rr$theta)#, p1, p2)
if(x$triv == FALSE && !is.null(rr$theta) &&  is.null(rr$CItheta) && tau == TRUE)  res <- data.frame(p12, theta = rr$theta, tau = rr$tau)#, p1, p2)
                     }

if(x$Model == "BPO0") res <- data.frame(p12)#, p1, p2)



if(x$triv == TRUE && type == "joint")  res <- data.frame(p123, theta12, theta13, theta23)#, p1, p2, p3)
if(x$triv == TRUE && type != "joint")  res <- data.frame(p123)#, p1, p2, p3)



if(dim(res)[1] == 1) dimnames(res)[[1]] <- ""



}

if(theta == FALSE){

pt <- which( grepl("theta", names(res)) )

if(identical(pt, integer(0)) == FALSE) res <- res[-pt]

}

#return(res)

res



}



