jc.probs <- function(x, y1, y2, y3 = NULL, newdata, type = "joint", cond = 0, intervals = FALSE, n.sim = 100, prob.lev = 0.05, 
                     epsilon = 1e-07, cumul = "no"){

######################################################################################################
# preliminary checks
######################################################################################################

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

if(missing(y1)) stop("You must provide a value for y1.")
if(missing(y2)) stop("You must provide a value for y2.")
if(x$triv == TRUE && missing(y3)) stop("You must provide a binary value for y3.")


if(!(type %in% c("joint","independence"))) stop("Error in parameter type value. It should be one of: joint, independence.")
if(!(cond %in% c(0,1,2, 3))) stop("Error in parameter cond value. It should be one of: 0, 1, 2 or 3 (for the trivariate model).")

if( type %in% c("independence") && x$VC$gamlssfit == FALSE && is.null(x$VC$K1)) stop("You need to re-fit the model and set gamlssfit = TRUE to obtain probabilities under independence.") # CopulaCLM


if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && (!is.wholenumber(y1) || y1 < 0)) stop("The value for y1 must be discrete and positive.")
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d) && (!is.wholenumber(y2) || y2 < 0)) stop("The value for y2 must be discrete and positive.")


if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !(y1 %in% c(0,1)) && is.null(x$VC$K1) ) stop("The value for y1 must be either 0 or 1.")

if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !is.null(x$VC$K1) && !(y1 %in% seq.int(1, x$VC$K1))) stop(paste(paste("The value for y1 must be an integer between 1 and", x$VC$K1, "")), ".") # CopulaCLM

if( x$VC$Cont == "NO" && x$margins[2] %in% bin.link){ if( !(y1 %in% c(0,1)) || !(y2 %in% c(0,1))   ) stop("The value for y1 and/or y2 must be either 0 or 1.") }

if(x$triv == FALSE) {if(!missing(newdata) && x$BivD %in% x$BivD2) stop("Prediction for models based on mixed copulae and a new dataset is not feasible.")}



if(x$triv == TRUE){

if( !(y1 %in% c(0,1)) ) stop("The value for y1 must be either 0 or 1.")
if( !(y2 %in% c(0,1)) ) stop("The value for y2 must be either 0 or 1.")
if( !(y3 %in% c(0,1)) ) stop("The value for y3 must be either 0 or 1.")

}






######################################################################################################
######################################################################################################
if(x$triv == FALSE){

if(x$VC$Cont == "YES" && x$surv == FALSE )              rr <- jc.probs1(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon, cumul)


if(x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) ){ 


                                                                                                                                                                                  #
 if( is.null(x$VC$K1)) rr <- jc.probs2(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon)                                        #
                                                                                                                                                                                  #
 if(!is.null(x$VC$K1)) {                                                                                                                                                          #
	if( type %in% c("joint")        && x$VC$ind.ord == "TRUE") stop("You need to provide the fitted joint model as input to obtain predictive probabilities.")                #
	if( type %in% c("independence") && x$VC$ind.ord != "TRUE") stop("You need to provide the fitted independence model as input to obtain probabilities under independence.") #
	rr <- jc.probs7(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon) # CopulaCLM                                           #
 }                                                                                                                                                                                # CopulaCLM

}

if(x$VC$Cont == "NO" &&   x$margins[2] %in% bin.link)   rr <- jc.probs3(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon)    

if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% bin.link && x$margins[2] %in% bin.link )             rr <- jc.probs4(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon)
if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% bin.link ) rr <- jc.probs5(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon)
# if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par,cont3par) && x$margins[2] %in% bin.link ) rr <- jc.probs5(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon)
}

if(x$triv == TRUE)  rr <- jc.probs6(x, y1, y2, y3, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link, epsilon)

# copulaReg
# SemiParBIV and copulaSampleSel
# SemiParBIV

######################################################################################################
######################################################################################################

p12s <- rr$p12s
p12  <- rr$p12
p1   <- rr$p1
p2   <- rr$p2

if(x$triv == TRUE){ 

p3 <- rr$p3 

theta12  <- rr$theta12
theta13  <- rr$theta13
theta23  <- rr$theta23 
theta12s <- rr$theta12s
theta13s <- rr$theta13s
theta23s <- rr$theta23s

}



######################################################################################################
######################################################################################################

if(intervals == TRUE){ 

CIp12 <- rowQuantiles(p12s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

if(x$triv == TRUE){ # triv

   if(type == "joint"){
     CItheta12 <- rowQuantiles(theta12s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     CItheta13 <- rowQuantiles(theta13s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)
     CItheta23 <- rowQuantiles(theta23s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
                      }else{
             CItheta12 <- CItheta13 <- CItheta23 <- cbind(0, 0)
                            }
                  } # triv


if(x$triv == FALSE){ # biv

  if(length(p12) > 1)  {res <- data.frame(p12, CIp12, p1, p2);       names(res)[2:3] <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2)))}
  if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, p1, p2))); names(res) <- c("p12",names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))),"p1","p2")}
                  
                   } # biv



if(x$triv == TRUE){

  if(length(p12) > 1){ p123 <- p12
                       CIp123 <- CIp12 
                       res <- data.frame(p123, CIp123, p1, p2, p3, theta12, theta13, theta23, CItheta12, CItheta13, CItheta23 )
                       names(res)[c(2:3, 10:15)] <- c(rep(names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))),4))
                     }
                      
                      
  if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, p1, p2, p3, theta12, theta13, theta23, CItheta12, CItheta13, CItheta23)))
                        names(res) <- c("p123",names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))),"p1","p2","p3",
                                        "theta12", "theta13", "theta23", c(rep(names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))),3))  )
                        }
                      
  rm(p12, CIp12)                      
                  }

}else{

if(x$triv == FALSE) res <- data.frame(p12, p1, p2)
if(x$triv == TRUE)  {p123 <- p12; res <- data.frame(p123, p1, p2, p3, theta12, theta13, theta23)}

}

res.n <- names(res) # CopulaCLM


if(x$triv == FALSE && !is.null(rr$tau) && !is.null(rr$CIkt)){ 

if(is.null(dim(rr$CIkt)))  {CItauLB <- rr$CIkt[1]; CItauUB <- rr$CIkt[2]}
if(!is.null(dim(rr$CIkt))) {CItauLB <- rr$CIkt[,1]; CItauUB <- rr$CIkt[,2]}


res <- data.frame(res, tau = rr$tau, CItauLB = CItauLB, CItauUB = CItauUB )     
	names(res)[1 : length(res.n)] <- res.n # CopulaCLM

}

if(x$triv == FALSE && !is.null(rr$tau) &&  is.null(rr$CIkt)) res <- data.frame(res, tau = rr$tau )     


return(res)



}



