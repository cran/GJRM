SATE <- function(x, trt, surv.t = NULL, int.var = NULL, joint = TRUE,
                 n.sim = 100, prob.lev = 0.05, ls = 10, plot.type = "none", ...){

####################

if(joint == TRUE)  type <- "joint"
if(joint == FALSE) type <- "univariate"
#if(x$Model == "ROY" && type == "univariate") stop("You need to fit a separate univariate model to obtain the effect(s) of interest.")
if(!(plot.type %in% c("none", "sate", "survival"))) stop("plot.type can be set to none, sate or survival.")
if(x$VC$margins[1] != "probit" || x$VC$margins[2] != "probit") stop("Currently, the margins have to be Gaussian for this measure to make sense from a causal point of view.") 
#if(x$surv == FALSE && x$end.surv == FALSE) stop("This function is not suitable for the fitted model.")
if(x$surv == FALSE) stop("This function is not suitable for the fitted model.")

if(missing(trt) && x$Model != "ROY") stop("You must provide the name of the binary treatment variable.")

if(!missing(trt)){ if(is.character(trt) == FALSE && x$Model != "ROY") stop("trt must be a character.")}

if(type == "univariate" && x$gamlssfit == FALSE && x$Model != "ROY") stop("You need to fit the univariate model to obtain the effect(s) of interest. \nRefit the model and set uni.fit = TRUE.")

d0 <- d1 <- d01 <- X2s <- X3s <- 1 # just for rm() later

####################



if( !is.null(int.var) && !is.list(int.var) ){

  if( length(int.var) != 2 )              stop("int.var must contain a name and a value for the interaction term.")
  if( is.character(int.var[1]) == FALSE ) stop("The first element of int.var must be the name of the interaction term.")

  int.var1 <- int.var[1]
  int.var2 <- as.numeric(int.var[2]) # as.numeric works for both numeric and factor vars

  if( !(int.var2 %in% c(0, 1)) ) stop("The interaction can only currently take value 0 or 1.")  

}


if( !is.null(int.var) && is.list(int.var) ){

int.var1 <- int.var2 <- NA

 for(i in 1:length(int.var)){

  if( length(int.var[[i]]) != 2 )              stop("Each of the int.var objects in the list must contain a name and a value for the interaction term.")
  if( is.character(int.var[[i]][1]) == FALSE ) stop("The first element of each of the int.var objects must be the name of the interaction term.")

  int.var1[i] <- int.var[[i]][1]
  int.var2[i] <- as.numeric(int.var[[i]][2]) # as.numeric works for both numeric and factor vars

  if( !(int.var2[i] %in% c(0, 1)) ) stop("The interaction term can only currently take value 0 or 1.")  
 }


}


##########################

SATEs <- S0s <- S1s <- NA
n.t <- as.character(x$formula[[2]][[2]])

if(is.null(surv.t)){ if(x$Model == "ROY") ry2 <- range(c(x$y2, x$y3)) else ry2 <- range(x$y2) # no need for 0.0001 correction here
                     surv.t <- seq(ry2[1], ry2[2], length.out = ls)

                    }

  if(!is.null(surv.t)){
  
    if(!is.numeric(surv.t)) stop("You must provide a numeric value or vector for surv.t.")
    surv.t <- ifelse(surv.t < 0.0001, 0.0001, surv.t)
  
  }
  
CIs <- CIs0 <- CIs1 <- matrix(NA, length(surv.t), 2)  
  

##########################



if(x$Model == "ROY"){

if( !( type %in% c("joint") ) ) stop("Error in parameter type value. You can only choose joint.")

    bs <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
    mono.sm.pos <- c(x$VC$X1.d2 + x$VC$mono.sm.pos2, x$VC$X1.d2 + x$VC$X2.d2 + x$VC$mono.sm.pos3)
    bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos])


    indp1   <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
    indp2   <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2)
    coefeq2 <- x$coef.t[indp1]
    coefeq3 <- x$coef.t[indp2]
    


    ff <- reformulate(all.vars(x$gam2$terms)[-1])
    d01 <- model.frame(ff, data = eval(x$mcd), subset = eval(x$sbs)) # this implies that x$call$data is always there, not needed as we must have it
    attr(d01,"terms") <- NULL

for(i in 1:length(surv.t)){

   d01[, n.t] <- surv.t[i]


    X2s <- try(predict(x$gam2, newdata = d01, type = "lpmatrix"), silent = TRUE)
    if(any(class(X2s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
      
    X3s <- try(predict(x$gam3, newdata = d01, type = "lpmatrix"), silent = TRUE)
    if(any(class(X3s)=="try-error")) stop("Check that the factor variables' levels\nin the selected sample for the second margin are the same as those in the complete dataset.\nRead the Details section in ?gjrm for more details.") 
  

if( !is.null(int.var) && !is.list(int.var)){

   if( any(grepl(int.var1, dimnames(X2s)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   #if( any(grepl(":", int.var1)) == FALSE      )           stop("Check the name provided for the interaction term.")

   if( int.var2 == 0) X2s[, int.var1] <- X3s[, int.var1] <- 0
   if( int.var2 == 1) X2s[, int.var1] <- X3s[, int.var1] <- 1

                       }



if( !is.null(int.var) && is.list(int.var)) {

 for(j in 1:length(int.var)){

   if( any(grepl(int.var1[j], dimnames(X2s)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   #if( any(grepl(":", int.var1[j])) == FALSE      )           stop("Check the name provided for the interaction term.")

if( int.var2[j] == 0) X2s[, int.var1[j]] <- X3s[, int.var1[j]] <- 0
if( int.var2[j] == 1) X2s[, int.var1[j]] <- X3s[, int.var1[j]] <- 1

 }

}

eta0 <- X2s%*% coefeq2
eta1 <- X3s%*% coefeq3

S0   <- probmS(eta0, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
S1   <- probmS(eta1, x$VC$margins[3], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr

SATEs[i] <- mean(S1, na.rm = TRUE) - mean(S0, na.rm = TRUE) 
S0s[i]   <- mean(S0, na.rm = TRUE) 
S1s[i]   <- mean(S1, na.rm = TRUE) 

eta0s <- X2s %*% t(bs[, indp1])
eta1s <- X3s %*% t(bs[, indp2])

S0sim   <- probmS(eta0s, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
S1sim   <- probmS(eta1s, x$VC$margins[3], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr

sim.SATE <- colMeans(S1sim, na.rm = TRUE) - colMeans(S0sim, na.rm = TRUE) 
sim.S0   <- colMeans(S0sim, na.rm = TRUE) 
sim.S1   <- colMeans(S1sim, na.rm = TRUE) 

 
CIs[i, ]  <- as.numeric(quantile(sim.SATE, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))
CIs0[i, ] <- as.numeric(quantile(sim.S0,   c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))
CIs1[i, ] <- as.numeric(quantile(sim.S1,   c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))


}

                                                                     

    } # end ROY
    
    
    
    
    ########################














if(x$Model != "ROY"){




###############################
# SET UP of some key quantities
###############################

if(type == "joint"){

    mono.sm.pos <- c(x$VC$mono.sm.pos1, x$VC$mono.sm.pos2 + x$VC$X1.d2)
    bs          <- rMVN(n.sim, mean = x$coef.t, sigma = x$Vb.t)
    indp        <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
    coefeq2     <- x$coef.t[indp]    
    
    }
    
if(type == "univariate"){
    
    mono.sm.pos <- x$gamlss2$VC$mono.sm.pos
    bs          <- rMVN(n.sim, mean = x$gamlss2$coef.t, sigma = x$gamlss2$Vb.t) 
    indp        <- 1:length(x$gamlss2$coefficients)
    coefeq2     <- x$gamlss2$coef.t   
    
    }
    
bs[, mono.sm.pos] <- ifelse(bs[, mono.sm.pos] < 0, 0, bs[, mono.sm.pos])


ff <- reformulate(all.vars(x$gam2$terms)[-1])
d0 <- d1 <- model.frame(ff, data = eval(x$mcd), subset = eval(x$sbs)) # this implies that x$call$data is always there, not needed as we must have it
attr(d0,"terms") <- attr(d1,"terms") <- NULL

if( is.logical(d0[, trt]) == TRUE) stop("The treatment variable must be a binary numeric or factor variable.")

d0[, trt] <- 0 
d1[, trt] <- 1 

##################################
# EFFECT CALCULATION and INTERVALS
##################################

for(i in 1:length(surv.t)){

d0[, n.t] <- d1[, n.t] <- surv.t[i]

d00 <- predict(x$gam2, newdata = d0, type = "lpmatrix")
d11 <- predict(x$gam2, newdata = d1, type = "lpmatrix")




if( !is.null(int.var) && !is.list(int.var)){

   if( any(grepl(int.var1, dimnames(d11)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1)) == FALSE      )           stop("Check the name provided for the interaction term.")

   if( int.var2 == 0) d11[, int.var1] <- 0
   if( int.var2 == 1) d11[, int.var1] <- 1

                       }



if( !is.null(int.var) && is.list(int.var)) {

 for(j in 1:length(int.var)){

   if( any(grepl(int.var1[j], dimnames(d11)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1[j])) == FALSE      )           stop("Check the name provided for the interaction term.")

if( int.var2[j] == 0) d11[, int.var1[j]] <- 0
if( int.var2[j] == 1) d11[, int.var1[j]] <- 1

 }

}







eta0 <- d00%*% coefeq2
eta1 <- d11%*% coefeq2

S0   <- probmS(eta0, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
S1   <- probmS(eta1, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr

SATEs[i] <- mean(S1, na.rm = TRUE) - mean(S0, na.rm = TRUE) 
S0s[i]   <- mean(S0, na.rm = TRUE) 
S1s[i]   <- mean(S1, na.rm = TRUE) 

eta0s <- d00 %*% t(bs[, indp])
eta1s <- d11 %*% t(bs[, indp])

S0sim   <- probmS(eta0s, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
S1sim   <- probmS(eta1s, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr

sim.SATE <- colMeans(S1sim, na.rm = TRUE) - colMeans(S0sim, na.rm = TRUE) 
sim.S0   <- colMeans(S0sim, na.rm = TRUE) 
sim.S1   <- colMeans(S1sim, na.rm = TRUE) 

 
CIs[i, ]  <- as.numeric(quantile(sim.SATE, c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))
CIs0[i, ] <- as.numeric(quantile(sim.S0,   c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))
CIs1[i, ] <- as.numeric(quantile(sim.S1,   c(prob.lev/2, 1 - prob.lev/2), na.rm = TRUE))


}




} # end != ROY


##########
## PLOT ##

if(plot.type == "survival"){


    Survival <- as.matrix(cbind(S1s, CIs1[, 1], CIs1[, 2], S0s, CIs0[, 1], CIs0[, 2]))
    if(dim(Survival)[1] < 2) stop("Survival plot not available for one time value.")
    attr(Survival, "dimnames")[[1]] <- surv.t

    matplot(surv.t, Survival, type = "l", col = c(rep("black", 3), rep("grey", 3)), lty = c(1, 2, 2, 1, 2, 2), ...)

                            }
 


 if(plot.type == "sate"){
 
     Sate <- as.matrix(cbind(CIs[, 1], CIs[, 2], SATEs))
     attr(Sate, "dimnames")[[1]] <- surv.t
 
     matplot(surv.t, Sate, pch = c(1, 1, 20), col = c("white", "white","black"), ...)
     for(i in 1:length(surv.t)) lines(c(surv.t[i], surv.t[i]), CIs[i, ])
 
                        }


######################################################################

rm(bs, d0, d1, d01, X2s, X3s) 

lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "")

res <- as.data.frame(cbind(surv.t, SATEs, CIs[, 1], CIs[, 2])); names(res) <- c("surv.t", "SATE", lbn, ubn)

out <- list(res = res, prob.lev = prob.lev, sim.SATE = sim.SATE, type = type)
 							 
class(out) <- "SATE"

out


}

 