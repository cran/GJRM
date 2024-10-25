SATE <- function(x, trt, surv.t = NULL, int.var = NULL, joint = TRUE,
                 n.sim = 100, prob.lev = 0.05, ls = 10, plot.type = "survival", ...){


####################
# CHECKS
####################

#if(is.null(surv.t)) full <- TRUE else full <- FALSE
#if(is.null(plot.type)) plot <- FALSE else{plot <- TRUE; full <- TRUE}
#if(plot == TRUE && plot.type == "sate")     plot.type <- 2
#if(plot == TRUE && plot.type == "survival") plot.type <- 1
#d0 <- x$gam2$model
#d0 <- d0[, !grepl("(weights)", names(d0))]
#attr(attributes(x$gam2$terms)$dataClasses, "names")[1] <- n.t
#names(d0)[1] <- n.t
#d0 <- d1 <- d0[, !grepl(n.t, names(d0))]
#names(x$gam2$model)[1:dim(d0)[2]] <- names(d0)
#ls(envir = sys.frame(-1))
#model.matrix(x$gam2), data = environment(object)
#eval(model.frame(ff, data = newdata, na.action = na.act), parent.frame())
#if( !( type %in% c("univariate","joint") ) ) stop("Error in parameter type value. It should be one of: univariate or joint.")


surv.tO <- surv.t 

ry2 <- range(x$y2)
Sates <- SATEs <- S0s <- S1s <- NA
CIs <- CIs0 <- CIs1 <- matrix(NA, ls, 2)

if(joint == TRUE)  type <- "joint"
if(joint == FALSE) type <- "univariate"


if(!(plot.type %in% c("sate", "survival"))) stop("plot.type can be set to sate or survival.")


if(x$VC$margins[1] != "probit" || x$VC$margins[2] != "probit") stop("Currently, the margins have to be probit and -probit for this measure to make sense from a causal point of view.") 
if(x$surv == FALSE && x$end.surv == FALSE) stop("This function is not suitable for the model you have fitted.")
if(missing(trt)) stop("You must provide the name of the binary treatment variable.")
if(is.character(trt) == FALSE) stop("trt must be a character.")



if( !is.null(int.var) ){

  if( length(int.var) != 2 )              stop("int.var must contain a name and a value for the interaction variable.")
  if( is.character(int.var[1]) == FALSE ) stop("The first element of int.var must be the name of the interaction.")

  int.var1 <- int.var[1]
  int.var2 <- as.numeric(int.var[2]) # as.numeric works for both numeric and factor vars

  if( !(int.var2 %in% c(0, 1)) ) stop("The interaction can only currently take value 0 or 1.")  

}



if(type == "univariate" && x$gamlssfit == FALSE) stop("You need to fit the univariate model to obtain the effect of interest. Refit the model and set uni.fit = TRUE.")


if(!is.null(surv.t)){

  if(length(surv.t) > 1) stop("You must provide a numeric value for surv.t.")
  
  surv.t <- ifelse(surv.t < 0.0001, 0.0001, surv.t)
  
  if( !(surv.t > ry2[1] &&  surv.t < ry2[2]) ) stop("The value of surv.t must be within the range of the observed time variable.")
}




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


if(is.null(surv.t)) surv.t <- seq(ry2[1], ry2[2], length.out = ls)
n.t <- as.character(x$formula[[2]][[2]])


ff <- reformulate(all.vars(x$gam2$terms)[-1])
d0 <- d1 <- model.frame(ff, data = get(x$call$data)) # this implies that x$call$data is always there, not needed as we must have it
attr(d0,"terms") <- attr(d1,"terms") <- NULL


if( is.logical(d0[, trt]) == TRUE) stop("The treatment variable must be a binary numeric or factor variable.")


#if( !is.null(int.var) && int.var2 == 0){ 
#
#    #int.varo <- d0[, int.var1] 
#
#    d0[, int.var1] <- int.var2 
#    d1[, int.var1] <- int.var2 
#
#}#p.d11 <- which(gT)



d0[, trt] <- 0 
d1[, trt] <- 1 


##################################
# EFFECT CALCULATION and INTERVALS
##################################

for(i in 1:length(surv.t)){

d0[, n.t] <- d1[, n.t] <- surv.t[i]

d00 <- predict(x$gam2, newdata = d0, type = "lpmatrix")
d11 <- predict(x$gam2, newdata = d1, type = "lpmatrix")

if( !is.null(int.var) ){

   if( any(grepl(int.var1, dimnames(d11)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1)) == FALSE      )           stop("Check the name provided for the interaction term.")

if( int.var2 == 0) d11[, int.var1] <- 0
if( int.var2 == 1) d11[, int.var1] <- 1

}





eta0 <- d00%*% coefeq2
eta1 <- d11%*% coefeq2

S0   <- probmS(eta0, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr
S1   <- probmS(eta1, x$VC$margins[2], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr

Sates[i] <- SATEs[i] <- mean(S1, na.rm = TRUE) - mean(S0, na.rm = TRUE) 
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

##########
## PLOT ##


#if(plot == TRUE && full == FALSE){ # not really a good option, will be removed as in ATE
#  
#  hist(sim.SATE, freq = FALSE, ylim = c(0, max(density(sim.SATE)$y, hist(sim.SATE, plot = FALSE)$density)), ...)
#  lines(density(sim.SATE))
#
#                    }


if(plot.type == "survival" && is.null(surv.tO)){

  #if(plot.type == 1){


    Survival <- as.matrix(cbind(S1s, CIs1[, 1], CIs1[, 2], S0s, CIs0[, 1], CIs0[, 2]))
    attr(Survival, "dimnames")[[1]] <- surv.t

    matplot(surv.t, Survival, type = "l", col = c(rep("black", 3), rep("grey", 3)), lty = c(1, 2, 2, 1, 2, 2), ...)

 #legend("bottomright", legend = c("treatment = 1", "treatment = 0"),
 #       col = c("black", "grey"), lty = c(1, 1), inset = c(0, 1), xpd = TRUE, horiz = TRUE, bty = "n") 
 
 }
 

 

 if(plot.type == "sate" && is.null(surv.tO)){
 
     SATEs <- as.matrix(cbind(CIs[, 1], CIs[, 2], SATEs))
     attr(SATEs, "dimnames")[[1]] <- surv.t
 
     matplot(surv.t, SATEs, pch = c(1, 1, 20), col = c("white", "white","black"), ...)
     for(i in 1:length(surv.t)) lines(c(surv.t[i], surv.t[i]), CIs[i, ])
 
 }


#}

######################################################################

rm(bs, d0, d1) 

lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "")

if(!is.null(surv.tO)) {res <-    c(CIs[1, 1], SATEs, CIs[1, 2]); names(res) <- c(lbn, "SATE", ubn)}
if(is.null(surv.tO))  {res <- as.data.frame(cbind(surv.t, Sates, CIs[, 1], CIs[, 2])); names(res) <- c("surv.t", "SATE", lbn, ubn)}

out <- list(res = res, prob.lev = prob.lev, sim.SATE = sim.SATE, type = type, full = is.null(surv.tO))
 							 

class(out) <- "SATE"

out


}

 