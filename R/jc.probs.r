jc.probs <- function(x, y1, y2, newdata, type = "bivariate", cond = 0, intervals = FALSE, n.sim = 100, prob.lev = 0.05){

######################################################################################################
# preliminary checks
######################################################################################################

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol

cont1par <- x$VC$m1d   
cont2par <- c(x$VC$m2,x$VC$m2d) 
cont3par <- x$VC$m3 
bin.link <- x$VC$bl

if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% c(cont2par) ) stop("This function is currently not suitable for survival models with\nparametric margins. Consider using semiparametric margins instead. ")


if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% bin.link ) y1 <- y2 <- 1

if(x$univar.gamlss == TRUE) stop("This function is not suitable for univariate models.")
if(x$triv == TRUE ) stop("This function is currently not suitable for trivariate binary models.")
#if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && x$margins[2] %in% c(x$VC$m1d, x$VC$m2d)) stop("This function is currently not suitable for models involving discrete margins.")
#if(x$VC$Cont == "YES" && x$margins[1] %in% c(x$VC$m1d, x$VC$m2d)) stop("This function is currently not suitable for models involving a discrete margin.")

if(missing(y1)) stop("You must provide a value for y1.")
if(missing(y2)) stop("You must provide a value for y2.")

if(!(type %in% c("bivariate","independence"))) stop("Error in parameter type value. It should be one of: bivariate, independence.")
if(!(cond %in% c(0,1,2))) stop("Error in parameter cond value. It should be one of: 0, 1, 2.")

if( type %in% c("independence") && x$VC$gamlssfit == FALSE) stop("You need to re-fit the model and set gamlssfit = TRUE to obtain probabilities under independence.")


if(x$margins[1] %in% c(x$VC$m1d, x$VC$m2d) && (!is.wholenumber(y1) || y1 < 0)) stop("The value for y1 must be discrete and positive.")
if(x$margins[2] %in% c(x$VC$m1d, x$VC$m2d) && (!is.wholenumber(y2) || y2 < 0)) stop("The value for y2 must be discrete and positive.")


if( x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) && !(y1 %in% c(0,1)) ) stop("The value for y1 must be either 0 or 1.")

if( x$VC$Cont == "NO" && x$margins[2] %in% bin.link){ if( !(y1 %in% c(0,1)) || !(y2 %in% c(0,1))   ) stop("The value for y1 and/or y2 must be either 0 or 1.") }

if(!missing(newdata) && x$BivD %in% x$BivD2) stop("Prediction for models based on mixed copulae and a new dataset is not feasible.")

######################################################################################################
######################################################################################################

if(x$VC$Cont == "YES" && x$surv == FALSE )              rr <- jc.probs1(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link)
if(x$VC$Cont == "NO" && !(x$margins[2] %in% bin.link) ) rr <- jc.probs2(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link) 
if(x$VC$Cont == "NO" &&   x$margins[2] %in% bin.link)   rr <- jc.probs3(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link)    

if(x$VC$Cont == "YES" && x$surv == TRUE && x$margins[1] %in% bin.link ) rr <- jc.probs4(x, y1, y2, newdata, type, cond, intervals, n.sim, prob.lev, cont1par, cont2par, cont3par, bin.link)


# copulaReg
# SemiParBIV and copulaSampleSel
# SemiParBIV

######################################################################################################
######################################################################################################

p12s <- rr$p12s
p12  <- rr$p12
p1   <- rr$p1
p2   <- rr$p2

######################################################################################################
######################################################################################################

if(intervals == TRUE){ 

CIp12 <- rowQuantiles(p12s, probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE)

if(length(p12) > 1)  {res <- data.frame(p12, CIp12, p1, p2);       names(res)[2:3] <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2)))}
if(length(p12) == 1) {res <- data.frame(t(c(p12, CIp12, p1, p2))); names(res) <- c("p12",names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2))),"p1","p2")}

}else{

res <- data.frame(p12, p1, p2)

}


return(res)



}



