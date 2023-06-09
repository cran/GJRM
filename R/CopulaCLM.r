CopulaCLM <- function(formula, data = list(), weights = NULL, subset = NULL,
                      Model = "B", BivD = "N", margins = c("probit","N"),
                      dof = 3, gamlssfit = FALSE,
                      fp = FALSE, hess = TRUE, infl.fac = 1, theta.fx = NULL, 
                      rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                      gc.l = FALSE, parscale, extra.regI = "t", intf = FALSE, knots = NULL,
                      drop.unused.levels = TRUE, #ind.ord = FALSE,
                      min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999){
  

##########################################################################################
# Model set up and starting values
##########################################################################################

##### Model set up #####

i.rho <- sp <- qu.mag <- n.sel <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- inde <- y2m <- K1 <- NULL  
end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- l.sp9 <- i.rho <- 0
gam1 <- gam2 <- gam3 <- gam4 <- gam5 <- gam6 <- gam7 <- gam8 <- gam9 <- gamlss2 <- dof.st <- NULL
gamlss2 <- NULL
Sl.sf <- NULL
sp.method <- "perf"
sel2 <- sel2.p <- sel2.m <- sel2.mm <- sel2.pm <- c2 <- D21 <- D22 <- NA
  
sp1 <- sp2 <- c.gam2 <- X2s <- X3s <- NULL
sp3 <- gp3 <-   gam3 <- X3         <- NULL  
sp4 <- gp4 <-   gam4 <- X4         <- NULL  
sp5 <- gp5 <-   gam5 <- X5         <- NULL   
sp6 <- gp6 <-   gam6 <- X6         <- NULL  
sp7 <- gp7 <-   gam7 <- X7         <- NULL  
sp8 <- gp8 <-   gam8 <- X8         <- NULL  
log.sig2 <- log.nu <- NULL
    
BivD2 <- c("C0C90", "C0C270", "C180C90", "C180C270",
           "J0J90", "J0J270", "J180J90", "J180J270",
           "G0G90", "G0G270", "G180G90", "G180G270",
           "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270")  # 13:16
  
opc  <- c("N", "C0", "C90", "C180", "C270", 
               "J0", "J90", "J180", "J270",
               "G0", "G90", "G180", "G270", "F", "AMH", "FGM", "T", "PL", "HO","GAL0", "GAL90", "GAL180", "GAL270") # family for GAL is 62:65
               
scc  <- c("C0" , "C180","GAL0" , "GAL180", "J0" , "J180", "G0","G180",BivD2)

sccn <- c("C90", "C270", "GAL90", "GAL270","J90", "J270", "G90", "G270")

mb   <- c("B", "BSS", "BPO", "BPO0")
m2   <- c("N", "GU", "rGU", "LO", "LN", "WEI", "iG", "GA", "BE", "FISK","GP","GPII","GPo")
m3   <- c("DAGUM", "SM","TW")
m1d  <- c("PO", "ZTP","DGP0") 
m2d  <- c("NBI", "NBII", "PIG","DGP","DGPII") 
bl   <- c("probit", "logit", "cloglog") # , "cauchit")   
M    <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, BivD = BivD, 
             opc = opc, extra.regI = extra.regI, margins = margins, bl = bl, intf = intf,
             theta.fx = theta.fx, Model = Model, mb = mb, BivD2 = BivD2, dof = dof)  
surv.flex <- FALSE

ct  <- data.frame(c(opc), c(1:14,55,56,57,60,61,62:65))

cta <- data.frame(c(opc), c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2,60,61,62:65))
                   
if(BivD %in% BivD2){
  
if(BivD %in% BivD2[1 : 4 ]) BivDt <- "C0" 
if(BivD %in% BivD2[5 : 12]) BivDt <- "J0"
if(BivD %in% BivD2[13 :16]) BivDt <- "C0" # useful for ass dep function but we calculate it differently, so ok like this

  
nC  <-  ct[which( ct[, 1] == BivDt), 2]
nCa <- cta[which(cta[, 1] == BivDt), 2]     
  
}
  
if(!(BivD %in% BivD2)){
    
nC  <-  ct[which( ct[, 1] == BivD), 2]
nCa <- cta[which(cta[, 1] == BivD), 2]     
    
}

# nCa not important for GAL as we do not use VineCopula

################################################## 
  
if(!is.list(formula)) stop("You must specify a list of equations.")

M$l.flist <- l.flist <- length(formula)
pream.wm(formula, margins, M, l.flist, type = "ord") 
form.check(formula, l.flist)   
  

################################################## 


cl <- match.call()
mf <- match.call(expand.dots = FALSE)
            
pred.varR <- pred.var(formula, l.flist) 
   
v1     <- pred.varR$v1  
v2     <- pred.varR$v2
pred.n <- pred.varR$pred.n  
  
fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
environment(fake.formula) <- environment(formula[[1]])

mf$formula <- fake.formula 


mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$ordinal <- mf$knots <- mf$dof <- mf$intf <- mf$theta.fx <- mf$Model <- mf$BivD <- mf$margins <- mf$fp <- mf$hess <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
mf$drop.unused.levels <- drop.unused.levels
mf$ind.ord <- NULL
  
#if(Model=="BSS") mf$na.action <- na.pass
  
mf[[1]] <- as.name("model.frame")
data    <- eval(mf, parent.frame())
  
if(gc.l == TRUE) gc()  

#if(Model=="BSS"){     
#   
#data[is.na(data[, v1[1]]), v1[1]] <- 0
#indS <- data[, v1[1]]    
#indS[is.na(indS)] <- 0   
#indS <- as.logical(indS)  
#data[indS == FALSE, v2[1]] <- 0  
#data <- na.omit(data)   
#   
#}
                                      

if(!("(weights)" %in% names(data))){

weights <- rep(1,dim(data)[1]) 
data$weights <- weights
names(data)[length(names(data))] <- "(weights)"

} else weights <- data[,"(weights)"]    
  

formula.eq1 <- formula[[1]]
formula.eq2 <- formula[[2]] 
  
  
if(Model == "B"){

if(v1[1] %in% v2[-1]) end <- 1
if(v2[1] %in% v1[-1]) end <- 2

}


##### Equation 1 ##### 

gam1.false <- eval(substitute(gam(formula.eq1, gamma = infl.fac, weights = weights, 
                              data = data, knots = knots, fit = FALSE, drop.unused.levels = drop.unused.levels), list(weights = weights))) 
y1.false <- gam1.false$y
X1.false <- gam1.false$X
K1 <- resp.CLM(y1.false)

if(dim(X1.false)[2] == 1 & all(X1.false == 1)) stop("The equation for the ordinal renspose must include at least one regressor alongside the intercept.")


M$K1 <- K1 # K1 computed and added to M


gam1 <- eval(substitute(gam(formula.eq1, family = ocat(R = K1), gamma = infl.fac, weights = weights, 
                            data = data, knots = knots, drop.unused.levels = drop.unused.levels), list(weights = weights))) 


X1    <- model.matrix(gam1); if(all(X1[, 1] == 1)) X1 <- as.matrix(X1[, -1]) 
X1.d2 <- dim(X1)[2]

l.sp1 <- length(gam1$sp)
if(l.sp1 != 0) sp1 <- gam1$sp
  
y1 <- gam1$y
n <- length(y1) 
  
gp1 <- K1 + gam1$nsdf - 2 # cut-points are added (K1 - 1) and the intercept removed (- 1)
    
inde <- rep(TRUE, n)


##### Equation 2: B and continuous response #####

if(Model=="B" && !(margins[2] %in% bl)){
  
is_ordcon <- TRUE
is_ordord <- FALSE

form.eq12R   <- form.eq12(formula.eq2, data, v2, margins[2], m1d, m2d)   
formula.eq2  <- form.eq12R$formula.eq1
formula.eq2r <- form.eq12R$formula.eq1r
y2           <- form.eq12R$y1
y2.test      <- form.eq12R$y1.test 
y2m          <- form.eq12R$y1m  
  
gam2         <- eval(substitute(gam(formula.eq2, gamma = infl.fac, weights = weights, data = data, knots = knots, drop.unused.levels = drop.unused.levels), list(weights = weights)))
gam2$formula <- formula.eq2r  
names(gam2$model)[1] <- as.character(formula.eq2r[2])

y2 <- y2.test  
if(margins[2] %in% c("LN")) y2 <- log(y2) 

X2 <- model.matrix(gam2)
X2.d2 <- dim(X2)[2]
l.sp2 <- length(gam2$sp) ; if(l.sp2 != 0) sp2 <- gam2$sp
    
#cy <- 1 - y1
    
gp2 <- gam2$nsdf

K2 <- NULL # Added for distinguishing between the ordinal-continuous model from the ordinal-ordinal one

} 


##### Equation 2: B and ordinal response #####

if(Model=="B" && margins[2] %in% bl){

is_ordord <- TRUE # This is needed to distinguis between a mixed model from an ordinal-ordinal one
is_ordcon <- FALSE

gam2.false <- eval(substitute(gam(formula.eq2, gamma = infl.fac, weights = weights, 
                              data = data, knots = knots, fit = FALSE, drop.unused.levels = drop.unused.levels), list(weights = weights))) 
y2.false <- gam2.false$y
X2.false <- gam2.false$X
K2 <- resp.CLM(y2.false)

if(dim(X2.false)[2] == 1 & all(X2.false == 1)) stop("The second equation for the ordinal renspose must include at least one regressor alongside the intercept.")


M$K2 <- K2 # K2 computed and added to M


gam2 <- eval(substitute(gam(formula.eq2, family = ocat(R = K2), gamma = infl.fac, weights = weights, 
                            data = data, knots = knots, drop.unused.levels = drop.unused.levels), list(weights = weights))) 


X2    <- model.matrix(gam2); if(all(X2[, 1] == 1)) X2 <- as.matrix(X2[, -1]) 
X2.d2 <- dim(X2)[2]

l.sp2 <- length(gam2$sp)
if(l.sp2 != 0) sp2 <- gam2$sp
  
y2 <- gam2$y

# In the ordinal-ordinal model, the parameter vector takes the form (c1, c2, other parameters)'. As such, gp1 will account for the cut points for both equations 1 and 2.

gp1 <- K1 + K2 + gam1$nsdf - 3 # cut-points are added (K1 - 1 + K2 - 1) and the intercept removed (- 1)
gp2 <- gam2$nsdf - 1           # intercept removed
   
} 


###

M$is_ordcon <- is_ordcon
M$is_ordord <- is_ordord

# Tests below prevent that CopulaCLM() is used to estimate a model where one (or two) binary responses are modelled
test_bin1 <- K1 == 2
test2     <- ifelse(is_ordord, K2 == 2, FALSE)
if(test_bin1) stop("The levels of the first response variable must be greater than 2.")
if(test2)     stop("The levels of the second response variable must be greater than 2.")

###

##### Starting values for cut points ##### 

# N.B. the intercept is removed from the coefficient's vector

if(names(gam1$coefficients)[1] == "(Intercept)"){

gam1.int <- gam1$coefficients[1]
gam1$coefficients <- gam1$coefficients[-1]

}

c1 <- gam1$family$getTheta(TRUE) - gam1.int #; c1.ind <- c1

c1.ti <- rep(0, K1 - 1)
c1.ti[1] <- c1[1] ; for(i in 2 : (K1 - 1)) {c1.ti[i] <- sqrt(c1[i] - c1[i - 1])}


n.num_1 <- seq(1, K1 - 1) 
names(c1.ti) <- paste(paste("c1", n.num_1, sep = ""), "star", sep = ".")


if (is_ordord) {
	if(names(gam2$coefficients)[1] == "(Intercept)") {
		gam2.int <- gam2$coefficients[1]
		gam2$coefficients <- gam2$coefficients[-1]
	}

	c2 <- gam2$family$getTheta(TRUE) - gam2.int #; c2.ind <- c2

	c2.ti <- rep(0, K2 - 1)
	c2.ti[1] <- c2[1] ; for(i in 2 : (K2 - 1)) {c2.ti[i] <- sqrt(c2[i] - c2[i - 1])}

	n.num_2 <- seq(1, K2 - 1) 
	names(c2.ti) <- paste(paste("c2", n.num_2, sep = ""), "star", sep = ".")
}# else {
#	c2.ind <- NULL
#}


##### Starting values for dependence parameter #####

if(is.null(theta.fx)){

    
if(Model == "B"){ 

res1 <- residuals(gam1)
res2 <- residuals(gam2)

ass.s <- cor(res1, res2, method = "kendall")
ass.s <- sign(ass.s) * ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))  
  
}  

                          
i.rho <- ass.dp(ass.s, BivD, scc, sccn, nCa)



}

names(i.rho) <- "theta.star"  


##### Starting values for whole parameter vector ##### 

if(margins[1] %in% bl && is_ordcon){ # Mixed ordinal-continuous case

start.snR <- startsn(margins[2], y2)    
log.sig2  <- start.snR$log.sig2.1; names(log.sig2) <- "sigma.star"

#if(margins[2] %in% c(m3    )){ log.nu <- start.snR$log.nu.1; names(log.nu) <- "nu.star"}  

if(margins[2] %in% c(m2))  start.v <- c(c1.ti, gam1$coefficients, gam2$coefficients, log.sig2, i.rho)    

#if(margins[2] %in%   m3     )  start.v <- c(c1.ti, gam1$coefficients, gam2$coefficients, log.sig2, log.nu, i.rho)                                  

} 


if (margins[1] %in% bl && is_ordord) { # Ordinal-ordinal case
	start.v <- c(c1.ti, c2.ti, gam1$coefficients, gam2$coefficients, i.rho)    
} 

################################################## 
  
if(l.flist > 2){  
 
vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2 = log.sig2, log.nu = log.nu, n = n, drop.unused.levels = drop.unused.levels)  
  
overall.svGR <- overall.svG(formula, data, ngc = 2, margins, M, vo, gam1, gam2, type = "biv", inde = inde, c.gam2 = c.gam2, knots = knots)
    
if (is_ordcon) { start.v = c(c1.ti, overall.svGR$start.v) } # cut points added
if (is_ordord) { start.v = c(c1.ti, c2.ti, overall.svGR$start.v) } # cut points added
X3 = overall.svGR$X3; X4 = overall.svGR$X4; X5 = overall.svGR$X5
X6 = overall.svGR$X6; X7 = overall.svGR$X7; X8 = overall.svGR$X8
X3.d2 = overall.svGR$X3.d2; X4.d2 = overall.svGR$X4.d2; X5.d2 = overall.svGR$X5.d2
X6.d2 = overall.svGR$X6.d2; X7.d2 = overall.svGR$X7.d2; X8.d2 = overall.svGR$X8.d2
gp3 = overall.svGR$gp3; gp4 = overall.svGR$gp4; gp5 = overall.svGR$gp5
gp6 = overall.svGR$gp6; gp7 = overall.svGR$gp7; gp8 = overall.svGR$gp8
gam3 = overall.svGR$gam3; gam4 = overall.svGR$gam4; gam5 = overall.svGR$gam5
gam6 = overall.svGR$gam6; gam7 = overall.svGR$gam7; gam8 = overall.svGR$gam8
l.sp3 = overall.svGR$l.sp3; l.sp4 = overall.svGR$l.sp4; l.sp5 = overall.svGR$l.sp5
l.sp6 = overall.svGR$l.sp6; l.sp7 = overall.svGR$l.sp7; l.sp8 = overall.svGR$l.sp8
sp3 = overall.svGR$sp3; sp4 = overall.svGR$sp4; sp5 = overall.svGR$sp5
sp6 = overall.svGR$sp6; sp7 = overall.svGR$sp7; sp8 = overall.svGR$sp8
X3s = overall.svGR$X3s; X4s = overall.svGR$X4s

}  

################################################## 

GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8, gam9 = NULL)   


if((l.sp1 != 0 || l.sp2 != 0 || l.sp3 != 0 || l.sp4 != 0 || 
    l.sp5 != 0 || l.sp6 != 0 || l.sp7 != 0 || l.sp8 != 0   ) && fp == FALSE){ 

L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), l.gam4 = length(gam4$coefficients),
              l.gam5 = length(gam5$coefficients), l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), l.gam8 = length(gam8$coefficients),
              l.gam9 = 0)

L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, l.sp9 = l.sp9 )

sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
if (is_ordcon) { qu.mag <- S.m(GAM, L.SP, L.GAM, K1 = K1) }
if (is_ordord) { qu.mag <- S.m(GAM, L.SP, L.GAM, K1 = K1, K2 = K2) }




}

##################################################

if(missing(parscale)) parscale <- 1   

respvec <- list(y1      = y1     ,
                y2      = y2     ,
                y1.y2   = y1.y2  , 
                y1.cy2  = y1.cy2 , 
                cy1.y2  = cy1.y2 , 
                cy1.cy2 = cy1.cy2, 
                cy1     = cy1    ,
                cy      = cy     , univ = 0)

### 1st ordinal equation

sel <- model.matrix(~ as.factor(y1) - 1)

sel.p <- as.matrix(sel[, 3 : K1])
sel.m <- as.matrix(sel[, 2 : (K1 - 1)])

sel.p[, (K1 - 2)] <- as.matrix(sel[, K1])
sel.m[, (K1 - 2)] <- as.matrix(sel[, (K1 - 1)])

sel.mm <- sel.pm <- matrix(nrow = n, ncol = K1 - 3, 0)

if( K1 > 3 ){
   for (l in 1:(K1 - 3)) {sel.pm [, l] <- rowSums(sel.p[, l : (K1 - 2)])} ; sel.p[, 1 : (K1 - 3)] <- sel.pm
   for (l in 1:(K1 - 3)) {sel.mm [, l] <- rowSums(sel.m[, l : (K1 - 2)])} ; sel.m[, 1 : (K1 - 3)] <- sel.mm
            }

c1 <- rep(0, K1 - 1)

D11 <- rowSums(sel[, 1 : (K1 - 1)])
D12 <- rowSums(sel[, 2 : K1])


### 2nd ordinal equation 

if (is_ordord) {
	sel2 <- model.matrix(~ as.factor(y2) - 1)

	sel2.p <- as.matrix(sel2[, 3 : K2])
	sel2.m <- as.matrix(sel2[, 2 : (K2 - 1)])

	sel2.p[, (K2 - 2)] <- as.matrix(sel2[, K2])
	sel2.m[, (K2 - 2)] <- as.matrix(sel2[, (K2 - 1)])

	sel2.mm <- sel2.pm <- matrix(nrow = n, ncol = K2 - 3, 0)

	if ( K2 > 3 ) {
   		for (l in 1:(K2 - 3)) {sel2.pm [, l] <- rowSums(sel2.p[, l : (K2 - 2)])} ; sel2.p[, 1 : (K2 - 3)] <- sel2.pm
   		for (l in 1:(K2 - 3)) {sel2.mm [, l] <- rowSums(sel2.m[, l : (K2 - 2)])} ; sel2.m[, 1 : (K2 - 3)] <- sel2.mm
	}

	c2 <- rep(0, K2 - 1)

	D21 <- rowSums(sel2[, 1 : (K2 - 1)])
	D22 <- rowSums(sel2[, 2 : K2])
}

### 

my.env <- new.env()
my.env$signind <- 1

lsgam1 <- length(gam1$smooth)
lsgam2 <- length(gam2$smooth)
lsgam3 <- length(gam3$smooth)
lsgam4 <- length(gam4$smooth)
lsgam5 <- length(gam5$smooth)
lsgam6 <- length(gam6$smooth)
lsgam7 <- length(gam7$smooth)
lsgam8 <- length(gam8$smooth)
lsgam9 <- length(gam9$smooth)


VC <- list(lsgam1 = lsgam1, robust = FALSE,  
           lsgam2 = lsgam2, Sl.sf = Sl.sf, sp.method = sp.method,
           lsgam3 = lsgam3, 
           lsgam4 = lsgam4,
           lsgam5 = lsgam5,
           lsgam6 = lsgam6,
           lsgam7 = lsgam7, 
           lsgam8 = lsgam8, lsgam9 = lsgam9, 
           sel1 = sel , sel1.p = sel.p , sel1.m = sel.m , sel1.mm = sel.mm , sel1.pm = sel.pm , c1 = c1, D11 = D11, D12 = D12, # added for CopulaCLM
           sel2 = sel2, sel2.p = sel2.p, sel2.m = sel2.m, sel2.mm = sel2.mm, sel2.pm = sel2.pm, c2 = c2, D21 = D21, D22 = D22, # added for CopulaCLM
           K1 = K1, # added for CopulaCLM
           K2 = K2, # added for CopulaCLM
           X1 = X1, inde = inde, my.env = my.env,
           X2 = X2, 
           X3 = X3,
           X4 = X4, 
           X5 = X5, 
           X6 = X6, 
           X7 = X7,
           X8 = X8,
           X1.d2 = X1.d2,
           X2.d2 = X2.d2,
           X3.d2 = X3.d2,
           X4.d2 = X4.d2,
           X5.d2 = X5.d2,
           X6.d2 = X6.d2,
           X7.d2 = X7.d2,
           X8.d2 = X8.d2,
           gp1 = gp1, 
           gp2 = gp2,
           gp3 = gp3,
           gp4 = gp4, 
           gp5 = gp5,
           gp6 = gp6,  
           gp7 = gp7, 
           gp8 = gp8, 
           l.sp1 = l.sp1, 
           l.sp2 = l.sp2,
           l.sp3 = l.sp3,
           l.sp4 = l.sp4,
           l.sp5 = l.sp5,
           l.sp6 = l.sp6,    
           l.sp7 = l.sp7,
           l.sp8 = l.sp8, 
           l.sp9 = 0,
           infl.fac = infl.fac,
           weights = weights,
           fp = fp, univ.gamls = FALSE,
           hess = hess, nCa = nCa,
           Model = Model, gamlssfit = gamlssfit,
           end = end,
           BivD = BivD, dof.st = log(dof - 2), dof = dof, 
           nC = nC, gc.l = gc.l, n = n, extra.regI = extra.regI,
           parscale = parscale, margins = margins,
           Cont = "NO", ccss = "no", m2 = m2, m3 = m3, m2d = m2d, m1d = m1d, m3d = NULL, bl = bl,
           X2s = X2s, X3s = X3s, triv = FALSE, y2m = y2m,
           theta.fx = theta.fx, i.rho = i.rho, 
           BivD2 = BivD2, cta = cta, ct = ct, zerov = -10, surv.flex = surv.flex, gp2.inf = NULL,
           informative = "no", sp.fixed = NULL,
           zero.tol = 1e-02,
           min.dn = min.dn, min.pr = min.pr, max.pr = max.pr,
           is_ordcon = is_ordcon, 
           is_ordord = is_ordord) # original n only needed in SemiParBIV.fit
           
if(gc.l == TRUE) gc()           
             
################################################## 

#if(gamlssfit == TRUE && !(margins[2] %in% bl)){ # This option is not available for the ordinal-ordinal model
#
#form.gamlR <- form.gaml(formula, l.flist, M, type = "biv")
#
#gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = subset,
#                 margin = margins[2], infl.fac = infl.fac,
#                 rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
#                 gc.l = gc.l, parscale = 1, extra.regI = extra.regI, drop.unused.levels = drop.unused.levels), list(weights = weights)))
#
#
## Updated starting values
#
#MM <- M; MM$BivD <- "N" # this is for T case, dof is never estimated...
#  
#SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
#gamls.upsvR <- gamls.upsv(gamlss1 = NULL, gamlss2, margins, MM, l.flist, nstv = NULL, VC, GAM, SP, type = "biv")
#sp <- gamls.upsvR$sp
#start.v <- c(c1.ti, gamls.upsvR$start.v) # cut points added
#
#}


##########################################################################################
# Model estimation
########################################################################################## 


func.opt <- func.OPT(margins, M, type = "biv")


############################## 

# Independence model #

# The part of the code below serves two purposes:
#    (i)  Estimating a bivariate model under independence (both OrdCon and OrdOrd); and
#    (ii) Producing better-calibrated start values bor the bivariate model without independence. This is achieved 
#         by fitting a gamlss for the OrdCon model, and using the estimated parameter vector under independence 
#         for the OrdOrd model.

if (gamlssfit == "TRUE") {
	VC.ind <- VC

	VC.ind$ind.ord <- TRUE
	VC.ind$BivD <- "J0" # This is just used to fool the fitting procedure: no matter which copula is used in the independence model
	VC.ind$nC  <- ct [which(ct [, 1] == VC.ind$BivD), 2]
	VC.ind$nCa <- cta[which(cta[, 1] == VC.ind$BivD), 2]

	if (is_ordcon) {
		if (is.null(X3)) {
			drop.ind <- K1 + X1.d2 + X2.d2 + 1
		} else {
			drop.ind <- (K1 + X1.d2 + X2.d2 + X3.d2) : (K1 + X1.d2 + X2.d2 + X3.d2 + X4.d2 - 1)
		}
	} else {
		if(is.null(X3)) {
			drop.ind <- K1 + K2 + X1.d2 + X2.d2 - 1
		} else {
			drop.ind <- (K1 + K2 + X1.d2 + X2.d2 - 1) : (K1 + K2 + X1.d2 + X2.d2 + X3.d2 - 2)
		}
	}

	VC.ind$drop.ind <- drop.ind
	start.v.ind <- start.v[-drop.ind]

	if (is_ordcon) { # The gamlss will be fitted only when the second marginal is continuous 
		form.gamlR <- form.gaml(formula, l.flist, M, type = "biv")

		gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = subset,
                 		margin = margins[2], infl.fac = infl.fac,
                 		rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                 		gc.l = gc.l, parscale = 1, extra.regI = extra.regI, drop.unused.levels = drop.unused.levels), list(weights = weights)))

		# Updated starting values

		MM <- M; MM$BivD <- "N" # this is for T case, dof is never estimated...
  
		SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)

		GAM.ind <- GAM
		GAM.ind$gam1$coefficients <- GAM.ind$gam4$coefficients <- NULL # When gamlssfit == "TRUE" only the continuous response is estimated under gamlss2
		
		gamls.upsvR <- gamls.upsv(gamlss1 = NULL, gamlss2, margins, MM, l.flist, nstv = NULL, VC, GAM.ind, SP, type = "biv")
		sp.ind <- gamls.upsvR$sp

		w.theta.star <- which(names(gamls.upsvR$start.v) == "theta.star")                        #
		if (length(w.theta.star) != 0) gamls.upsvR$start.v <- gamls.upsvR$start.v[-w.theta.star] # theta.star is removed fom the estimated parameters

		if (is.null(X3)) {
			start.v.ind[(K1 + X1.d2) : (K1 + X1.d2 + X2.d2)] <- gamls.upsvR$start.v # When X3 is null, sigma.star is estimated as a scalar
		} else {
			start.v.ind[(K1 + X1.d2) : (K1 + X1.d2 + X2.d2 + X3.d2 - 1)] <- gamls.upsvR$start.v
		}
	}

	qu.mag.ind <- qu.mag
	w.off.ind <- which(qu.mag.ind$off > length(start.v.ind))

	if (length(w.off.ind) > 0) {
		qu.mag.ind$rank <- qu.mag.ind$rank[-w.off.ind]
		qu.mag.ind$off  <- qu.mag.ind$off [-w.off.ind]
		qu.mag.ind$Ss   <- qu.mag.ind$Ss  [-w.off.ind]
		
		sp.ind <- sp[-w.off.ind] # QUESTION: does this create an issue given that I define sp.ind also above when is_ordcon?

		l.spvec <- c(VC.ind$l.sp1, VC.ind$l.sp2, VC.ind$l.sp3, VC.ind$l.sp4, VC.ind$l.sp5, VC.ind$l.sp6, VC.ind$l.sp7, VC.ind$l.sp8)
			l.spvec[length(formula) : length(l.spvec)] <- rep(0) #l.spvec[w.off : length(l.spvec)] <- rep(0) # The equation for theta is the last input in formula.
		VC.ind$l.sp1 <- l.spvec[1]
		VC.ind$l.sp2 <- l.spvec[2]
		VC.ind$l.sp3 <- l.spvec[3]
		VC.ind$l.sp4 <- l.spvec[4]
		VC.ind$l.sp5 <- l.spvec[5]
		VC.ind$l.sp6 <- l.spvec[6]
		VC.ind$l.sp7 <- l.spvec[7]
		VC.ind$l.sp8 <- l.spvec[8]
	} else {
		sp.ind <- sp # NULL
	}

	# The independence model is fitted

	fit_ind <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v.ind,
                              	  rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                                  respvec = respvec, VC = VC.ind, sp = sp.ind, qu.mag = qu.mag.ind)

	coeff.ind <- fit_ind$fit$argument

	if (is_ordcon) { # The start values of the bivariate model are updated
		start.v[1 : (min(drop.ind) - 1)] <- coeff.ind
	} else {
		start.v[1 : (K1 + K2 + X1.d2 + X2.d2 - 2)] <- coeff.ind
	}

	# Cut points are transformed back and added to the parameter vector

	infty <- 1e+25

	c1.ti <- coeff.ind[1 : (K1 - 1)]
	c1 <- rep(0, K1 - 1) ; c1[1] <- c1.ti[1] ; for (i in 2 : (K1 - 1)) c1[i] <- c1[i - 1] + c1.ti[i]^2
	c1 <- sign(c1) * pmin(10000 * infty, abs(c1))

	coeff.ind[1 : (K1 - 1)] <- c1
	names(coeff.ind)[1 : (K1 - 1)] <- paste("c1", n.num_1, sep = "")

	if (is_ordord) {
		c2.ti <- coeff.ind[K1 : (K1 + K2 - 2)]
		c2 <- rep(0, K2 - 1) ; c2[1] <- c2.ti[1] ; for (i in 2 : (K2 - 1)) c2[i] <- c2[i - 1] + c2.ti[i]^2
		c2 <- sign(c2) * pmin(10000 * infty, abs(c2))

		coeff.ind[K1 : (K1 + K2 - 2)] <- c2
		names(coeff.ind)[K1 : (K1 + K2 - 2)] <- paste("c2", n.num_2, sep = "")
	}

	# Post-estimation

	CopulaCLMFit.p.ind <- SemiParBIV.fit.post(SemiParFit = fit_ind, Model = Model, VC = VC.ind, GAM)     

	# The coefficients of the independence model are stored in a more user-friendly way

	if (is_ordcon) { 
		coefficients.ind <- list(c1 = coeff.ind[1 : (K1 - 1)],
                                         beta1 = coeff.ind[K1 : (K1 + X1.d2 - 1)],
                                         beta2 = coeff.ind[(K1 + X1.d2) : (K1 + X1.d2 + X2.d2 - 1)])
		if (!is.null(X3)){
			coefficients.ind$sigma2 <- coeff.ind[(K1 + X1.d2 + X2.d2) : (K1 + X1.d2 + X2.d2 + X3.d2 - 1)]
		}
	} else {
		coefficients.ind <- list(c1 = coeff.ind[1 : (K1 - 1)],
                                         c2 = coeff.ind[K1 : (K1 + K2 - 2)],
                                         beta1 = coeff.ind[(K1 + K2 - 1) : (K1 + K2 + X1.d2 - 2)],
                                         beta2 = coeff.ind[(K1 + K2 + X1.d2 - 1) : (K1 + K2 + X1.d2 + X2.d2 - 2)])
	}

	# When the independence model is estimated, the bivariate model is subsequently estimated hence I should set ind.ord to FALSE.

	VC$ind.ord <- FALSE
	Vb.ind <- CopulaCLMFit.p.ind$Vb
} else {
	VC$ind.ord <- FALSE
	coefficients.ind <- NULL
	Vb.ind <- NULL
}

########################################################################################## 

CopulaCLMFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v,
                               rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                               respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag)


# Cut points are transformed back and added to the parameter vector

infty <- 1e+25

c1.ti <- CopulaCLMFit$fit$argument[1 : (K1 - 1)]
c1 <- rep(0, K1 - 1) ; c1[1] <- c1.ti[1] ; for (i in 2 : (K1 - 1)) c1[i] <- c1[i - 1] + c1.ti[i]^2
c1 <- sign(c1) * pmin(10000 * infty, abs(c1))

CopulaCLMFit$fit$argument[1 : (K1 - 1)] <- c1
names(CopulaCLMFit$fit$argument)[1 : (K1 - 1)] <- paste("c1", n.num_1, sep = "")

if (is_ordord) {
	c2.ti <- CopulaCLMFit$fit$argument[K1 : (K1 + K2 - 2)]
	c2 <- rep(0, K2 - 1) ; c2[1] <- c2.ti[1] ; for (i in 2 : (K2 - 1)) c2[i] <- c2[i - 1] + c2.ti[i]^2
	c2 <- sign(c2) * pmin(10000 * infty, abs(c2))

	CopulaCLMFit$fit$argument[K1 : (K1 + K2 - 2)] <- c2
	names(CopulaCLMFit$fit$argument)[K1 : (K1 + K2 - 2)] <- paste("c2", n.num_2, sep = "")

} 


########################################################################################## 
# Post estimation
##########################################################################################

CopulaCLMFit.p <- SemiParBIV.fit.post(SemiParFit = CopulaCLMFit, Model = Model, VC = VC, GAM)                                      
CopulaCLMFit   <- CopulaCLMFit.p$SemiParFit # useful for SS models, eta2 calculatons etc.
 
y2.m <- y2
if(margins[2] == "LN") y2.m <- exp(y2)


################################################## 

if(gc.l == TRUE) gc()


##################################################

cov.c(CopulaCLMFit)


##################################################

gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- 
gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 


##################################################

# Bit below useful for AT calculations when end is continuous

if(!(Model == "B" && !(margins[2] %in% bl) && end == 2)) {dataset <- NULL; rm(data)} else {attr(data, "terms") <- NULL; dataset <- data; rm(data)} 

L <- list(fit = CopulaCLMFit$fit, dataset = dataset, formula = formula, CopulaCLMFit = CopulaCLMFit, robust = FALSE,
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6,
          gam7 = gam7, gam8 = gam8,
          coefficients = CopulaCLMFit$fit$argument, coef.t = NULL, iterlimsp = iterlimsp,
          weights = weights,
          sp = CopulaCLMFit.p$sp, iter.sp = CopulaCLMFit$iter.sp,
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3,
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6,
          l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,l.sp9 = l.sp9, gam9 = gam9,
          fp = fp,
          iter.if = CopulaCLMFit$iter.if, iter.inner = CopulaCLMFit$iter.inner,
          theta = CopulaCLMFit.p$theta,
          theta.a = CopulaCLMFit.p$theta.a,
          OR = CopulaCLMFit.p$OR,
          GM = CopulaCLMFit.p$GM,
          n = n, n.sel = n.sel,
          #c1.ind = c1.ind, c2.ind = c2.ind,
          K1 = K1, K2 = K2,
          X1 = X1, X2 = X2, X3 = X3, X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2,
          X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X4.d2 = X4.d2, X5.d2 = X5.d2,
          X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,
          He = CopulaCLMFit.p$He, HeSh = CopulaCLMFit.p$HeSh, Vb = CopulaCLMFit.p$Vb, Ve = CopulaCLMFit.p$Ve,
          F = CopulaCLMFit.p$F, F1 = CopulaCLMFit.p$F1,
          t.edf = CopulaCLMFit.p$t.edf, edf = CopulaCLMFit.p$edf,
          edf11 = CopulaCLMFit.p$edf11,
          edf1 = CopulaCLMFit.p$edf1, edf2 = CopulaCLMFit.p$edf2, edf3 = CopulaCLMFit.p$edf3,
          edf4 = CopulaCLMFit.p$edf4, edf5 = CopulaCLMFit.p$edf5, edf6 = CopulaCLMFit.p$edf6,
          edf7 = CopulaCLMFit.p$edf7, edf8 = CopulaCLMFit.p$edf8,
          edf1.1 = CopulaCLMFit.p$edf1.1, edf1.2 = CopulaCLMFit.p$edf1.2, edf1.3 = CopulaCLMFit.p$edf1.3,
          edf1.4 = CopulaCLMFit.p$edf1.4, edf1.5 = CopulaCLMFit.p$edf1.5, edf1.6 = CopulaCLMFit.p$edf1.6,
          edf1.7 = CopulaCLMFit.p$edf1.7, edf1.8 = CopulaCLMFit.p$edf1.8,
          R = CopulaCLMFit.p$R,
          bs.mgfit = CopulaCLMFit$bs.mgfit, conv.sp = CopulaCLMFit$conv.sp,
          wor.c = CopulaCLMFit$wor.c,
          p11 = CopulaCLMFit$fit$p11, p10 = CopulaCLMFit$fit$p10, p01 = CopulaCLMFit$fit$p01, p00 = CopulaCLMFit$fit$p00,
          p1 = CopulaCLMFit$fit$p1, p2 = CopulaCLMFit$fit$p2,
          eta1 = CopulaCLMFit$fit$eta1, eta2 = CopulaCLMFit$fit$eta2, etad = CopulaCLMFit$fit$etad,
          etas = CopulaCLMFit$fit$etas, etan = CopulaCLMFit$fit$etan,
          y1 = y1, y2 = y2.m,
          BivD = BivD, margins = margins,
          logLik = CopulaCLMFit.p$logLik,
          nC = nC, hess = hess,
          respvec = respvec, inde = inde,
          qu.mag = qu.mag, sigma = CopulaCLMFit.p$sigma, sigma.a = CopulaCLMFit.p$sigma.a,
          sigma2 = CopulaCLMFit.p$sigma2, sigma2.a = CopulaCLMFit.p$sigma2.a,
          nu = CopulaCLMFit.p$nu, nu.a = CopulaCLMFit.p$nu.a, Vb.t = CopulaCLMFit.p$Vb.t,
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          X2s = X2s, X3s = X3s, p1n = CopulaCLMFit.p$p1n , p2n = CopulaCLMFit.p$p2n, 
          VC = VC, Model = Model, magpp = CopulaCLMFit$magpp,
          gamlssfit = gamlssfit, Cont = "NO", tau = CopulaCLMFit.p$tau, 
          tau.a = CopulaCLMFit.p$tau.a, l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE,
          gamlss = gamlss2, BivD2 = BivD2, dof = dof, dof.a = dof, call = cl,
          surv = FALSE, surv.flex = surv.flex, ordinal = TRUE, drop.unused.levels = drop.unused.levels, Model = Model,
          is_ordcon = is_ordcon,
          is_ordord = is_ordord,
          coefficients.ind = coefficients.ind, Vb.ind = Vb.ind)


if(BivD %in% BivD2){

L$teta1     <- CopulaCLMFit$fit$teta1
L$teta.ind1 <- CopulaCLMFit$fit$teta.ind1
L$teta2     <- CopulaCLMFit$fit$teta2
L$teta.ind2 <- CopulaCLMFit$fit$teta.ind2
L$Cop1      <- CopulaCLMFit$fit$Cop1
L$Cop2      <- CopulaCLMFit$fit$Cop2

}


class(L) <- c("CopulaCLM", "SemiParBIV", "gjrm")

L

}
