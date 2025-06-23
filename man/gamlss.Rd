\name{gamlss}
\alias{gamlss}

\title{Generalised Additive Models for Location, Scale and Shape and Beyond}

\description{ 
\code{gamlss} fits flexible univariate regression models for several continuous and discrete distributions as well as survival outcomes, and types of covariate 
effects. When first designed, the purpose of this function was only to provide, in some cases, starting values 
for the simultaneous models in the package. At a later stage, it 
was made available in the form of a proper function should the user wish to fit 
univariate models using the general estimation approach of this package. The continuous and discrete distributions used here 
are parametrised according to Rigby and Stasinopoulos (2005).
 
}



\usage{
gamlss(formula, data = list(), weights = NULL, subset = NULL, offset = NULL, 
       family = "N", cens = NULL, type.cens = "R", ub.t = NULL, left.trunc = 0,
       robust = FALSE, rc = 3, lB = NULL, uB = NULL, infl.fac = 1, 
       rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
       gc.l = FALSE, parscale, gev.par = -0.25,
       chunk.size = 10000, knots = NULL,
       informative = "no", inform.cov = NULL, family2 = "-cloglog", 
       fp = FALSE, sp = NULL,
       drop.unused.levels = TRUE, siginit = NULL, shinit = NULL,
       sp.method = "perf", hrate = NULL, d.lchrate = NULL, d.rchrate = NULL,
       d.lchrate.td = NULL, d.rchrate.td = NULL, truncation.time = NULL,
       min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.9999999, ygrid.tol = 1e-08)
}



\arguments{ 

\item{formula}{List of equations. This should contain one or more equations.} 

\item{data}{A data frame.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{offset}{Optional vector specifying an offset for use in fitting. Option introduced for dealing with offset with discrete distributions.}

\item{family}{Possible choices are Gaussian ("N"), truncated Gaussian ("tN"), Tweedie ("TW"),
               log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), generelised Pareto ("GP"), 
               generelised Pareto II ("GPII") where the shape parameter is forced to be > -0.5, 
               generelised Pareto (with orthogonal parametrisation) ("GPo") where the shape parameter is forced to be > -0.5,
               discrete generelised Pareto ("DGP"),
               discrete generelised Pareto II ("DGPII") where the shape parameter is forced to be positive, discrete generelised Pareto derived
               under the scenario in which shape = 0 ("DGP0"), logistic ("LO"), Weibull ("WEI"), Inverse Gaussian ("IG"), gamma ("GA"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic), Poisson ("P"), truncated 
               Poisson ("tP"), negative binomial - type I ("NBI"), negative 
               binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG"), truncated negative binomial - type I ("tNBI"), truncated negative 
               binomial - type II ("tNBII"), truncated Poisson inverse Gaussian ("tPIG"), generalised extreme value link function ("GEVlink", this 
               is used for binary responses and is more stable and faster than the \code{R} package \code{bgeva}). For survival models, \code{family} can 
               be "-cloglog" (similar to generalised proportional hazards), "-logit" (similar to generalised proportional odds), "-probit" (generalised probit).}

\item{cens}{This is required for a survival model. When \code{type.cens} is different from \code{mixed}, this variable can be equal to 1 if the event occurred 
            and 0 otherwise. If \code{type.cens = "mixed"} then \code{cens} is a mixed factor variable (made up of four possible 
            categories: \code{I} for interval, \code{L} for left, \code{R} for right, and \code{U} for uncensored.}
 
\item{type.cens}{Type of censoring mechanism. This can be "R", "L", "I" or "mixed".}
 
\item{ub.t}{Variable name of right/upper bound when \code{type.cens = "I"} or \code{type.cens = "mixed"} and interval censoring is present.}
 
 \item{left.trunc}{Value of truncation at left. Currently done for count distributions only.}
 
\item{robust}{If \code{TRUE} then the robust version of the model is fitted.}
             
\item{rc}{Robust constant.}       

\item{lB, uB}{Bounds for integral in robust case.}  

\item{infl.fac}{Inflation factor for the model degrees of freedom in the approximate AIC. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds.}

\item{rmax}{Maximum allowed trust region radius. This may be set very large. If set small, the algorithm traces a steepest 
            descent path.} 

\item{iterlimsp}{A positive integer specifying the maximum number of loops to be performed before the smoothing parameter estimation 
                 step is terminated.} 

\item{tolsp}{Tolerance to use in judging convergence of the algorithm when automatic smoothing parameter estimation is used.} 

\item{gc.l}{This is relevant when working with big datasets. If \code{TRUE} then the garbage collector is called more often than it is 
usually done. This keeps the memory footprint down but it will slow down the routine.}
      
\item{parscale}{The algorithm will operate as if optimizing objfun(x / parscale, ...) where parscale is a scalar. If missing then no 
rescaling is done. See the 
documentation of \code{trust} for more details.}    

\item{gev.par}{GEV link parameter.}
      
\item{chunk.size}{This is used for discrete robust models.}  

\item{knots}{Optional list containing user specified knot values to be used for basis construction.}
  
\item{informative}{If "yes" then informative censoring is assumed when using a survival model.}
  
\item{inform.cov}{If above is "yes" then a set of informative covariates must be provided.}
  
\item{family2}{In the informative survival case, the family for the censored equation can be different from that of the survival equation. Choices are
               "-cloglog" (siilar to generalised proportional hazards), "-logit" (similar to generalised proportional odds), "-probit" (generalised probit).}
  
\item{fp}{If \code{TRUE} then a fully parametric model with unpenalised regression splines if fitted.}              
  
\item{sp}{A vector of smoothing parameters can be provided here. Smoothing parameters must be supplied in the order that the smooth 
          terms appear in the model equation(s).}  
 
\item{drop.unused.levels}{By default unused levels are dropped from factors before fitting. For some smooths involving factor variables 
                          this may have to be turned off (only use if you know what you are doing).}  

\item{siginit, shinit}{For the GP and DGP distributions, initial values for sigma and shape may be provided.}

\item{sp.method}{Multiple smoothing automatic parameter selection is perf. efs is an alternative and only sensible option
                 for robust models.}
                 
\item{hrate}{Vector of population hazard rates computed at time of death of each uncensored patient. The length of \code{hrate} should be equal to the number of uncensored observations in the dataset. Needed in the context of excess hazard modelling when uncensored observations are present. Note that this includes left truncated uncensored observations as well.}          

\item{d.lchrate}{Vector of differences of population cumulative excess hazards computed at the age of the patient when the left 
                 censoring occurred and at the initial age of the patient. The length of \code{d.lchrate} should be equal to the number 
                 of left and/or interval censored observations in the dataset. Needed in the context of excess hazard modelling 
                 if left censored and/or interval censored observations are present. In the latter case, \code{d.rchrate} also need be provided.}

\item{d.rchrate}{Vector of differences of population cumulative excess hazards computed at the age of the patient when the at the right 
                 interval censoring time and at the initial age of the patient. The length of \code{d.rchrate} should be equal to the number 
                 of right censored and/or interval censored observations in the dataset. Needed in the context of excess hazard modelling 
                 if right censored and/or interval censored observations are present. In the latter case, \code{d.lchrate} also need be provided.}

\item{d.lchrate.td}{Vector of differences of population cumulative excess hazards computed at the age of the patient when the left 
                   censoring occurred and at the age of the patient when the truncation occurred. The length of \code{d.lchrate.td} should be 
                   equal to the number of left truncated left censored and/or left truncated interval censored observations in 
                   the dataset. Needed in the context of excess hazard modelling if left truncated left censored and/or left truncated 
                   interval censored observations are present. In the latter case, \code{d.rchrate.td} also need be provided.}

\item{d.rchrate.td}{Vector of differences of population cumulative excess hazards computed at the age of the patient when the right 
                   censoring occurred and at the age of the patient when the truncation occurred. The length of \code{d.rchrate.td} should be 
                   equal to the number of left truncated right censored and/or left truncated interval censored observations in the dataset. Needed in 
                   the context of excess hazard modelling if left truncated right censored and/or left truncated interval censored 
                   observations are present. In the latter case, \code{d.lchrate.td} also need be provided.}

\item{truncation.time}{Variable name of truncation time.}

            
\item{min.dn, min.pr, max.pr}{These values are used to set, depending on the model used for modelling, the minimum and maximum allowed 
                              for the densities and probabilities. These
                              parameters are employed to avoid potential overflows/underflows in the calculations and the default 
                              values seem to offer a good compromise. Function \code{conv.check()} provides some relevant  
                              diagnostic information which can be used, for example, to check whether the lower bounds 
                              of \code{min.dn} and \code{min.pr} have been reached. So based on this or if the user wishes to do some sensitivity 
                              analysis then this can be easily carried out using these three arguments.
                              However, the user has to be cautious. For instance, it would not make much sense to choose for \code{min.dn} and \code{min.pr} 
                              values bigger than the default ones. Bear in mind that the bounds can be reached for ill-defined models. For 
                              certain distributions/models, if convergence failure occurs and the bounds have been reached then the user
                              can try a sensitivity analysis as mentioned above.}

\item{ygrid.tol}{Tolerance used to choose grid of response values for robust discrete models. Values smaller than 1e-160 are not allowed for.}


}


\details{

The underlying algorithm is described in ?gjrm. 

There are many continuous/discrete distributions to choose from and we plan to include more 
options. Get in touch if you are interested in a particular distribution.

The \code{"GEVlink"} option is used for binary response additive models and is more stable and faster than the \code{R} package \code{bgeva}.
This model has been incorporated into this package to take advantage of the richer set of smoother choices, and of the 
estimation approach. Details on the model can be found in Calabrese, Marra and Osmetti (2016). 

}


\value{
The function returns an object of class \code{gamlss} as described in \code{gamlssObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Aeberhard W.H., Cantoni E., Marra G., Radice R. (2021), Robust Fitting for Generalized Additive Models for Location, Scale and Shape. \emph{Statistics and Computing}, 31(11), 1-16.

Eletti A., Marra G., Quaresma M., Radice R., Rubio F.J. (2022), A Unifying Framework for Flexible Excess Hazard Modeling with Applications in Cancer Epidemiology. \emph{Journal of the Royal Statistical Society Series C}, 71(4), 1044-1062.

Marra G., Farcomeni A., Radice R. (2021), Link-Based Survival Additive Models under Mixed Censoring to Assess Risks of Hospital-Acquired Infections. \emph{Computational Statistics and Data Analysis}, 155, 107092.

Marra G., Radice R. (2017), Bivariate Copula Additive Models for Location, Scale and Shape. \emph{Computational Statistics and Data Analysis}, 112, 99-113.

Ranjbar S., Cantoni E., Chavez-Demoulin V., Marra G., Radice R., Jaton-Ogay K. (2022), Modelling the Extremes of Seasonal Viruses and Hospital Congestion: The Example of Flu in a Swiss Hospital. \emph{Journal of the Royal Statistical Society Series C}, 71(4), 884-905.

Calabrese R., Marra G., Osmetti SA (2016), Bankruptcy Prediction of Small and Medium Enterprises Using a Flexible Binary Generalized Extreme Value Model. \emph{Journal of the Operational Research Society}, 67(4), 604-615.

Marincioni V., Marra G., Altamirano-Medina H. (2018), Development of Predictive Models for the Probabilistic Moisture Risk Assessment of Internal Wall Insulation. \emph{Building and Environment}, 137, 5257-267. 


}



\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{gamlss()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model. In such a situation, the user may use rescaling (see \code{parscale}). However, the user should especially consider
re-specifying/simplifying the model, and/or checking that the chosen distribution fits the response well.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size is low compared to the complexity of the model. 
It is also worth bearing in mind that the use of three parameter distributions requires the data
to be more informative than a situation in which two parameter distributions are used instead.

} 

\seealso{
\code{\link{GJRM-package}}, \code{\link{gamlssObject}}, \code{\link{conv.check}}, \code{\link{summary.gamlss}}
}

\examples{

\dontrun{ 

library(GJRM)

set.seed(0)

n <- 400

x1 <- round(runif(n))
x2 <- runif(n)
x3 <- runif(n)
f1 <- function(x) cos(pi*2*x) + sin(pi*x)
y1 <- -1.55 + 2*x1 + f1(x2) + rnorm(n)

dataSim <- data.frame(y1, x1, x2, x3)
resp.check(y1, "N")

eq.mu <- y1 ~ x1 + s(x2) + s(x3)
eq.s  <-    ~ s(x3)
fl    <- list(eq.mu, eq.s)

out <- gamlss(fl, data = dataSim)

conv.check(out)
res.check(out)

plot(out, eq = 1, scale = 0, pages = 1, seWithMean = TRUE)
plot(out, eq = 2, seWithMean = TRUE)

summary(out)

AIC(out)
BIC(out)


################
# Robust example
################

eq.mu <- y1 ~ x1 + x2 + x3
fl    <- list(eq.mu)

out <- gamlss(fl, data = dataSim, family = "N", robust = TRUE, 
                  rc = 3, lB = -Inf, uB = Inf)

conv.check(out)
summary(out)
rob.const(out, 100)

##

eq.s  <-    ~ x3
fl    <- list(eq.mu, eq.s)

out <- gamlss(fl, data = dataSim, family = "N", robust = TRUE)

conv.check(out)
summary(out)

##

eq.mu <- y1 ~ x1 + s(x2) + s(x3)
eq.s  <-    ~ s(x3)
fl    <- list(eq.mu, eq.s)

out1 <- gamlss(fl, data = dataSim, family = "N", robust = TRUE, 
               sp.method = "efs")

conv.check(out1)
summary(out1)
AIC(out, out1)

plot(out1, eq = 1, all.terms = TRUE, pages = 1, seWithMean = TRUE)
plot(out1, eq = 2, seWithMean = TRUE)


##########################
## GEV link binary example
##########################
# this incorporates the bgeva
# model implemented in the bgeva package
# however this implementation is more general, 
# stable and efficient

set.seed(0)

n <- 400

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1 <- function(x) cos(pi*2*x) + sin(pi*x)
f2 <- function(x) x+exp(-30*(x-0.5)^2)   

y  <- ifelse(-3.55 + 2*x1 + f1(x2) + rnorm(n) > 0, 1, 0)

dataSim <- data.frame(y, x1, x2, x3)

out1 <- gamlss(list(y ~ x1 + x2 + x3), family = "GEVlink", data = dataSim)
out2 <- gamlss(list(y ~ x1 + s(x2) + s(x3)), family = "GEVlink", data = dataSim)

conv.check(out1)
conv.check(out2)
summary(out1)
summary(out2)
AIC(out1, out2)
BIC(out1, out2)

plot(out2, eq = 1, all.terms = TRUE, pages = 1, seWithMean = TRUE)

##################
# prediction of Pr
##################

# Calculate eta (that is, X*model.coef)
# For a new data set the argument newdata should be used

eta <- predict(out2, eq = 1, type = "link")

# extract gev tail parameter

gev.par <- out2$gev.par

# multiply gev tail parameter by eta

gevpeta <- gev.par*eta 
  
# establish for which values the model is defined   

gevpetaIND <- ifelse(gevpeta < -1, FALSE, TRUE) 
gevpeta <- gevpeta[gevpetaIND]
    
# estimate probabilities  

pr <- exp(-(1 + gevpeta)^(-1/gev.par))


###################################
## Flexible survival model examples
###################################

## Simulate proportional hazards data ##

set.seed(0)
n  <- 2000
c  <- runif(n, 3, 8)
u  <- runif(n, 0, 1)
z1 <- rbinom(n, 1, 0.5)
z2 <- runif(n, 0, 1)
t  <- rep(NA, n)

beta_0 <- -0.2357
beta_1 <- 1

f <- function(t, beta_0, beta_1, u, z1, z2){ 
  S_0 <- 0.7 * exp(-0.03*t^1.9) + 0.3*exp(-0.3*t^2.5)
  exp(-exp(log(-log(S_0))+beta_0*z1 + beta_1*z2))-u
}


for (i in 1:n){
   t[i] <- uniroot(f, c(0, 8), tol = .Machine$double.eps^0.5, 
                   beta_0 = beta_0, beta_1 = beta_1, u = u[i], 
                   z1 = z1[i], z2 = z2[i], extendInt = "yes" )$root
}

delta   <- ifelse(t < c, 1, 0)
u       <- apply(cbind(t, c), 1, min)
dataSim <- data.frame(u, delta, z1, z2)
1-mean(delta) # average censoring rate

# log(u) helps obtaining smoother hazards

out <- gamlss(list(u ~ s(log(u), bs = "mpi") + z1 + s(z2) ), data = dataSim, 
              family = "-cloglog", cens = delta)
res.check(out)
summary(out)
AIC(out)
BIC(out)
plot(out, eq = 1, scale = 0, pages = 1)
haz.surv(out, newdata = data.frame(z1 = 0, z2 = 0), shade = TRUE, 
        n.sim = 1000, baseline = TRUE)
haz.surv(out, type = "haz", newdata = data.frame(z1 = 0, z2 = 0), 
        shade = TRUE, n.sim = 1000, baseline = TRUE)

# library(mgcv)
# out1 <- mgcv::gam(u ~ z1 + s(z2), family = cox.ph(), 
#                   data = dataSim, weights = delta)
# summary(out1)
# estimates of z1 and s(z2) are
# nearly identical between out and out1 

#####################################
## Simulate proportional odds data ##
#####################################

set.seed(0)

n <- 2000
c <- runif(n, 4, 8)
u <- runif(n, 0, 1)
z <- rbinom(n, 1, 0.5)
beta_0 <- -1.05
t <- rep(NA, n)

f <- function(t, beta_0, u, z){ 
  S_0 <- 0.7 * exp(-0.03*t^1.9) + 0.3*exp(-0.3*t^2.5)
  1/(1 + exp(log((1-S_0)/S_0)+beta_0*z))-u
}


for (i in 1:n){
    t[i] <- uniroot(f, c(0, 8), tol = .Machine$double.eps^0.5, 
                    beta_0 = beta_0, u = u[i], z = z[i], 
                    extendInt="yes" )$root
}

delta   <- ifelse(t < c,1, 0)
u       <- apply(cbind(t, c), 1, min)
dataSim <- data.frame(u, delta, z)
1-mean(delta) # average censoring rate

out <- gamlss(list(u ~ s(log(u), bs = "mpi") + z ), data = dataSim, 
              family = "-logit", cens = delta)
res.check(out)
summary(out)
AIC(out)
BIC(out)
plot(out, eq = 1, scale = 0)
haz.surv(out, newdata = data.frame(z = 0), shade = TRUE, n.sim = 1000,
        baseline = TRUE)
haz.surv(out, type = "haz", newdata = data.frame(z = 0), 
        shade = TRUE, n.sim = 1000)
             
                          
#############################
## Mixed censoring example ##
#############################             
             
f1 <- function(t, u, z1, z2, z3, z4, s1, s2){ 

    S_0 <- 0.7 * exp(-0.03*t^1.8) + 0.3*exp(-0.3*t^2.5)
   
    exp( -exp(log(-log(S_0)) + 1.3*z1 + 0.5*z2 + s1(z3) + s2(z4)  ) ) - u   
            
  }
    
datagen <- function(n, z1, z2, z3, z4, s1, s2, f1){
  
  u <- runif(n, 0, 1)
  t <- rep(NA, n)
  
  for (i in 1:n) t[i] <- uniroot(f1, c(0, 100), tol = .Machine$double.eps^0.5, 
                                 u = u[i], s1 = s1, s2 = s2, z1 = z1[i], z2 = z2[i], 
                                 z3 = z3[i], z4 = z4[i], extendInt = "yes")$root
 
  c1 <-      runif(n, 0, 2)
  c2 <- c1 + runif(n, 0, 6) 
  
  df <- data.frame(u1 = t, u2 = t, cens = character(n), stringsAsFactors = FALSE)

for (i in 1:n){

  if(t[i] <= c1[i]) {
        df[i, 1] <- c1[i]
        df[i, 2] <- NA
        df[i, 3] <- "L"
       
  }else if(c1[i] < t[i] && t[i] <= c2[i]){
        df[i, 1] <- c1[i]
        df[i, 2] <- c2[i]
        df[i, 3] <- "I"
        
  }else if(t[i] > c2[i]){
        df[i, 1] <- c2[i]
        df[i, 2] <- NA
        df[i, 3] <- "R"}

}

uncens <- (df[, 3] \%in\% c("L", "I")) + (rbinom(n, 1, 0.2) == 1) == 2 

df[uncens, 1] <- t[uncens]
df[uncens, 2] <- NA
df[uncens, 3] <- "U"

dataSim <- data.frame(u1 = df$u1, u2 = df$u2, cens = as.factor(df$cens), z1, z2, z3, z4, t)
dataSim
  
}

set.seed(0)

n      <- 1000
SigmaC <- matrix(0.5, 4, 4); diag(SigmaC) <- 1
cov    <- rMVN(n, rep(0,4), SigmaC)
cov    <- pnorm(cov)
z1     <- round(cov[, 1])
z2     <- round(cov[, 2])
z3     <- cov[, 3]
z4     <- cov[, 4]
s1     <- function(x) -0.075*exp(3.2 * x) 
s2     <- function(x) sin(2*pi*x) 
 
eq1    <- u1 ~ s(log(u1), bs = "mpi") + z1 + z2 + s(z3) + s(z4)

dataSim <- datagen(n, z1, z2, z3, z4, s1, s2, f1)

out <- gamlss(list(eq1), data = dataSim, family = "-cloglog", 
              cens = cens, type.cen = "mixed", ub.t = "u2")

conv.check(out)
summary(out)
plot(out, eq = 1, scale = 0, pages = 1)             
  
ndf <- data.frame(z1 = 1, z2 = 0, z3 = 0.2, z4 = 0.5)

haz.surv(out, eq = 1, newdata = ndf, type = "surv")
haz.surv(out, eq = 1, newdata = ndf, type = "haz", n.sim = 1000)         
}

}

\keyword{regression modelling} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{distribution}
\keyword{survival data}


