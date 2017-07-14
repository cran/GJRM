\name{copulaReg}
\alias{copulaReg}

\title{Copula Models with Continuous/Discrete/Survival Margins}

\description{ 
\code{copulaReg} fits flexible copula models with continuous/discrete/survival margins with several types of covariate 
effects, copula and marginal distributions. 
 
}



\usage{
copulaReg(formula, data = list(), weights = NULL, subset = NULL,  
                 BivD = "N", margins = c("N","N"), dof = 3,  
                 surv = FALSE, cens1 = NULL, cens2 = NULL,
                 gamlssfit = FALSE, fp = FALSE, infl.fac = 1, 
                 rinit = 1, rmax = 100, 
                 iterlimsp = 50, tolsp = 1e-07,
                 gc.l = FALSE, parscale, extra.regI = "t")
}



\arguments{ 

\item{formula}{In the basic setup this will be a list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms 
               are used to specify smooth functions of 
                  predictors. For the case of more than two equations see the example below and the documentation 
                  of \code{SemiParBIV()} for more details. When one outcome is binary and the other continuous/discrete 
                  then the first equation MUST refer to the binary outcome whereas 
                  the second to the continuous/discrete one. When one outcome is discrete and the other continuous 
                  then the first equation MUST refer to the discrete one.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{copulaReg} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}


\item{margins}{It indicates the distributions used for the two margins. Possible distributions are normal ("N"), normal where 
               sigma2 corresponds to 
               the standard deviation instead of the variance ("N2"),  
               log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution), Poisson ("PO"), zero truncated 
               Poisson ("ZTP"), negative binomial - type I ("NBI"), negative 
               binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG"). When the first equation is binary then
               possible link functions are "probit", "logit", "cloglog"".}
  
\item{dof}{If \code{BivD = "T"} then the degrees of freedom can be set to a value greater than 2 and smaller than 249. Only for continuous margins,
           this will be taken as a starting value and the dof estiamted from the data.}   
 
\item{surv}{If \code{TRUE} then a bivariate survival model is fitted.}    

\item{cens1}{Binary censoring indicator 1. This is required when \code{surv = TRUE}.}

\item{cens2}{Binary censoring indicator 2. This is required when \code{surv = TRUE}.}

\item{gamlssfit}{If \code{gamlssfit = TRUE} then gamlss univariate models are also fitted. This is useful for obtaining 
                 starting values, for instance.}               

\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "AMH", "FGM", "T", "PL", "HO" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), 
            survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, Ali-Mikhail-Haq,
             Farlie-Gumbel-Morgenstern, Student-t with \code{dof}, Plackett, Hougaard. Each of the Clayton, Joe and Gumbel copulae is allowed to be mixed with a rotated version of the same
family. The options are: "C0C90", "C0C270", "C180C90", "C180C270", "G0G90", "G0G270", "G180G90",
"G180G270", "J0J90", "J0J270", "J180J90" and "J180J270". This allows the user to model negative and positive tail dependencies.}            

\item{fp}{If \code{TRUE} then a fully parametric model with unpenalised regression splines if fitted. See the example below.}              

\item{infl.fac}{Inflation factor for the model degrees of freedom in the approximate AIC. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{rinit}{Starting trust region radius. The trust region radius is adjusted as the algorithm proceeds. See the documentation
             of \code{trust} for further details.}

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

\item{extra.regI}{If "t" then regularization as from \code{trust} is applied to the information matrix if needed. 
                  If different from "t" then extra regularization is applied via the options "pC" (pivoted Choleski - this
                  will only work when the information matrix is semi-positive or positive definite) and "sED" (symmetric eigen-decomposition).}  
     
                      
}


\details{

The underlying algorithm is described in ?SemiParBIV. 

There are many continuous/discrete/survival distributions and copula functions to choose from and we plan to include more 
options. Get in touch if you are interested in a particular distribution.


}


\value{
The function returns an object of class \code{copulaReg} as described in \code{copulaRegObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G. and Radice R. (2017), Bivariate Copula Additive Models for Location, Scale and Shape. \emph{Computational Statistics and Data Analysis}, 112, 99-113.

}



\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{copulaReg()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model.
In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). Using \code{gamlssfit = TRUE} is typically more effective than the first two options as
this will provide better calibrated starting values as compared to those obtained from the default starting value procedure.
The default option is, however, \code{gamlssfit = FALSE} only because it tends to be computationally cheaper and because the 
default starting value procedure has typically been found to do a satisfactory job in most cases. 
(The results obtained when using 
\code{gamlssfit = FALSE} and \code{gamlssfit = TRUE} could also be compared to check if starting values make any difference.)

The above suggestions may help, especially the latter option. However, the user should also consider
re-specifying/simplifying the model, and/or using a diferrent dependence structure and/or checking that the chosen marginal 
distributions fit the responses well.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size is low compared to the complexity of the model. Examples
of misspecification include using a Clayton copula rotated by 90 degrees when a positive
association between the margins is present instead, using marginal distributions that do not fit
the responses, and 
employing a copula which does not accommodate the type and/or strength of
the dependence between the margins (e.g., using AMH when the association between the margins is strong). When using 
smooth functions, if the covariate's values are too sparse then convergence may be affected by this.
It is also worth bearing in mind that the use of three parameter marginal distributions requires the data
to be more informative than a situation in which two parameter distributions are used instead.


} 

\seealso{
\code{\link{adjCov}}, \code{\link{VuongClarke}}, \code{\link{GJRM-package}}, \code{\link{copulaRegObject}}, \code{\link{conv.check}}, \code{\link{summary.copulaReg}}
}

\examples{

library(GJRM)

\dontrun{

############
## EXAMPLE 1
## Generate data
## Correlation between the two equations 0.5 - Sample size 400 

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- -1.55 + 2*x1    + f1(x2) + u[,1]
y2 <- -0.25 - 1.25*x1 + f2(x2) + u[,2]

dataSim <- data.frame(y1, y2, x1, x2, x3)

resp.check(y1, "N")
resp.check(y2, "N")

eq.mu.1     <- y1 ~ x1 + s(x2) + s(x3)
eq.mu.2     <- y2 ~ x1 + s(x2) + s(x3)
eq.sigma2.1 <-    ~ 1
eq.sigma2.2 <-    ~ 1
eq.theta    <-    ~ x1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma2.1, eq.sigma2.2, eq.theta)

# the order above is the one to follow when
# using more than two equations

out  <- copulaReg(fl, data = dataSim)

conv.check(out)
post.check(out)
summary(out)
AIC(out)
BIC(out)
jc.probs(out, 1.4, 2.3, intervals = TRUE)[1:4,]

############
## EXAMPLE 2
############
## Generate data with one endogenous binary variable 
## and continuous outcome

set.seed(0)

n <- 1000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

cov   <- rMVN(n, rep(0,2), Sigma)
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <-        -0.25 - 1.25*y1 + f2(x2) + u[,2] 

dataSim <- data.frame(y1, y2, x1, x2)


## RECURSIVE Model

rc <- resp.check(y2, margin = "N", print.par = TRUE, loglik = TRUE)
AIC(rc); BIC(rc)

out <- copulaReg(list(y1 ~ x1 + x2, 
                      y2 ~ y1 + x2), 
                 data = dataSim, margins = c("probit","N"))
conv.check(out)                        
summary(out)
post.check(out)

## SEMIPARAMETRIC RECURSIVE Model

eq.mu.1   <- y1 ~ x1 + s(x2) 
eq.mu.2   <- y2 ~ y1 + s(x2)
eq.sigma2 <-    ~ 1
eq.theta  <-    ~ 1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma2, eq.theta)

out <- copulaReg(fl, data = dataSim, 
                 margins = c("probit","N"), gamlssfit = TRUE)
conv.check(out)                        
summary(out)
post.check(out)
jc.probs(out, 1, 1.5, intervals = TRUE)[1:4,]
AT(out, nm.end = "y1")
AT(out, nm.end = "y1", type = "univariate")


#
#

############
## EXAMPLE 3
############
## Generate data with one endogenous continuous exposure 
## and binary outcome

set.seed(0)

n <- 1000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

cov   <- rMVN(n, rep(0,2), Sigma)
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2) 

y1 <-        -0.25 - 2*x1    + f2(x2) + u[,2] 
y2 <- ifelse(-0.25 - 0.25*y1 + f1(x2) + u[,1] > 0, 1, 0)

dataSim <- data.frame(y1, y2, x1, x2)

eq.mu.1   <- y2 ~ y1 + s(x2) 
eq.mu.2   <- y1 ~ x1 + s(x2)
eq.sigma2 <-    ~ 1
eq.theta  <-    ~ 1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma2, eq.theta)

out <- copulaReg(fl, data = dataSim, 
                 margins = c("probit","N"))
conv.check(out)                        
summary(out)
post.check(out)
AT(out, nm.end = "y1")
AT(out, nm.end = "y1", type = "univariate")
RR(out, nm.end = "y1", rr.plot = TRUE)
RR(out, nm.end = "y1", type = "univariate")
OR(out, nm.end = "y1", or.plot = TRUE)
OR(out, nm.end = "y1", type = "univariate")

#
#

############
## EXAMPLE 4
############
## Survival model

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

delta1  <- ifelse(t < c, 1, 0)
u1      <- apply(cbind(t, c), 1, min)
dataSim <- data.frame(u1, delta1, z1, z2)


c <- runif(n, 4, 8)
u <- runif(n, 0, 1)
z <- rbinom(n, 1, 0.5)
beta_0 <- -1.05
t      <- rep(NA, n)

f <- function(t, beta_0, u, z){ 
  S_0 <- 0.7 * exp(-0.03*t^1.9) + 0.3*exp(-0.3*t^2.5)
  1/(1 + exp(log((1-S_0)/S_0)+beta_0*z))-u
}



for (i in 1:n){
    t[i] <- uniroot(f, c(0, 8), tol = .Machine$double.eps^0.5, 
                    beta_0 = beta_0, u = u[i], z = z[i], 
                    extendInt="yes" )$root
}

delta2 <- ifelse(t < c,1, 0)
u2     <- apply(cbind(t, c), 1, min)
dataSim$delta2 <- delta2
dataSim$u2     <- u2
dataSim$z      <- z



eq1 <- u1 ~ z1 + s(z2) + s(u1, bs = "mpi")
eq2 <- u2 ~ z + s(u2, bs = "mpi") 
eq3 <-    ~ s(z2)

out <- copulaReg(list(eq1, eq2), data = dataSim, surv = TRUE,
                 margins = c("PH", "PO"), 
                 cens1 = delta1, cens2 = delta2)
                 
# PH margin fit can also be compared with cox.ph from mgcv
                 
conv.check(out)
post.check(out)
summary(out)
AIC(out); BIC(out)
plot(out, eq = 1, scale = 0, pages = 1)
plot(out, eq = 2, scale = 0, pages = 1)

hazsurv.plot(out, eq = 1, newdata = data.frame(z1 = 0, z2 = 0), 
             shade = TRUE, n.sim = 1000)
hazsurv.plot(out, eq = 1, newdata = data.frame(z1 = 0, z2 = 0), 
             shade = TRUE, n.sim = 1000, type = "hazard")
hazsurv.plot(out, eq = 2, newdata = data.frame(z = 0), 
             shade = TRUE, n.sim = 1000)
hazsurv.plot(out, eq = 2, newdata = data.frame(z = 0), 
             shade = TRUE, n.sim = 1000, type = "hazard")
 
jc.probs(out, type = "bivariate", intervals = TRUE)[1:5,]
 
newd0 <- newd1 <- data.frame(z = 0, z1 = mean(dataSim$z1), 
                             z2 = mean(dataSim$z2), 
                             u1 =  mean(dataSim$u1) + 1, 
                             u2 =  mean(dataSim$u2) + 1) 
newd1$z <- 1                   

jc.probs(out, type = "bivariate", newdata = newd0, intervals = TRUE)
jc.probs(out, type = "bivariate", newdata = newd1, intervals = TRUE)

out1 <- copulaReg(list(eq1, eq2, eq3), data = dataSim, surv = TRUE,
                  margins = c("PH", "PO"), 
                  cens1 = delta1, cens2 = delta2, gamlssfit = TRUE) 

eq1 <- u1 ~ z1 + s(z2)
eq2 <- u2 ~ z  
eq3 <-    ~ s(z2)    

# note that Weibull is implemented as AFT model
out2 <- copulaReg(list(eq1, eq2, ~ 1, ~ 1, eq3), data = dataSim, surv = TRUE,
                  margins = c("WEI", "WEI"), 
                  cens1 = delta1, cens2 = delta2)                   
                   
}
}

\keyword{generalised joint regression modelling}
\keyword{smooth} 
\keyword{regression spline}
\keyword{copula}
\keyword{marginal distribution}


