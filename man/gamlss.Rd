\name{gamlss}
\alias{gamlss}

\title{Generalised Additive Models for Location, Scale and Shape}

\description{ 
\code{gamlss} fits flexible univariate regression models with several continuous and discrete distributions, and types of covariate 
effects. The purpose of this function was only to provide, in some cases, starting values 
for the simultaneous models in the package, but it has now been made available in the form of a proper function should the user wish to fit 
univariate models using the general estimation approach of this package. The distributions implemented here 
have been parametrised according to Rigby and Stasinopoulos (2005).
 
}



\usage{
gamlss(formula, data = list(), weights = NULL, subset = NULL,  
       margin = "N", surv = FALSE, cens = NULL,
       robust = FALSE, rc = 3, lB = NULL, uB = NULL, infl.fac = 1, 
       rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
       gc.l = FALSE, parscale, extra.regI = "t", gev.par = -0.25,
       chunk.size = 10000, k.tvc = 0, knots = NULL)
}



\arguments{ 

\item{formula}{List of equations. This should contain one or more equations.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{gamlss} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}


\item{margin}{Possible distributions are normal ("N"), normal where 
               sigma2 corresponds to 
               the standard deviation instead of the variance ("N2"),  
               log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution), Poisson ("PO"), zero truncated 
               Poisson ("ZTP"), negative binomial - type I ("NBI"), negative 
               binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG"), generalised extreme value link function ("GEVlink", this 
               is used for binary responses and is more stable and faster than the \code{R} package \code{bgeva}).}
 
\item{surv}{If \code{TRUE} then a survival model is fitted. Here margin can be "PH" (generalised proportional hazards), "PO" (generalised proportional odds), 
            "probit" (generalised probit).}

\item{cens}{Binary censoring indicator. This is required when \code{surv = TRUE}. This variable has to be equal to 1 if the event occurred and 0 otherwise.}
 
\item{robust}{If \code{TRUE} then the robust version of the model is fitted.}
             
\item{rc}{Robust constant.}       

\item{lB, uB}{Lower and upper bounds for integrals, when robust estimation is employed. If not provided then standard bounds are used.}  
 
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

\item{extra.regI}{If "t" then regularization as from \code{trust} is applied to the information matrix if needed. 
                  If different from "t" then extra regularization is applied via the options "pC" (pivoted Choleski - this
                  will only work when the information matrix is semi-positive or positive definite) and "sED" (symmetric eigen-decomposition).}  

\item{gev.par}{GEV link parameter.}
      
\item{chunk.size}{This is used for discrete robust models.}  

\item{k.tvc}{Only used for tvc ps smoothers when using survival models.}  

\item{knots}{Optional list containing user specified knot values to be used for basis construction.}
  
}


\details{

The underlying algorithm is described in ?SemiParBIV. 

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

Marra G. and Radice R. (2017), Bivariate Copula Additive Models for Location, Scale and Shape. \emph{Computational Statistics and Data Analysis}, 112, 99-113.

Rigby R.A. and Stasinopoulos D.M. (2005). Generalized additive models for location, scale and shape (with discussion). \emph{Journal of the Royal Statistical Society, Series C}, 54(3), 507-554.

Calabrese R., Marra G. and Osmetti SA (2016), Bankruptcy Prediction of Small and Medium Enterprises Using a Flexible Binary Generalized Extreme Value Model. \emph{Journal of the Operational Research Society}, 67(4), 604-615.

}



\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{gamlss()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model. In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). However, the user should especially consider
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
eq.s2 <-    ~ s(x3)
fl    <- list(eq.mu, eq.s2)

out <- gamlss(fl, data = dataSim)

conv.check(out)
post.check(out)

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

out <- gamlss(fl, data = dataSim, margin = "N", robust = TRUE, 
                  rc = 3, lB = -Inf, uB = Inf)

conv.check(out)
summary(out)

##

eq.s2 <-    ~ x3
fl    <- list(eq.mu, eq.s2)

out <- gamlss(fl, data = dataSim, margin = "N", robust = TRUE)

conv.check(out)
summary(out)

##

eq.mu <- y1 ~ x1 + s(x2) + s(x3)
eq.s2 <-    ~ s(x3)
fl    <- list(eq.mu, eq.s2)

out1 <- gamlss(fl, data = dataSim, margin = "N", robust = TRUE)

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
# however this implementation is more general 
# stable and efficient

set.seed(0)

n <- 400

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1 <- function(x) cos(pi*2*x) + sin(pi*x)
f2 <- function(x) x+exp(-30*(x-0.5)^2)   

y  <- ifelse(-3.55 + 2*x1 + f1(x2) + rnorm(n) > 0, 1, 0)

dataSim <- data.frame(y, x1, x2, x3)

out1 <- gamlss(list(y ~ x1 + x2 + x3), margin = "GEVlink", data = dataSim)
out2 <- gamlss(list(y ~ x1 + s(x2) + s(x3)), margin = "GEVlink", data = dataSim)

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

##################################
## Flexible survival model example
##################################

library(GJRM)

########################################
## Simulate proportional hazards data ##
########################################

set.seed(0)
n  <- 2000
c  <- runif(n, 3, 8)
u  <- runif(n, 0, 1)
z1 <- rbinom(n, 1, 0.5)
z2 <- runif(n, 0, 1)
t <- rep(NA, n)

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

out <- gamlss(list(u ~ z1 + s(z2) + s(u, bs = "mpi") ), data = dataSim, 
              surv = TRUE, margin = "PH", cens = delta)
post.check(out)
summary(out)
AIC(out)
BIC(out)
plot(out, eq = 1, scale = 0, pages = 1)
hazsurv.plot(out, newdata = data.frame(z1 = 0, z2 = 0), shade = TRUE, n.sim = 1000)
hazsurv.plot(out, type = "hazard", newdata = data.frame(z1 = 0, z2 = 0), 
             shade = TRUE, n.sim = 1000)

out1 <- gam(u ~ z1 + s(z2), family = cox.ph(), 
            data = dataSim, weights = delta)
summary(out1)
# estimates of z1 and s(z2) are
# nearly identical between out and out1 

# note that the Weibull is implemented as AFT
# as using the PH parametrisation makes
# computation unstable
out2 <- gamlss(list(u ~ z1 + s(z2) ), data = dataSim, surv = TRUE, 
               margin = "WEI", cens = delta)
 
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

out <- gamlss(list(u ~ z + s(u, bs = "mpi") ), data = dataSim, surv = TRUE, 
              margin = "PO", cens = delta)
post.check(out)
summary(out)
AIC(out)
BIC(out)
plot(out, eq = 1, scale = 0)
hazsurv.plot(out, newdata = data.frame(z =0), shade = TRUE, n.sim = 1000)
hazsurv.plot(out, type = "hazard", newdata = data.frame(z = 0), 
             shade = TRUE, n.sim = 1000)

# note that the Fisk is implemented as AFT
# as using the PH parametrisation makes
# computation unstable
out1 <- gamlss(list(u ~ z), data = dataSim, surv = TRUE, 
               margin = "FISK", cens = delta)

}

}

\keyword{regression modelling} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{distribution}


