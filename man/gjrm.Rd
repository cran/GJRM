\name{gjrm}
\alias{gjrm}

\title{Generalised Joint Regression Models with Binary/Continuous/Discrete/Survival Margins}

\description{ 
\code{gjrm} fits flexible joint models with binary/continuous/discrete/survival margins, with several types of covariate 
effects, copula and marginal distributions. 
 
}



\usage{
gjrm(formula, data = list(), weights = NULL, subset = NULL,  
     copula = "N", copula2 = "N", margins, model, dof = 3, dof2 = 3, ordinal = FALSE,  
     surv = FALSE, cens1 = NULL, cens2 = NULL, cens3 = NULL, dep.cens = FALSE,
     upperBt1 = NULL, upperBt2 = NULL,   
     gamlssfit = FALSE, fp = FALSE, infl.fac = 1, 
     rinit = 1, rmax = 100, 
     iterlimsp = 50, tolsp = 1e-07,
     gc.l = FALSE, parscale, extra.regI = "t",
     k1.tvc = 0, k2.tvc = 0, knots = NULL,
     penCor = "unpen", sp.penCor = 3, 
     Chol = FALSE, gamma = 1, w.alasso = NULL,
     drop.unused.levels = TRUE, 
     min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999)
}



\arguments{ 

\item{formula}{In the basic setup this will be a list of two (or three) formulas, one for equation 1, the other for equation 2 and another one
               for equation 3 if a trivariate model is fitted to the data. Otherwise, more equations can be used depending on the 
               number of distributional parameters. \code{s} terms 
               are used to specify smooth functions of 
                  predictors; see the documentation of \code{mgcv} for further 
                  details on formula specifications. Note that if a selection model is employed (that is, \code{model = "BSS"}
                  or \code{model = "TSS"}) then the first formula (and the second as well for trivariate models) MUST refer to the selection equation(s).
                  When one outcome is binary and the other continuous/discrete 
                  then the first equation should refer to the binary outcome whereas 
                  the second to the continuous/discrete one. When one outcome is discrete and the other continuous 
                  then the first equation has to be the discrete one.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{gjrm} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{copula}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "GAL0", "GAL90", "GAL180", "GAL270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "AMH", "FGM", "T", "PL", "HO" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), 
            survival Clayton,
             rotated Clayton (270 degrees), Galambos, rotated Galambos (90 degrees), 
            survival Galambos,
             rotated Galambos (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, Ali-Mikhail-Haq,
             Farlie-Gumbel-Morgenstern, Student-t with \code{dof}, Plackett, Hougaard. Each of the Clayton, Galambos, Joe and Gumbel copulae is allowed to be mixed with a rotated version of the same
family. The options are: "C0C90", "C0C270", "C180C90", "C180C270",  "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270", "G0G90", "G0G270", "G180G90",
"G180G270", "J0J90", "J0J270", "J180J90" and "J180J270". This allows the user to model negative and positive tail dependencies.}  

\item{copula2}{As above but used only for Roy models.}

\item{margins}{It indicates the distributions used for the two or three margins. Possible distributions are normal ("N"),  
               log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), generelised Pareto ("GP"), 
               generelised Pareto II ("GPII") where the shape parameter is forced to be > -0.5,
               generelised Pareto (with orthogonal parametrisation) ("GPo") where the shape parameter is forced to be > -0.5,
               discrete generelised Pareto ("DGP"), 
               discrete generelised Pareto II ("DGPII") where the shape parameter is forced to be positive, discrete generelised Pareto derived
               under the scenario in which shape = 0 ("DGP0"), logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution), Poisson ("PO"), zero truncated 
               Poisson ("ZTP"), negative binomial - type I ("NBI"), negative 
               binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG"). If the responses are binary then
               possible link functions are "probit", "logit", "cloglog". For survival models, the margins can be proportional hazars ("PH"), odds ("PO")
               or "probit".}

\item{model}{Possible values are "B" (bivariate model), "T" (trivariate model) 
             "BSS" (bivariate model with non-random sample selection), "TSS" (trivariate model with double non-random sample selection),
             "TESS" (trivariate model with endogeneity and non-random sample selection),
             "BPO" (bivariate model with partial observability)
             and "BPO0" (bivariate model with partial observability and zero correlation).
             Options "T", "TESS" and "TSS" are currently for trivariate binary models only. "BPO" and "BPO0" are for bivariate binary models only.    
             "ROY" is for the Roy switching regression model.
             } 
  
\item{dof}{If \code{copula = "T"} then the degrees of freedom can be set to a value greater than 2 and smaller than 249. Only for continuous margins,
                 this will be taken as a starting value and the dof estimated from the data.}   
                 
\item{dof2}{As above but used only for Roy models.} 
 
\item{ordinal}{If \code{TRUE} then the ordinal model is employed.}    
 
 
\item{surv}{If \code{TRUE} then a bivariate survival model is fitted.}    

\item{cens1}{Censoring indicator for the first equation. This is required when \code{surv = TRUE}. For the case of right censored data only, this variable can be equal to 1 if the event occurred 
            and 0 otherwise. However, if there are several censoring mechanisms then \code{cens} will have to be specified as a factor 
            variable made up of four possible categories: \code{I} for interval, \code{L} for left, \code{R} for right, and \code{U} for uncensored.}

\item{cens2}{Same as above but for the second equation.}


\item{cens3}{Binary censoring indicator employed only when \code{surv = TRUE}, \code{dep.cens = TRUE} and administrative censoring is present.}

\item{dep.cens}{If TRUE then the dependence censored model is employed.} 

\item{upperBt1, upperBt2}{Variable names of right/upper bounds when interval censoring is present.}

\item{gamlssfit}{If \code{gamlssfit = TRUE} then gamlss univariate models are also fitted. This is useful for obtaining 
                 starting values, for instance.}               


\item{fp}{If \code{TRUE} then a fully parametric model with unpenalised regression splines if fitted. See the Example 2 below.}              

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
     
\item{k1.tvc, k2.tvc}{Experimental. Only used for tvc ps smoothers when using survival models.}  
 
\item{knots}{Optional list containing user specified knot values to be used for basis construction.}
  

\item{penCor}{This and the arguments below are only for trivariate binary models. Type of penalty for 
              correlation coefficients. Possible values are "unpen", "lasso", "ridge", "alasso".}

\item{sp.penCor}{Starting value for smoothing parameter of \code{penCor}.}

\item{Chol}{If \code{TRUE} then the Cholesky method instead of the eigenvalue method is employed for the correlation matrix.}
                   
\item{gamma}{Inflation factor used only for the alasso penalty.}

\item{w.alasso}{When using the alasso penalty a weight vector made up of three values must be provided.} 
 
\item{drop.unused.levels}{By default unused levels are dropped from factors before fitting. For some smooths involving factor variables 
                          this may have to be turned off (only use if you know what you are doing).} 
 
 \item{min.dn, min.pr, max.pr}{These values are used to set, depending on the model used for modelling, the minimum and maximum allowed 
                              for the densities and probabilities; recall that the margins of copula models have to be in the range (0,1). These
                              parameters are employed to avoid potential overflows/underflows in the calculations and the default 
                              values seem to offer a good compromise. Function \code{conv.check()} provides some relevant  
                              diagnostic information which can be used, for example, to check whether the lower bounds 
                              of \code{min.dn} and \code{min.pr} have been reached. So based on this or if the user wishes to do some sensitivity 
                              analysis then this can be easily carried out using these three arguments.
                              However, the user has to be cautious. For instance, it would not make much sense to choose for \code{min.dn} and \code{min.pr} 
                              values bigger than the default ones. Bear in mind that the bounds can be reached for ill-defined models. For 
                              certain distributions/models, if convergence failure occurs and the bounds have been reached then the user
                              can try a sensitivity analysis as mentioned above.}

 
}


\details{

The joint models considered by this function consist of two or three model equations which depend on flexible linear predictors and
whose dependence between the responses is modelled through one or more parameters of a chosen multivariate distribution. The additive predictors of the 
equations are flexibly specified using 
parametric components and smooth functions of covariates. The same can be done for the dependence parameter(s) if it makes sense.
Estimation is achieved within a penalized likelihood framework with integrated automatic multiple smoothing parameter selection. The use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support 
from the data. The trade-off between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimise an approximate AIC.

For sample selection models, if there are factors in the model then before fitting the user has to ensure 
that the numbers of factor variables' levels in the selected sample 
are the same as those in the complete dataset. Even if a model could be fitted in such a situation,
the model may produce fits which are 
not coherent with the nature of the correction sought. As an example consider the 
situation in which the complete dataset contains a factor variable with five levels and that only three of them   
appear in the selected sample. For the outcome equation (which is the one of interest) only three levels of such variable 
exist in the population, but their effects will be corrected for non-random selection using a selection equation 
in which five levels exist instead. Having differing numbers of factors' levels between complete and selected samples will also  
make prediction not feasible (an aspect which may be particularly important for selection models);
clearly it is not possible to predict the response of interest for the missing entries using a 
dataset that contains all levels of a factor variable but using an outcome model 
estimated using a subset of these levels. 

There are many continuous/discrete/survival distributions and copula functions to choose from and we plan to include more 
options. Get in touch if you are interested in a particular distribution.


}


\value{
The function returns an object of class \code{gjrm} as described in \code{gjrmObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

See help("GJRM-package").

}



\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{gjrm()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model.
In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). Using \code{gamlssfit = TRUE} is typically more effective than the first two options as
this will provide better calibrated starting values as compared to those obtained from the default starting value procedure.
The default option is, however, \code{gamlssfit = FALSE} only because it tends to be computationally cheaper and because the 
default procedure has typically been found to do a satisfactory job in most cases. 
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


In the contexts of endogeneity and non-random sample selection, extra attention is required when specifying
the dependence parameter as a function of covariates. This is because in these situations the dependence parameter mainly models the 
association between the unobserved confounders in the two equations. Therefore, this option would make sense when it 
is believed that the  
strength of the association between the unobservables in the two equations varies based on some grouping factor or across geographical 
areas, for instance. In any case, a clear rationale is typically needed in such cases.

} 

\seealso{
\code{\link{adjCov}}, \code{\link{VuongClarke}}, \code{\link{GJRM-package}}, \code{\link{gjrmObject}}, \code{\link{conv.check}}, \code{\link{summary.gjrm}}
}

\examples{

library(GJRM)

####################################
####################################
####################################
# JOINT MODELS WITH BINARY MARGINS #
####################################
####################################
####################################

############
## Example 1
############

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*x1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1, y2, x1, x2, x3)

## CLASSIC BIVARIATE PROBIT

out  <- gjrm(list(y1 ~ x1 + x2 + x3, 
                       y2 ~ x1 + x2 + x3), 
                       data = dataSim, 
                       margins = c("probit", "probit"),
                       model = "B")
conv.check(out)
summary(out)
AIC(out)
BIC(out)

\dontrun{ 


## BIVARIATE PROBIT with Splines

out  <- gjrm(list(y1 ~ x1 + s(x2) + s(x3), 
                  y2 ~ x1 + s(x2) + s(x3)),  
                  data = dataSim,
                  margins = c("probit", "probit"),
                  model = "B")
conv.check(out)
summary(out)
AIC(out)


## estimated smooth function plots - red lines are true curves

x2 <- sort(x2)
f1.x2 <- f1(x2)[order(x2)] - mean(f1(x2))
f2.x2 <- f2(x2)[order(x2)] - mean(f2(x2))
f3.x3 <- rep(0, length(x3))

par(mfrow=c(2,2),mar=c(4.5,4.5,2,2))
plot(out, eq = 1, select = 1, seWithMean = TRUE, scale = 0)
lines(x2, f1.x2, col = "red")
plot(out, eq = 1, select = 2, seWithMean = TRUE, scale = 0)
lines(x3, f3.x3, col = "red")
plot(out, eq = 2, select = 1, seWithMean = TRUE, scale = 0)
lines(x2, f2.x2, col = "red")
plot(out, eq = 2, select = 2, seWithMean = TRUE, scale = 0)
lines(x3, f3.x3, col = "red")


## BIVARIATE PROBIT with Splines and 
## varying dependence parameter

eq.mu.1  <- y1 ~ x1 + s(x2)
eq.mu.2  <- y2 ~ x1 + s(x2)
eq.theta <-    ~ x1 + s(x2)

fl <- list(eq.mu.1, eq.mu.2, eq.theta)

outD <- gjrm(fl, data = dataSim,
             margins = c("probit", "probit"),
             model = "B")
             
conv.check(outD)  
summary(outD)
outD$theta

plot(outD, eq = 1, seWithMean = TRUE)
plot(outD, eq = 2, seWithMean = TRUE)
plot(outD, eq = 3, seWithMean = TRUE)
graphics.off()

############
## Example 2
############
## Generate data with one endogenous variable 
## and exclusion restriction

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

cov   <- rMVN(n, rep(0,2), Sigma)
cov   <- pnorm(cov)
x1 <- round(cov[,1]); x2 <- cov[,2]

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2)   

y1 <- ifelse(-1.55 + 2*x1    + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*y1 + f2(x2) + u[,2] > 0, 1, 0)

dataSim <- data.frame(y1, y2, x1, x2)

#
## Testing the hypothesis of absence of endogeneity... 
#

LM.bpm(list(y1 ~ x1 + s(x2), y2 ~ y1 + s(x2)), dataSim, model = "B")


## CLASSIC RECURSIVE BIVARIATE PROBIT

out <- gjrm(list(y1 ~ x1 + x2, 
                 y2 ~ y1 + x2), 
                 data = dataSim,
                 margins = c("probit", "probit"),
                 model = "B")
conv.check(out)                        
summary(out)
AIC(out); BIC(out)

## FLEXIBLE RECURSIVE BIVARIATE PROBIT

out <- gjrm(list(y1 ~ x1 + s(x2), 
                 y2 ~ y1 + s(x2)), 
                 data = dataSim,
                 margins = c("probit", "probit"),
                 model = "B")
conv.check(out)                        
summary(out)
AIC(out); BIC(out)

#
## Testing the hypothesis of absence of endogeneity post estimation... 

gt.bpm(out)

#
## treatment effect, risk ratio and odds ratio with CIs

mb(y1, y2, model = "B")
AT(out, nm.end = "y1", hd.plot = TRUE) 
RR(out, nm.end = "y1") 
OR(out, nm.end = "y1") 
AT(out, nm.end = "y1", type = "univariate") 
re.imp <- imputeCounter(out, m = 10, "y1")
re.imp$AT

## try a Clayton copula model... 

outC <- gjrm(list(y1 ~ x1 + s(x2), 
                  y2 ~ y1 + s(x2)), 
                  data = dataSim, copula = "C0",
                  margins = c("probit", "probit"),
                  model = "B")
conv.check(outC)                         
summary(outC)
AT(outC, nm.end = "y1") 
re.imp <- imputeCounter(outC, m = 10, "y1")
re.imp$AT

## try a Joe copula model... 

outJ <- gjrm(list(y1 ~ x1 + s(x2), 
                  y2 ~ y1 + s(x2)), 
                  data = dataSim, copula = "J0",
                  margins = c("probit", "probit"),
                  model = "B")
conv.check(outJ)
summary(outJ)
AT(outJ, "y1") 
re.imp <- imputeCounter(outJ, m = 10, "y1")
re.imp$AT

VuongClarke(out, outJ)

#
## recursive bivariate probit modelling with unpenalized splines 
## can be achieved as follows

outFP <- gjrm(list(y1 ~ x1 + s(x2, bs = "cr", k = 5), 
                   y2 ~ y1 + s(x2, bs = "cr", k = 6)), 
                   fp = TRUE, data = dataSim,
                   margins = c("probit", "probit"),
                   model = "B")
conv.check(outFP)                            
summary(outFP)

# in the above examples a third equation could be introduced
# as illustrated in Example 1

#
#################
## See also ?meps
#################

############
## Example 3
############
## Generate data with a non-random sample selection mechanism 
## and exclusion restriction

set.seed(0)

n <- 2000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

SigmaC <- matrix(0.5, 3, 3); diag(SigmaC) <- 1
cov    <- rMVN(n, rep(0,3), SigmaC)
cov    <- pnorm(cov)
bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + u[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +         + u[, 2] > 0
yo <- y*(ys > 0)
  
dataSim <- data.frame(y, ys, yo, bi, x1, x2)

#
## Testing the hypothesis of absence of non-random sample selection... 

LM.bpm(list(ys ~ bi + s(x1) + s(x2), yo ~ bi + s(x1)), dataSim, model = "BSS")

# p-value suggests presence of sample selection, hence fit a bivariate model

#
## SEMIPARAMETRIC SAMPLE SELECTION BIVARIATE PROBIT
## the first equation MUST be the selection equation

out <- gjrm(list(ys ~ bi + s(x1) + s(x2), 
                 yo ~ bi + s(x1)), 
                 data = dataSim, model = "BSS",
                 margins = c("probit", "probit"))
conv.check(out)                          
gt.bpm(out)                        

## compare the two summary outputs
## the second output produces a summary of the results obtained when
## selection bias is not accounted for

summary(out)
summary(out$gam2)

## corrected predicted probability that 'yo' is equal to 1

mb(ys, yo, model = "BSS")
prev(out, hd.plot = TRUE)
prev(out, type = "univariate", hd.plot = TRUE)

## estimated smooth function plots
## the red line is the true curve
## the blue line is the univariate model curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)]-mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1.65,0.95)); lines(x1.s, f21.x1, col="red")
par(new = TRUE)
plot(out$gam2, se = FALSE, col = "blue", ylim = c(-1.65,0.95), 
     ylab = "", rug = FALSE)

#
#
## try a Clayton copula model... 

outC <- gjrm(list(ys ~ bi + s(x1) + s(x2), 
                  yo ~ bi + s(x1)), 
                  data = dataSim, model = "BSS", copula = "C0",
                  margins = c("probit", "probit"))
conv.check(outC)
summary(outC)
prev(outC)


#######################
# Impute using Mice
#######################

library(mice)

ys <- dataSim$ys

dataSim$yo[dataSim$ys == FALSE] <- NA  
dataSim <- dataSim[, -c(1:2)]

imp <- mice(dataSim, m = 1, method = c("copulaSS", "", "", ""),  
            mice.formula = outC$mice.formula, margins = outC$margins, 
            copula = outC$BivD, maxit = 1)

comp.yo <- dataSim$yo
comp.yo[ys == 0] <- imp$imp$yo[[1]]
mean(comp.yo)

#
################
## See also ?hiv
################

############
## Example 4
############
## Generate data with partial observability

set.seed(0)

n <- 10000

Sigma <- matrix(0.5, 2, 2); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,2), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

y1 <- ifelse(-1.55 + 2*x1 + x2 + u[,1] > 0, 1, 0)
y2 <- ifelse( 0.45 - x3        + u[,2] > 0, 1, 0)
y  <- y1*y2

dataSim <- data.frame(y, x1, x2, x3)


## BIVARIATE PROBIT with Partial Observability

out  <- gjrm(list(y ~ x1 + x2, 
                  y ~ x3), 
                  data = dataSim, model = "BPO",
                  margins = c("probit", "probit"))
conv.check(out)
summary(out)

# first ten estimated probabilities for the four events from object out

cbind(out$p11, out$p10, out$p00, out$p01)[1:10,]


# case with smooth function 
# (more computationally intensive)  

f1 <- function(x) cos(pi*2*x) + sin(pi*x)

y1 <- ifelse(-1.55 + 2*x1 + f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse( 0.45 - x3            + u[,2] > 0, 1, 0)
y  <- y1*y2

dataSim <- data.frame(y, x1, x2, x3)

out  <- gjrm(list(y ~ x1 + s(x2), 
                  y ~ x3), 
                  data = dataSim, model = "BPO",
                  margins = c("probit", "probit"))

conv.check(out)
summary(out)


# plot estimated and true functions

x2 <- sort(x2); f1.x2 <- f1(x2)[order(x2)] - mean(f1(x2))
plot(out, eq = 1, scale = 0); lines(x2, f1.x2, col = "red")

#
################
## See also ?war
################
}

\dontrun{

######################################################
######################################################
######################################################
# JOINT MODELS WITH BINARY AND/OR CONTINUOUS MARGINS #
######################################################
######################################################
######################################################

library(GJRM)

############
## Example 5
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
eq.sigma1   <-    ~ 1
eq.sigma2   <-    ~ 1
eq.theta    <-    ~ x1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma1, eq.sigma2, eq.theta)

# the order above is the one to follow when
# using more than two equations

out  <- gjrm(fl, data = dataSim, margins = c("N", "N"),
             model = "B")

conv.check(out)
post.check(out)
summary(out)
AIC(out)
BIC(out)
jc.probs(out, 1.4, 2.3, intervals = TRUE)[1:4,]

############
## Example 6
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

out <- gjrm(list(y1 ~ x1 + x2, 
                 y2 ~ y1 + x2), 
                 data = dataSim, margins = c("probit","N"),
                 model = "B")
conv.check(out)                        
summary(out)
post.check(out)

## SEMIPARAMETRIC RECURSIVE Model

eq.mu.1   <- y1 ~ x1 + s(x2) 
eq.mu.2   <- y2 ~ y1 + s(x2)
eq.sigma  <-    ~ 1
eq.theta  <-    ~ 1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma, eq.theta)

out <- gjrm(fl, data = dataSim, 
            margins = c("probit","N"), gamlssfit = TRUE,
            model = "B")
conv.check(out)                        
summary(out)
post.check(out)
jc.probs(out, 1, 1.5, intervals = TRUE)[1:4,]
AT(out, nm.end = "y1")
AT(out, nm.end = "y1", type = "univariate")

#
#

############
## Example 7
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
eq.sigma  <-    ~ 1
eq.theta  <-    ~ 1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma, eq.theta)

out <- gjrm(fl, data = dataSim, 
            margins = c("probit","N"),
            model = "B")
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
## Example 8
##################
## Survival models
##################

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



eq1 <- u1 ~ s(log(u1), bs = "mpi") + z1 + s(z2)
eq2 <- u2 ~ s(log(u2), bs = "mpi") + z 
eq3 <-    ~ s(z2)

out <- gjrm(list(eq1, eq2), data = dataSim, surv = TRUE,
            margins = c("PH", "PO"), 
            cens1 = delta1, cens2 = delta2, model = "B")
                 
# PH margin fit can also be compared with cox.ph from mgcv
                 
conv.check(out)
res <- post.check(out)

## martingale residuals
mr1 <- out$cens1 - res$qr1
mr2 <- out$cens2 - res$qr2

# can be plotted against covariates
# obs index, survival time, rank order of
# surv times

# to determine func form, one may use
# res from null model against covariate

# to test for PH, use:
# library(survival)
# fit <- coxph(Surv(u1, delta1) ~ z1 + z2, data = dataSim) 
# temp <- cox.zph(fit) 
# print(temp)                  
# plot(temp, resid = FALSE)     


summary(out)
AIC(out); BIC(out)
plot(out, eq = 1, scale = 0, pages = 1)
plot(out, eq = 2, scale = 0, pages = 1)

hazsurv.plot(out, eq = 1, newdata = data.frame(z1 = 0, z2 = 0), 
             shade = TRUE, n.sim = 100, baseline = TRUE)
hazsurv.plot(out, eq = 1, newdata = data.frame(z1 = 0, z2 = 0), 
             shade = TRUE, n.sim = 100, type = "hazard", baseline = TRUE, 
             intervals = FALSE)
hazsurv.plot(out, eq = 2, newdata = data.frame(z = 0), 
             shade = FALSE, n.sim = 100, baseline = TRUE)
hazsurv.plot(out, eq = 2, newdata = data.frame(z = 0), 
             shade = TRUE, n.sim = 100, type = "hazard", baseline = TRUE)
 
jc.probs(out, type = "joint", intervals = TRUE)[1:5,]
 
newd0 <- newd1 <- data.frame(z = 0, z1 = mean(dataSim$z1), 
                             z2 = mean(dataSim$z2), 
                             u1 =  mean(dataSim$u1) + 1, 
                             u2 =  mean(dataSim$u2) + 1) 
newd1$z <- 1                   

jc.probs(out, type = "joint", newdata = newd0, intervals = TRUE)
jc.probs(out, type = "joint", newdata = newd1, intervals = TRUE)

out1 <- gjrm(list(eq1, eq2, eq3), data = dataSim, surv = TRUE,
                  margins = c("PH", "PO"), 
                  cens1 = delta1, cens2 = delta2, gamlssfit = TRUE,
                  model = "B") 


# eq1 <- u1 ~ z1 + s(z2)
# eq2 <- u2 ~ z  
# eq3 <-    ~ s(z2)    
# note that Weibull is implemented as AFT model (test case)
# out2 <- gjrm(list(eq1, eq2, ~ 1, ~ 1, eq3), data = dataSim, surv = TRUE,
#                   margins = c("WEI", "WEI"), 
#                   cens1 = delta1, cens2 = delta2,
#                   model = "B")    
                  
#########################################
## Joint continuous and survival outcomes
#########################################  
# work in progress
#
# eq1 <- z2 ~ z1
# eq2 <- u2 ~ s(u2, bs = "mpi") + z  
# eq3 <-    ~ s(z2)                  
# eq4 <-    ~ s(z2)                  
#                   
# f.l <- list(eq1, eq2, eq3, eq4)                  
#   
# out3 <- gjrm(f.l, data = dataSim, surv = TRUE,
#                   margins = c("N", "PO"), 
#                   cens1 = NULL, cens2 = delta2, 
#                   gamlssfit = TRUE, model = "B")   
# 
# conv.check(out3)
# post.check(out3)
# summary(out3)
# AIC(out3); BIC(out3)
# plot(out3, eq = 2, scale = 0, pages = 1)
# plot(out3, eq = 3, scale = 0, pages = 1)  
# plot(out3, eq = 4, scale = 0, pages = 1)                  
# 
# newd <- newd1 <- data.frame(z = 0, z1 = mean(dataSim$z1), 
#                              z2 = mean(dataSim$z2), 
#                              u2 =  mean(dataSim$u2) + 1) 
# 
# jc.probs(out3, y1 = 0.6, type = "joint", newdata = newd, intervals = TRUE)                
}

\dontrun{ 

##########################################
##########################################
##########################################
# JOINT MODELS WITH THREE BINARY MARGINS #
##########################################
##########################################
##########################################

library(GJRM)

############
## Example 9
############
## Generate data
## Correlation between the two equations 0.5 - Sample size 400 

set.seed(0)

n <- 400

Sigma <- matrix(0.5, 3, 3); diag(Sigma) <- 1
u     <- rMVN(n, rep(0,3), Sigma)

x1 <- round(runif(n)); x2 <- runif(n); x3 <- runif(n)

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2) 

y1 <- ifelse(-1.55 +    2*x1 - f1(x2) + u[,1] > 0, 1, 0)
y2 <- ifelse(-0.25 - 1.25*x1 + f2(x2) + u[,2] > 0, 1, 0)
y3 <- ifelse(-0.75 + 0.25*x1          + u[,3] > 0, 1, 0)

dataSim <- data.frame(y1, y2, y3, x1, x2)

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1)  

out  <- gjrm(f.l, data = dataSim, model = "T",
             margins = c("probit", "probit", "probit"))
out1 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T",
             margins = c("probit", "probit", "probit"))

conv.check(out)
summary(out)
plot(out, eq = 1)
plot(out, eq = 2)
AIC(out)
BIC(out)

out  <- gjrm(f.l, data = dataSim, model = "T", 
             margins = c("probit","logit","cloglog"))
out1 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T",
             margins = c("probit","logit","cloglog"))                    
conv.check(out)
summary(out)
plot(out, eq = 1)
plot(out, eq = 2)
AIC(out)
BIC(out)

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ 1, ~ 1, ~ 1) 
               
out1 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T", 
             margins = c("probit", "probit", "probit"))
   
f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ 1, ~ s(x2), ~ 1) 
               
out2 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T", 
             margins = c("probit", "probit", "probit"))   

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1, ~ s(x2), ~ x1 + s(x2)) 
               
out2 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T", 
             margins = c("probit", "probit", "probit"))   

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1, ~ x1, ~ s(x2)) 
               
out2 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T", 
             margins = c("probit", "probit", "probit")) 

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1, ~ x1 + x2, ~ s(x2)) 
               
out2 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T", 
             margins = c("probit", "probit", "probit")) 

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1 + x2, ~ x1 + x2, ~ x1 + x2) 
               
out2 <- gjrm(f.l, data = dataSim, Chol = TRUE, model = "T", 
             margins = c("probit", "probit", "probit")) 
       
       
jcres1 <- jc.probs(out2, 1, 1, 1, type = "joint", cond = 0, 
                   intervals = TRUE, n.sim = 100)
       
nw <- data.frame( x1 = 0, x2 = seq(0, 1, length.out = 100) )   
       
jcres2 <- jc.probs(out2, 1, 1, 1, newdata = nw, type = "joint", 
                   cond = 0, intervals = TRUE, n.sim = 100)
      
#############
## Example 10
#############
## Generate data
## with double sample selection

set.seed(0)

n <- 5000

Sigma <- matrix(c(1,   0.5, 0.4,
                  0.5,   1, 0.6,
                  0.4, 0.6,   1 ), 3, 3)

u <- rMVN(n, rep(0,3), Sigma)

f1   <- function(x) cos(pi*2*x) + sin(pi*x)
f2   <- function(x) x+exp(-30*(x-0.5)^2) 

x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- runif(n)
  
y1 <-  1    + 1.5*x1 -     x2 + 0.8*x3 - f1(x4) + u[, 1] > 0
y2 <-  1    - 2.5*x1 + 1.2*x2 +     x3          + u[, 2] > 0
y3 <-  1.58 + 1.5*x1 - f2(x2)                   + u[, 3] > 0

dataSim <- data.frame(y1, y2, y3, x1, x2, x3, x4)

f.l <- list(y1 ~ x1 + x2 + x3 + s(x4), 
            y2 ~ x1 + x2 + x3, 
            y3 ~ x1 + s(x2))   
          
out <- gjrm(f.l, data = dataSim, model = "TSS",
            margins = c("probit", "probit", "probit"))
conv.check(out)
summary(out)
plot(out, eq = 1)
plot(out, eq = 3)
prev(out)
prev(out, type = "univariate")
prev(out, type = "naive")
}

\dontrun{ 

###################################################
###################################################
###################################################
# JOINT MODELS WITH BINARY AND CONTINUOUS MARGINS #
# WITH SAMPLE SELECTION                           #
###################################################
###################################################
###################################################

library(GJRM)

######################################################################
## Generate data
## Correlation between the two equations and covariate correlation 0.5 
## Sample size 2000 
######################################################################
#############
## Example 11
#############
set.seed(0)

n <- 2000

rh <- 0.5      

sigmau <- matrix(c(1, rh, rh, 1), 2, 2)
u      <- rMVN(n, rep(0,2), sigmau)

sigmac <- matrix(rh, 3, 3); diag(sigmac) <- 1
cov    <- rMVN(n, rep(0,3), sigmac)
cov    <- pnorm(cov)

bi <- round(cov[,1]); x1 <- cov[,2]; x2 <- cov[,3]
  
f11 <- function(x) -0.7*(4*x + 2.5*x^2 + 0.7*sin(5*x) + cos(7.5*x))
f12 <- function(x) -0.4*( -0.3 - 1.6*x + sin(5*x))  
f21 <- function(x) 0.6*(exp(x) + sin(2.9*x)) 

ys <-  0.58 + 2.5*bi + f11(x1) + f12(x2) + u[, 1] > 0
y  <- -0.68 - 1.5*bi + f21(x1) +           u[, 2]
yo <- y*(ys > 0)
  
dataSim <- data.frame(ys, yo, bi, x1, x2)

## CLASSIC SAMPLE SELECTION MODEL
## the first equation MUST be the selection equation

resp.check(yo[ys > 0], "N")

out <- gjrm(list(ys ~ bi + x1 + x2, 
                 yo ~ bi + x1), 
                 data = dataSim, model = "BSS",
                 margins = c("probit", "N"))
conv.check(out)
post.check(out)
summary(out)

AIC(out)
BIC(out)


## SEMIPARAMETRIC SAMPLE SELECTION MODEL

out <- gjrm(list(ys ~ bi + s(x1) + s(x2), 
                 yo ~ bi + s(x1)), 
                 data = dataSim, model = "BSS",
                 margins = c("probit", "N"))
conv.check(out) 
post.check(out)
AIC(out)

## compare the two summary outputs
## the second output produces a summary of the results obtained when only 
## the outcome equation is fitted, i.e. selection bias is not accounted for

summary(out)
summary(out$gam2)

## estimated smooth function plots
## the red line is the true curve
## the blue line is the naive curve not accounting for selection bias

x1.s <- sort(x1[dataSim$ys>0])
f21.x1 <- f21(x1.s)[order(x1.s)] - mean(f21(x1.s))

plot(out, eq = 2, ylim = c(-1, 0.8)); lines(x1.s, f21.x1, col = "red")
par(new = TRUE)
plot(out$gam2, se = FALSE, lty = 3, lwd = 2, ylim = c(-1, 0.8), 
     ylab = "", rug = FALSE)


## IMPUTE MISSING VALUES

n.m <- 10
res <- imputeSS(out, n.m)
bet <- NA

for(i in 1:n.m){

dataSim$yo[dataSim$ys == 0] <- res[[i]]

outg <- gamlss(list(yo ~ bi + s(x1)), data = dataSim)
bet[i] <- coef(outg)["bi"]
print(i)
}

mean(bet)

##


## SEMIPARAMETRIC SAMPLE SELECTION MODEL with association 
## and dispersion parameters 
## depending on covariates as well

eq.mu.1   <- ys ~ bi + s(x1) + s(x2)
eq.mu.2   <- yo ~ bi + s(x1)
eq.sigma  <-    ~ bi
eq.theta  <-    ~ bi + x1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma, eq.theta)

out <- gjrm(fl, data = dataSim, model = "BSS",
                 margins = c("probit", "N"))
conv.check(out)   
post.check(out)
summary(out)
out$sigma
out$theta

jc.probs(out, 0, 0.3, intervals = TRUE)[1:4,]

outC0 <- gjrm(fl, data = dataSim, copula = "C0", model = "BSS",
              margins = c("probit", "N"))
conv.check(outC0)
post.check(outC0)
AIC(out, outC0)
BIC(out, outC0)

## IMPUTE MISSING VALUES

n.m <- 10
res <- imputeSS(outC0, n.m)

###############
# or using MICE
###############

library(mice)

ys <- dataSim$ys

dataSim$yo[dataSim$ys == FALSE] <- NA  
dataSim <- dataSim[, -1]

imp <- mice(dataSim, m = 1, method = c("copulaSS", "", "", ""),  
            mice.formula = outC0$mice.formula, margins = outC0$margins, 
            copula = outC0$BivD, maxit = 1)

comp.yo <- dataSim$yo
comp.yo[ys == 0] <- imp$imp$yo[[1]]
mean(comp.yo)


#
#
#######################################################
## example using Gumbel copula and normal-gamma margins
#######################################################
#############
## Example 12
#############
set.seed(1)

y  <- exp(-0.68 - 1.5*bi + f21(x1) + u[, 2])
yo <- y*(ys > 0)
    
dataSim <- data.frame(ys, yo, bi, x1, x2)


out <- gjrm(list(ys ~ bi + s(x1) + s(x2), 
                 yo ~ bi + s(x1)), 
                 data = dataSim, copula = "G0", 
                 margins = c("probit", "GA"),
                 model = "BSS")
conv.check(out)
post.check(out)
summary(out)


ATE <- NA
n.m <- 10
res <- imputeSS(out, n.m)

for(i in 1:n.m){

dataSim$yo[dataSim$ys == 0] <- res[[i]]

outg <- gamlss(list(yo ~ bi + s(x1)), margin = "GA", data = dataSim)

out$gamlss <- outg
ATE[i] <- AT(out, nm.end = "bi", type = "univariate")$res[2]

print(i)

}

AT(out, nm.end = "bi")
mean(ATE)
}

}

\keyword{generalised joint regression modelling}
\keyword{smooth} 
\keyword{regression spline}
\keyword{copula}
\keyword{marginal distribution}
\keyword{non-random sample selection}
\keyword{flexible copula regression modelling} 
\keyword{partial observability} 
\keyword{endogeneity} 
\keyword{survival data}







