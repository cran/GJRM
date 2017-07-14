\name{SemiParBIV}
\alias{SemiParBIV}

\title{Copula Bivariate Binary Models}

\description{ 
\code{SemiParBIV} fits copula bivariate binary models with several types of covariate 
effects, copula distributions and link functions. During the model fitting process, the 
possible presence of associated error equations, endogeneity, non-random sample selection or partial observability is accounted for.

}


\usage{
SemiParBIV(formula, data = list(), weights = NULL, subset = NULL,  
           Model = "B", BivD = "N", 
           margins = c("probit","probit"), dof = 3, gamlssfit = FALSE,
           fp = FALSE, hess = TRUE, infl.fac = 1, theta.fx = NULL, 
           rinit = 1, rmax = 100, 
           iterlimsp = 50, tolsp = 1e-07,
           gc.l = FALSE, parscale, extra.regI = "t", intf = FALSE)
}



\arguments{ 

\item{formula}{In the basic setup this will be a list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms 
               are used to specify smooth functions of 
                  predictors. \code{SemiParBIV} supports the use shrinkage smoothers for variable
                  selection purposes and more. See the examples below and the documentation of \code{mgcv} for further 
                  details on formula specifications. Note that 
                  if \code{Model = "BSS"} then the first formula MUST refer to the selection equation. Furthermore,
                  if it makes sense, a third equation for the dependence parameter can be specified (see Example 1 below).} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParBIV} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{Model}{It indicates the type of model to be used in the analysis. Possible values are "B" (bivariate model), 
             "BSS" (bivariate model with non-random sample selection), "BPO" (bivariate model with partial observability)
             and "BPO0" (bivariate model with partial observability and zero correlation).}

\item{margins}{It indicates the link functions used for the two margins. Possible choices are "probit", "logit", "cloglog".}
           
\item{dof}{If \code{BivD = "T"} then the degrees of freedom can be set to a value greater than 2 and smaller than 249.}          
           
\item{gamlssfit}{This is for internal purposes only.}               

\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "AMH", "FGM", "T", "PL", "HO" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, Ali-Mikhail-Haq,
             Farlie-Gumbel-Morgenstern, Student-t with fixed \code{dof}, Plackett, Hougaard. Note that Clayton, Joe and Gumbel are somewhat similar. Also, 
             there might be situations in which the use of a specific copula may result in more stable computations. If \code{Model = "B"} then
             each of the Clayton, Joe and Gumbel copulae is allowed to be mixed with a rotated version of the same
family. The options are: "C0C90", "C0C270", "C180C90", "C180C270", "G0G90", "G0G270", "G180G90",
"G180G270", "J0J90", "J0J270", "J180J90" and "J180J270". This allows the user to model negative and positive tail dependencies.}            

\item{fp}{If \code{TRUE} then a fully parametric model with unpenalised regression splines if fitted. See the example below.}              

\item{hess}{If \code{FALSE} then the expected/Fisher (rather than observed) information matrix is employed. The Fisher information matrix
            is not available for cases different from binary treatment and binary outcome.}

\item{infl.fac}{Inflation factor for the model degrees of freedom in the approximate AIC. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{theta.fx}{If \code{Model = "B"} and \code{BivD = "N"} then the theta parameter can be fixed in estimation.}

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
     
\item{intf}{This is for internal use.} 

}


\details{

The joint models considered by this function consist of two model equations which depend on flexible linear predictors and
whose association between the responses is modelled through parameter \eqn{\theta}{\theta} of a standardised bivariate normal 
distribution or that of a bivariate copula distribution. The linear predictors of the two equations are flexibly specified using 
parametric components and smooth functions of covariates. The same can be done for the dependence parameter if it makes sense.
Estimation is achieved within a penalized likelihood framework with integrated automatic multiple smoothing parameter selection. The use of 
penalty matrices allows for the suppression of that part of smooth term complexity which has no support 
from the data. The trade-off between smoothness 
and fitness is controlled by smoothing parameters associated with the penalty matrices. Smoothing parameters are chosen to 
minimise an approximate AIC.

For sample selection models, if there are factors in the model, before fitting, the user has to ensure 
that the numbers of factor variables' levels in the selected sample 
are the same as those in the complete dataset. Even if a model could be fitted in such a situation,
the model may produce fits which are 
not coherent with the nature of the correction sought. As an example consider the 
situation in which the complete dataset contains a factor variable with five levels and that only three of them   
appear in the selected sample. For the outcome equation (which is the one of interest) only three levels of such variable 
exist in the population, but their effects will be corrected for non-random selection using a selection equation 
in which five levels exist instead.
Having differing numbers of factors' levels between complete and selected samples will also  
make prediction not feasible (an aspect which may be particularly important for selection models);
clearly it is not possible to predict the response of interest for the missing entries using a 
dataset that contains all levels of a factor variable but using an outcome model 
estimated using a subset of these levels. 


}


\value{
The function returns an object of class \code{SemiParBIV} as described in \code{SemiParBIVObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian 
Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G., Radice R., Barnighausen T., Wood S.N. and McGovern M.E. (2017), A Simultaneous Equation Approach to Estimating HIV Prevalence with Non-Ignorable Missing Responses. \emph{Journal of the American Statistical Association}, 112(518), 484-496.

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.

Radice R., Marra G. and Wojtys M. (2016), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}, 26(5), 981-995. 

Poirier D.J. (1980), Partial Observability in Bivariate Probit Models. \emph{Journal of Econometrics}, 12, 209-217.


}



\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{SemiParBIV()} will produce some warnings when there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model. In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). These suggestions may help, especially the latter option. However, the user should also consider
re-specifying/simplifying the model and/or using a diferrent dependence structure and/or using different link functions.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size (and/or number of selected observations in selection models) is low 
compared to the complexity of the model. Examples
of misspecification include using a Clayton copula rotated by 90 degrees when a positive
association between the margins is present instead, using marginal distributions that are not adequate, and 
employing a copula which does not accommodate the type and/or strength of
the dependence between the margins (e.g., using AMH when the association between the margins is strong). When using 
smooth functions, if the covariate's values are too sparse then convergence may be affected by this.

In the contexts of endogeneity and non-random sample selection, extra attention is required when specifying
the dependence parameter as a function of covariates. This is because in these situations the dependence parameter mainly models the 
association between the unobserved confounders in the two equations. Therefore, this option would make sense when it 
is believed that the  
strength of the association between the unobservables in the two equations varies based on some grouping factor or across geographical 
areas, for instance.

} 

\seealso{
\code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{SemiParTRIV}}, \code{\link{AT}}, \code{\link{OR}}, \code{\link{RR}}, \code{\link{adjCov}}, \code{\link{prev}}, \code{\link{gt.bpm}}, \code{\link{LM.bpm}}, \code{\link{VuongClarke}}, \code{\link{plot.SemiParBIV}}, \code{\link{GJRM-package}}, \code{\link{SemiParBIVObject}}, \code{\link{conv.check}}, \code{\link{summary.SemiParBIV}}, \code{\link{predict.SemiParBIV}}
}

\examples{

library(GJRM)

############
## EXAMPLE 1
############
## Generate data
## Correlation between the two equations 0.5 - Sample size 400 

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

#
#

## CLASSIC BIVARIATE PROBIT

out  <- SemiParBIV(list(y1 ~ x1 + x2 + x3, 
                        y2 ~ x1 + x2 + x3), 
                        data = dataSim)
conv.check(out)
summary(out)
AIC(out)
BIC(out)

\dontrun{ 

## SEMIPARAMETRIC BIVARIATE PROBIT

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2
## For COPULA models use BivD argument 

out  <- SemiParBIV(list(y1 ~ x1 + s(x2, bs = "tp", k = 10, m = 2) + s(x3), 
                        y2 ~ x1 + s(x2) + s(x3)),  
                        data = dataSim)
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

## p-values suggest to drop x3 from both equations, with a stronger 
## evidence for eq. 2. This can be also achieved using shrinkage smoothers

outSS <- SemiParBIV(list(y1 ~ x1 + s(x2, bs = "ts") + s(x3, bs = "cs"), 
                         y2 ~ x1 + s(x2, bs = "cs") + s(x3, bs = "ts")), 
                         data = dataSim)
conv.check(outSS)                          

plot(outSS, eq = 1, select = 1, scale = 0, shade = TRUE)
plot(outSS, eq = 1, select = 2, ylim = c(-0.1,0.1))
plot(outSS, eq = 2, select = 1, scale = 0, shade = TRUE)
plot(outSS, eq = 2, select = 2, ylim = c(-0.1,0.1))

## SEMIPARAMETRIC BIVARIATE PROBIT with association parameter 
## depending on covariates as well

eq.mu.1  <- y1 ~ x1 + s(x2)
eq.mu.2  <- y2 ~ x1 + s(x2)
eq.theta <-    ~ x1 + s(x2)

fl <- list(eq.mu.1, eq.mu.2, eq.theta)

outD <- SemiParBIV(fl, data = dataSim)
conv.check(outD)  
summary(outD)
outD$theta

plot(outD, eq = 1, seWithMean = TRUE)
plot(outD, eq = 2, seWithMean = TRUE)
plot(outD, eq = 3, seWithMean = TRUE)
graphics.off()

#
#

############
## EXAMPLE 2
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

LM.bpm(list(y1 ~ x1 + s(x2), y2 ~ y1 + s(x2)), dataSim, Model = "B")

# p-value suggests presence of endogeneity, hence fit a bivariate model


## CLASSIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIV(list(y1 ~ x1 + x2, 
                       y2 ~ y1 + x2), 
                       data = dataSim)
conv.check(out)                        
summary(out)
AIC(out); BIC(out)

## SEMIPARAMETRIC RECURSIVE BIVARIATE PROBIT

out <- SemiParBIV(list(y1 ~ x1 + s(x2), 
                       y2 ~ y1 + s(x2)), 
                       data = dataSim)
conv.check(out)                        
summary(out)
AIC(out); BIC(out)

#

## Testing the hypothesis of absence of endogeneity post estimation... 

gt.bpm(out)

#
## reatment effect, risk ratio and odds ratio with CIs

mb(y1, y2, Model = "B")
AT(out, nm.end = "y1", hd.plot = TRUE) 
RR(out, nm.end = "y1") 
OR(out, nm.end = "y1") 
AT(out, nm.end = "y1", type = "univariate") 


## try a Clayton copula model... 

outC <- SemiParBIV(list(y1 ~ x1 + s(x2), 
                        y2 ~ y1 + s(x2)), 
                        data = dataSim, BivD = "C0")
conv.check(outC)                         
summary(outC)
AT(outC, nm.end = "y1") 

## try a Joe copula model... 

outJ <- SemiParBIV(list(y1 ~ x1 + s(x2), 
                        y2 ~ y1 + s(x2)), 
                        data = dataSim, BivD = "J0")
conv.check(outJ)
summary(outJ)
AT(outJ, "y1") 


VuongClarke(out, outJ)

#
## recursive bivariate probit modelling with unpenalized splines 
## can be achieved as follows

outFP <- SemiParBIV(list(y1 ~ x1 + s(x2, bs = "cr", k = 5), 
                         y2 ~ y1 + s(x2, bs = "cr", k = 6)), 
                         fp = TRUE, data = dataSim)
conv.check(outFP)                            
summary(outFP)

# in the above examples a third equation could be introduced
# as illustrated in Example 1

#
#################
## See also ?meps
#################

############
## EXAMPLE 3
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

## Testing the hypothesis of absence of non-random sample selection... 

LM.bpm(list(ys ~ bi + s(x1) + s(x2), yo ~ bi + s(x1)), dataSim, Model = "BSS")

# p-value suggests presence of sample selection, hence fit a bivariate model

#
## SEMIPARAMETRIC SAMPLE SELECTION BIVARIATE PROBIT
## the first equation MUST be the selection equation

out <- SemiParBIV(list(ys ~ bi + s(x1) + s(x2), 
                       yo ~ bi + s(x1)), 
                       data = dataSim, Model = "BSS")
conv.check(out)                          
gt.bpm(out)                        

## compare the two summary outputs
## the second output produces a summary of the results obtained when
## selection bias is not accounted for

summary(out)
summary(out$gam2)

## corrected predicted probability that 'yo' is equal to 1

mb(ys, yo, Model = "BSS")
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

outC <- SemiParBIV(list(ys ~ bi + s(x1) + s(x2), 
                        yo ~ bi + s(x1)), 
                        data = dataSim, Model = "BSS", BivD = "C0")
conv.check(outC)
summary(outC)
prev(outC)

# in the above examples a third equation could be introduced
# as illustrated in Example 1

#
################
## See also ?hiv
################

############
## EXAMPLE 4
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

out  <- SemiParBIV(list(y ~ x1 + x2, 
                        y ~ x3), 
                        data = dataSim, Model = "BPO")
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

out  <- SemiParBIV(list(y ~ x1 + s(x2), 
                        y ~ x3), 
                        data = dataSim, Model = "BPO")

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
}


\keyword{non-random sample selection}
\keyword{flexible copula regression modelling} 
\keyword{partial observability} 
\keyword{endogeneity} 
\keyword{smooth} 
\keyword{regression spline}








