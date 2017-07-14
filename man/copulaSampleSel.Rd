\name{copulaSampleSel}
\alias{copulaSampleSel}

\title{Copula Regression Models with Non-Random Sample Selection}

\description{ 
\code{copulaSampleSel} fits copula sample selection models with several types of covariate 
effects, copula and marginal distributions.

}



\usage{
copulaSampleSel(formula, data = list(), weights = NULL, subset = NULL,  
                 BivD = "N", margins = c("probit", "N"), dof = 3, 
                 fp = FALSE, infl.fac = 1, rinit = 1, rmax = 100, 
                 iterlimsp = 50, tolsp = 1e-07,
                 gc.l = FALSE, parscale, extra.regI = "t")
}



\arguments{ 

\item{formula}{In the basic setup, this will be a list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms 
               are used to specify smooth functions of 
                  predictors. For the case of more than two equations see the example below and the documentation 
                  of \code{SemiParBIV()} for more details. Note that 
                  the first formula MUST refer to the selection equation.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{copulaSampleSel} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}


\item{margins}{It indicates the distributions used for the two margins. The first is one of "probit", "logit", "cloglog" which refer 
               to the link function of the first equation whose response is always assumed to be binary. The response 
               for the second equation can be normal ("N"), normal where sigma2 corresponds to 
               the standard deviation instead of the variance ("N2"), log-normal ("LN"), Gumbel ("GU"), reverse Gumbel ("rGU"), 
               logistic ("LO"), Weibull ("WEI"), inverse Gaussian ("iG"), gamma ("GA"), Dagum ("DAGUM"), 
               Singh-Maddala ("SM"), beta ("BE"), Fisk ("FISK", also known as log-logistic distribution), Poisson ("PO"), zero truncated 
               Poisson ("ZTP"), negative binomial - type I ("NBI"), negative 
               binomial - type II ("NBII"), Poisson inverse Gaussian ("PIG").}
               
\item{dof}{If \code{BivD = "T"} then the degrees of freedom can be set to a value greater than 2 and smaller than 249.}          
                 
\item{BivD}{Type of bivariate error distribution employed. Possible choices are "N", "C0", "C90", "C180", "C270", "J0", "J90", "J180", "J270", 
            "G0", "G90", "G180", "G270", "F", "AMH", "FGM", "T", "PL", "HO" which stand for bivariate normal, Clayton, rotated Clayton (90 degrees), survival Clayton,
             rotated Clayton (270 degrees), Joe, rotated Joe (90 degrees), survival Joe, rotated Joe (270 degrees),
             Gumbel, rotated Gumbel (90 degrees), survival Gumbel, rotated Gumbel (270 degrees), Frank, Ali-Mikhail-Haq, 
             Farlie-Gumbel-Morgenstern, Student-t with fixed \code{dof}, Plackett, Hougaard.}            

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

This function works as \code{SemiParSampleSel()} in \code{SemiParSampleSel} and has been included in \code{GJRM} (which
already included sample selection models for binary outcomes) for the user's convenience (given several 
requests). \code{copulaSampleSel()} allows for several continuous/discrete distributions and link functions 
for the selection equation. \code{SemiParSampleSel()} allows for a probit link, and normal, gamma and several discrete distributions. 

If there are factors in the model, before fitting, the user has to ensure 
that the numbers of factor variables' levels in the selected sample 
are the same as those in the complete dataset. Even if a model could be fitted in such a situation,
the model may produce fits which are 
not coherent with the nature of the correction sought. For more details see \code{?SemiParBIV}. 

There are many continuous/discrete distributions and copula functions to choose from and we plan to include more 
options. Get in touch if you are interested in a particular distribution.

}


\value{
The function returns an object of class \code{copulaSampleSel} as described in \code{copulaSampleSelObject}.
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Marra G. and Radice R. (2013), Estimation of a Regression Spline Sample Selection Model. \emph{Computational Statistics and Data Analysis}, 61, 158-173.

Marra G. and Wyszynski K. (2016), Semi-Parametric Copula Sample Selection Models for Count Responses. \emph{Computational Statistics and Data Analysis}, 104, 110-129.

Wojtys M. and Marra G. (submitted). Copula-Based Generalized Additive Models with Non-Random Sample Selection.

}


\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{copulaSampleSel()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model. In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). However, the user should also consider
re-specifying/simplifying the model and/or using a diferrent dependence structure and/or checking that the chosen marginal 
distributions are adequate.
In our experience, we found that convergence failure typically occurs 
when the model has been misspecified and/or the sample size/number of selected observations is low compared 
to the complexity of the model. Examples
of misspecification include using a Clayton copula rotated by 90 degrees when a positive
association between the margins is present instead, using marginal distributions that do not fit
the responses or that contain many parameters, and 
employing a copula which does not accommodate the type and/or strength of
the dependence between the margins (e.g., using AMH when the association between the margins is strong). When using 
smooth functions, if the covariate's values are too sparse then convergence may be affected by this.
It is also worth bearing in mind that the use of a three parameter marginal distribution requires the data
to be more informative than a situation in which a two parameter distribution is used instead.


Extra attention is required when specifying
the dependence parameter as a function of covariates. This is because in these situations the dependence parameter mainly models the 
association between the unobserved confounders in the two equations. Therefore, this option would make sense when it 
is believed that the  
strength of the association between the unobservables in the two equations varies based on some grouping factor or across geographical 
areas, for instance.

} 

\seealso{
\code{\link{copulaReg}}, \code{\link{SemiParBIV}}, \code{\link{adjCov}}, \code{\link{VuongClarke}}, \code{\link{GJRM-package}}, \code{\link{copulaSampleSelObject}}, \code{\link{conv.check}}, \code{\link{summary.copulaSampleSel}}
}

\examples{

\dontrun{ 

library(GJRM)

######################################################################
## Generate data
## Correlation between the two equations and covariate correlation 0.5 
## Sample size 2000 
######################################################################

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

out <- copulaSampleSel(list(ys ~ bi + x1 + x2, 
                            yo ~ bi + x1), 
                       data = dataSim)
conv.check(out)
post.check(out)
summary(out)

AIC(out)
BIC(out)


## SEMIPARAMETRIC SAMPLE SELECTION MODEL

## "cr" cubic regression spline basis      - "cs" shrinkage version of "cr"
## "tp" thin plate regression spline basis - "ts" shrinkage version of "tp"
## for smooths of one variable, "cr/cs" and "tp/ts" achieve similar results 
## k is the basis dimension - default is 10
## m is the order of the penalty for the specific term - default is 2

out <- copulaSampleSel(list(ys ~ bi + s(x1, bs = "tp", k = 10, m = 2) + s(x2), 
                            yo ~ bi + s(x1)), 
                       data = dataSim)
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
eq.sigma2 <-    ~ bi
eq.theta  <-    ~ bi + x1

fl <- list(eq.mu.1, eq.mu.2, eq.sigma2, eq.theta)

out <- copulaSampleSel(fl, data = dataSim)
conv.check(out)   
post.check(out)
summary(out)
out$sigma2
out$theta

jc.probs(out, 0, 0.3, intervals = TRUE)[1:4,]

outC0 <- copulaSampleSel(fl, data = dataSim, BivD = "C0")
conv.check(outC0)
post.check(outC0)
AIC(out, outC0)
BIC(out, outC0)

## IMPUTE MISSING VALUES

n.m <- 10
res <- imputeSS(outC0, n.m)

#
#

#######################################################
## example using Gumbel copula and normal-gamma margins
#######################################################

set.seed(1)

y  <- exp(-0.68 - 1.5*bi + f21(x1) + u[, 2])
yo <- y*(ys > 0)
    
dataSim <- data.frame(ys, yo, bi, x1, x2)


out <- copulaSampleSel(list(ys ~ bi + s(x1) + s(x2), 
                            yo ~ bi + s(x1)), 
                        data = dataSim, BivD = "G0", 
                        margins = c("probit", "GA"))
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


#
#
}
}


\keyword{non-random sample selection}
\keyword{generalised joint regression modelling} 
\keyword{smooth} 
\keyword{regression spline}
\keyword{copula}
\keyword{marginal distribution}







