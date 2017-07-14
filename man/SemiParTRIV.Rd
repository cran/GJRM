\name{SemiParTRIV}
\alias{SemiParTRIV}

\title{Semiparametric Trivariate Binary Models}

\description{ 
\code{SemiParTRIV} fits flexible trivariate binary models with several types of covariate effects. 

}



\usage{
SemiParTRIV(formula, data = list(), weights = NULL, subset = NULL,  
            Model = "T", margins = c("probit", "probit", "probit"),
            penCor = "unpen", sp.penCor = 3, 
            approx = FALSE, Chol = FALSE, infl.fac = 1, 
            gamma = 1, w.alasso = NULL, rinit = 1, rmax = 100, 
            iterlimsp = 50, tolsp = 1e-07,
            gc.l = FALSE, parscale, extra.regI = "t")
}


\arguments{ 

\item{formula}{In the basic setup this will be a list of three formulas. \code{s} terms 
               are used to specify smooth functions of 
                  predictors. See the examples below and the documentation of \code{SemiParBIV} for further 
                  details.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}, typically the environment from which \code{SemiParTRIV} is called.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{Model}{It indicates the type of model to be used in the analysis. Possible values are "T" (trivariate model), 
             "TSS" (trivariate model with double sample selection).}

\item{margins}{It indicates the link functions used for the three margins. Possible choices are "probit", "logit", "cloglog"".}

\item{penCor}{Type of penalty for correlation coefficients. Possible values are "unpen", "lasso", "ridge", "alasso".}

\item{sp.penCor}{Starting value for smoothing parameter of \code{penCor}.}

\item{approx}{If \code{TRUE} then an approximation of the trivariate normal integral is employed. This may speed up computations but make them 
              unstable at the same time (especially for highly correlated responses).}

\item{Chol}{If \code{TRUE} then the Cholesky method instead of the eigenvalue method is employed for the correlation matrix.}
                    
\item{infl.fac}{Inflation factor for the model degrees of freedom in the approximate AIC. Smoother models can be obtained setting 
             this parameter to a value greater than 1.}

\item{gamma}{Inflation factor used only for the alasso penalty.}

\item{w.alasso}{When using the alasso penalty a weight vector made up of three values must be provided.}

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

This function fits trivariate binary models.


For sample selection models, if there are factors in the model, before fitting, the user has to ensure 
that the numbers of factor variables' levels in the selected sample 
are the same as those in the complete dataset. Even if a model could be fitted in such a situation,
the model may produce fits which are 
not coherent with the nature of the correction sought. For more details see \code{?SemiParBIV}. 


}


\value{
The function returns an object of class \code{SemiParTRIV} as described in \code{SemiParTRIVObject}.
}


\author{

Authors: Panagiota Filippou, Giampiero Marra and Rosalba Radice

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Filippou P., Marra G. and Radice R. (in press), Penalized Likelihood Estimation of a Trivariate Additive Probit Model. \emph{Biostatistics}.


}





\section{WARNINGS}{

Convergence can be checked using \code{\link{conv.check}} which provides some 
information about 
the score and information matrix associated with the fitted model. The former should be close to 0 and the latter positive definite.
\code{SemiParTRIV()} will produce some warnings if there is a convergence issue.

Convergence failure may sometimes occur. This is not necessarily a bad thing as it may indicate specific problems 
with a fitted model. In such a situation, the user may use some extra regularisation (see \code{extra.regI}) and/or
rescaling (see \code{parscale}). Penalising the correlations using argument \code{penCor} may help a lot as
in our experience in hard situations the correlation coefficients are typically the most difficult to estimate. 
The user should also consider re-specifying/simplifying the model. Moreover, when using 
smooth functions, if the covariate's values are too sparse then convergence may be affected by this. It is also helpful to 
look into the proportions of 1 and 0 available for each event of the trivariate model; it may 
be the case that certain events do not have many observations associated with them, in which case estimation may be more challenging.

} 

\seealso{
\code{\link{SemiParBIV}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{GJRM-package}}, \code{\link{SemiParTRIVObject}}, \code{\link{conv.check}}, \code{\link{summary.SemiParTRIV}}
}

\examples{

\dontrun{ 

library(GJRM)

############
## EXAMPLE 1
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

out  <- SemiParTRIV(f.l, data = dataSim)
out1 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE)

conv.check(out)
summary(out)
plot(out, eq = 1)
plot(out, eq = 2)
AIC(out)
BIC(out)

out  <- SemiParTRIV(f.l, data = dataSim, 
                    margins = c("probit","logit","cloglog"))
out1 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE,
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
               
out1 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE)
   
f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ 1, ~ s(x2), ~ 1) 
               
out2 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE)   

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1, ~ s(x2), ~ x1 + s(x2)) 
               
out2 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE)   

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1, ~ x1, ~ s(x2)) 
               
out2 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE) 

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1, ~ x1 + x2, ~ s(x2)) 
               
out2 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE) 

f.l <- list(y1 ~ x1 + s(x2), 
            y2 ~ x1 + s(x2),
            y3 ~ x1,
               ~ x1 + x2, ~ x1 + x2, ~ x1 + x2) 
               
out2 <- SemiParTRIV(f.l, data = dataSim, Chol = TRUE) 
         
############
## EXAMPLE 2
############
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
          
out <- SemiParTRIV(f.l, data = dataSim, Model = "TSS")
conv.check(out)
summary(out)
plot(out, eq = 1)
plot(out, eq = 3)
prev(out)
prev(out, type = "univariate")
prev(out, type = "naive")

}

}

\keyword{generalised joint regression modelling}  
\keyword{smooth} 
\keyword{regression spline}








