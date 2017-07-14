\name{SemiParTRIVObject}
\alias{SemiParTRIVObject}

\title{Fitted SemiParTRIV object}

\description{A fitted semiparametric trivariate probit object returned by function \code{SemiParTRIV} and of class "SemiParTRIV".}


\value{ 


\item{fit}{List of values and diagnostics extracted from the output of the algorithm. For instance, \code{fit$gradient}, \code{fit$Fisher} and \code{fit$S.h} return
           the gradient vector, Fisher information (when used) and overall penalty matrix scaled 
           by its smoothing parameters, for the fitted bivariate probit model. See the 
           documentation of \code{trust} for details on the diagnostics provided.}

\item{gam1}{Univariate fit for equation 1. See the documentation of \code{mgcv} for full details.}

\item{gam2, gam3, ...}{Univariate fit for equation 2 and equations 3 and 4 (these are available when the dispersion and association parameters 
                       are modelled 
                       as functions of covariates).}

\item{coefficients}{The coefficients of the fitted model. They are given in the following 
                    order: parametric and regression spline (if present)
                    coefficients for the first equation, parametric and regression spline
                    coefficients for the second equation, and dispersion parameter (or coefficients for the third equation) and 
                    association coefficient (or coefficients for the fourth equation).}

\item{weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.if}{Number of iterations performed in the initial step of the algorithm.}

\item{iter.inner}{Number of iterations performed within the smoothing parameter estimation step.}

\item{theta}{Estimated dependence parameter linking the two equations.}

\item{n}{Sample size.}

\item{X1, X2, X3, ...}{Design matrices associated with the linear predictors.}

\item{X1.d2, X2.d2, X3.d2, ...}{Number of columns of \code{X1}, \code{X2}, \code{X3}, etc.}

\item{l.sp1, l.sp2, l.sp3, ...}{Number of smooth components in the equations.}

\item{He}{Penalized -hessian/Fisher. This is the same as \code{HeSh} for unpenalized models.}

\item{HeSh}{Unpenalized -hessian/Fisher.}

\item{Vb}{Inverse of \code{He}. This corresponds to the Bayesian variance-covariance matrix 
          used for confidence/credible interval calculations.}
 
\item{F}{This is obtained multiplying Vb by HeSh.}  
 
\item{t.edf}{Total degrees of freedom of the estimated bivariate model. It is calculated as \code{sum(diag(F))}.}

\item{edf1, edf2, edf3, ...}{Degrees of freedom for the two equations of the fitted bivariate model (and for the third and fourth
                        equations if present. They  
                        are calculated when splines are used.}

\item{bs.mgfit}{List of values and diagnostics extracted from \code{magic} in \code{mgcv}.}

\item{conv.sp}{If \code{TRUE} then the smoothing parameter selection algorithm stopped before reaching the maximum number of iterations allowed.}
                
\item{wor.c}{Working model quantities.}                

\item{p111, p110,...}{Model probabilities.}
                       
\item{eta1, eta2, eta3, ...}{Estimated linear predictors for the two equations (as well as the third and fourth equations if present).}

\item{y1, y2, y3}{Responses of the two equations.}

\item{logLik}{Value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}

\item{respvec}{List containing response vectors.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{SemiParTRIV}}, \code{\link{summary.SemiParTRIV}}
}











