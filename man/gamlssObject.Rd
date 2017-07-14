\name{gamlssObject}
\alias{gamlssObject}

\title{Fitted gamlssObject object}

\description{A fitted gamlss object returned by function \code{gamlss} and of class "gamlss" and "SemiParBIV".}


\value{ 


\item{fit}{List of values and diagnostics extracted from the output of the algorithm.}

\item{gam1, gam2, gam3}{Univariate starting values' fits.}

\item{coefficients}{The coefficients of the fitted model.}

\item{weights}{Prior weights used during model fitting.}

\item{sp}{Estimated smoothing parameters of the smooth components.}

\item{iter.sp}{Number of iterations performed for the smoothing parameter estimation step.}

\item{iter.if}{Number of iterations performed in the initial step of the algorithm.}

\item{iter.inner}{Number of iterations performed within the smoothing parameter estimation step.}

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

\item{edf1, edf2, edf3, ...}{Degrees of freedom for the model's equations.}
                
\item{wor.c}{Working model quantities.}                
                  
\item{eta1, eta2, eta3, ...}{Estimated linear predictors.}

\item{y1}{Response.}

\item{logLik}{Value of the (unpenalized) log-likelihood evaluated at the (penalized or unpenalized) parameter 
            estimates.}


}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{gamlss}}, \code{\link{summary.gamlss}}
}











