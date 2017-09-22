\name{predict.SemiParBIV}
\alias{predict.SemiParBIV}


\title{Prediction function}

\description{

It takes a fitted \code{gjrm} object and, 
for each equation, produces predictions 
for a new set of values of the model covariates or the original values used for the model fit. 
Standard errors of predictions can be produced and are based on the posterior distribution of the model coefficients. This function is a 
wrapper for \code{predict.gam()} in \code{mgcv}. Please see the documentation of \code{predict.gam()} for full details. 
}

\usage{

\method{predict}{SemiParBIV}(object, eq, ...)


}


\arguments{ 

\item{object}{A fitted \code{gjrm} object.}

\item{eq}{The equation to be considered for prediction.}

\item{...}{Other arguments as in \code{predict.gam()} in \code{mgcv}.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\section{WARNINGS}{

When \code{type = "response"} (which gives predictions on the scale of the response variable), for 
the case of continuous responses this function will NOT produce correct predictions for the outcome variable (except for some the Gaussian case). This
is because for all distributions (except the Gaussian) implemented in this package the distributional parameters determine the mean and variance
through functions of them.

When predicting based on a new data set, this function can not return correct predictions for models based 
on a BivD value of "C0C90", "C0C270", "C180C90", "C180C270", "G0G90", "G0G270", "G180G90",
"G180G270", "J0J90", "J0J270", "J180J90" or "J180J270".


} 






\seealso{
\code{\link{gjrm}}
}


\keyword{prediction}







