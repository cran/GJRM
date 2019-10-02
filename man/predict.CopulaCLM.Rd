\name{predict.CopulaCLM}
\alias{predict.CopulaCLM}


\title{Prediction function}

\description{

It takes a fitted \code{gjrm} object for the ordinal-continuous case and, 
for each equation, produces predictions 
for a new set of values of the model covariates or the original values used for the model fit. 
Standard errors of predictions can be produced and are based on the posterior distribution of the model coefficients.  
}

\usage{

\method{predict}{CopulaCLM}(object, eq, type = "link", ...)


}


\arguments{ 

\item{object}{A fitted \code{gjrm} object.}

\item{eq}{The equation to be considered for prediction.}

\item{type}{Type of prediction.}

\item{...}{Other arguments as in \code{predict.gam()} in \code{mgcv}.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}






\seealso{
\code{\link{gjrm}}
}


\keyword{prediction}







