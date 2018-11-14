\name{cv.inform}
\alias{cv.inform}

\title{Cross validation for informative censoring univariate survival models}

\description{ 
\code{cv.inform} carries out cross validation to help choosing the set of informative covariates.}



\usage{

cv.inform(x, K = 5, data, informative = "yes")

}

\arguments{ 

\item{x}{A fitted \code{gamlss} object as produced by the respective fitting function.}

\item{K}{No. of folds.} 

\item{data}{Data.} 

\item{informative}{If no then cv is carried out for the case of no informative censoring. This is useful for comparison purposes.} 



}


\details{

\code{cv.inform} carries out cross validation to help choosing the set of informative covariates.

}



\value{

\item{sl}{Overall sum of predicted likelihood contributions.}
          
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\seealso{
\code{\link{GJRM-package}}, \code{\link{gamlss}} 
}


\keyword{regression modelling} 
\keyword{survival data}








