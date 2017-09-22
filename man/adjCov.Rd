\name{adjCov}
\alias{adjCov}
\alias{aCov}

\title{Adjustment for the covariance matrix from a fitted gjrm model}

\description{ 

\code{adjCov} can be used to adjust the covariance matrix of a fitted \code{gjrm} object.}


\usage{

adjCov(x, id)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object as produced by the respective fitting function.}
 
\item{id}{Cluster identifier.} 




}


\details{

This adjustment can be made when dealing with clustered data and the cluster structure is neglected when fitting the model. The basic idea is that 
the model is fitted as though observations were independent, and subsequently adjust the covariance matrix of the parameter estimates. Using the
terminology of Liang and Zeger (1986), this would correspond to using an independence structure within 
the context of generalized estimating equations. The parameter estimators are still consistent but are inefficient as compared
to a model which accounts for the correct cluster dependence structure. The covariance matrix of the independence estimators can be adjusted 
as described in Liang and Zeger (1986, Section 2). 

}



\value{

This function returns a fitted object which is identical to that supplied in \code{adjCov} but with adjusted covariance matrix.


}



\section{WARNINGS}{

This correction may not be appropriate for models fitted using penalties.

} 


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\references{

Liang K.-Y. and Zeger S. (1986), Longitudinal Data Analysis Using Generalized Linear Models. \emph{Biometrika}, 73(1), 13-22.

}




\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}
}


\keyword{covariance matrix adjustment} 

