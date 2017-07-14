\name{adjCovSD}
\alias{adjCovSD}

\title{Adjustment for the covariance matrix from a SemiParBIV/copulaReg/copulaSampleSel/SemiParTRIV model fitted to complex survey data.}

\description{ 

\code{adjCovSD} can be used to adjust the covariance matrix of a fitted \code{SemiParBIV}, \code{copulaReg}, \code{copulaSampleSel}, \code{SemiParTRIV} object.}


\usage{

adjCovSD(x, design)

}

\arguments{ 

\item{x}{A fitted \\code{SemiParBIV}, \code{copulaReg}, \code{copulaSampleSel}, \code{SemiParTRIV} object as produced by the respective fitting function.}
 
\item{design}{A \code{svydesign} object as produced by \code{svydesign()} from the \code{survey} package.} 




}


\details{

This function has been extracted from the \code{survey} package and adapted to the class of this package's models. It computes the sandwich 
variance estimator for a copula model fitted to data from a complex sample survey (Lumley, 2004).






}



\value{


This function returns a fitted object which is identical to that supplied in \code{adjCovSD} but with adjusted covariance matrix.



}



\section{WARNINGS}{

This correction may not be appropriate for models fitted using penalties.

} 


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\references{

Lumley T. (2004), Analysis of Complex Survey Samples. \emph{Journal of Statistical Software}, 9(8), 1-19.


}




\seealso{
\code{\link{GJRM-package}}, \code{\link{SemiParBIV}}, \code{\link{summary.SemiParBIV}},
\code{\link{SemiParTRIV}}, \code{\link{summary.SemiParTRIV}},
\code{\link{copulaReg}}, \code{\link{summary.copulaReg}},
\code{\link{copulaSampleSel}}, \code{\link{summary.copulaSampleSel}}
}


\keyword{covariance matrix adjustment} 
\keyword{complex survey design} 










