\name{logLik.SemiParBIV}
\alias{logLik.SemiParBIV}

\title{Extract the log likelihood for a fitted copula model}


\description{It extracts the log-likelihood for a fitted \code{SemiParBIV}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV} model.

}

\usage{
\method{logLik}{SemiParBIV}(object, ...)
}


\arguments{ 
  
\item{object}{A fitted \code{SemiParBIV}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV} object.}

\item{...}{Un-used for this function.}


}


\details{

Modification of the classic \code{logLik} which accounts for the estimated degrees of freedom used in \code{SemiParBIV}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV} objects.
This function is provided so that information criteria work correctly by using the correct number of degrees 
of freedom. 

}



\value{

Standard \code{\link{logLik}} object.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{AIC}}, \code{\link{BIC}}
}



\examples{
## see examples for SemiParBIV
}

\keyword{logLik} 
\keyword{AIC} 
\keyword{BIC}






