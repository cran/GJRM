\name{vis.SemiParBIV}
\alias{vis.SemiParBIV}

\title{Visualization function}


\description{

It takes a fitted \code{SemiParBIV}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV}/\code{gamlss} object produced 
by \code{SemiParBIV()}, \code{copulaReg()}, \code{copulaSampleSel()}, \code{SemiParTRIV()}, \code{gamlss()} and 
produces perspective or contour plot views of model predictions. This function is a 
wrapper of \code{vis.gam()} in \code{mgcv}. Please see 
the documentation of \code{vis.gam()} for full details. 

}

\usage{
vis.SemiParBIV(x, eq, ...)
}


\arguments{ 
  
\item{x}{A fitted \code{SemiParBIV}/\code{copulaReg}/\code{copulaSampleSel}/\code{SemiParTRIV()} object.}

\item{eq}{The equation from which smooth terms should be considered for printing.}

\item{...}{Other graphics parameters to pass on to plotting commands, as described for \code{vis.gam()} in \code{mgcv}.}


}


\value{
The function generates plots.
}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{SemiParBIV}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{SemiParTRIV}}
}




\keyword{smooth} 
\keyword{regression} 
\keyword{hplot}






