\name{vis.gjrm}
\alias{vis.gjrm}

\title{Visualization function}


\description{

It takes a fitted \code{gjrm} object produced 
by \code{gjrm()} and 
produces perspective or contour plot views of model predictions. This function is a 
wrapper of \code{vis.gam()} in \code{mgcv}. Please see 
the documentation of \code{vis.gam()} for full details. 

}

\usage{
vis.gjrm(x, eq, fun = NULL, ...)
}


\arguments{ 
  
\item{x}{A fitted \code{gjrm} object.}

\item{eq}{The equation number.}

\item{fun}{Either mean or variance. If left as equal to NULL then predictions on the scale of the predictor will be produced.}

\item{...}{Other graphics parameters to pass on to plotting commands, as described for \code{vis.gam()} in \code{mgcv}.}


}


\value{
The function generates plots.
}








\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gjrm}}
}




\keyword{smooth} 
\keyword{regression} 
\keyword{hplot}






