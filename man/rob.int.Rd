\name{rob.int}
\alias{rob.int}
\alias{rob.int}

\title{Tool for tuning bounds of integral in robust models}

\description{ 

Tool for tuning bounds of integral in robust GAMLSS. 

}

\usage{

rob.int(x, rc, l.grid = 1000, tol = 1e-4, var.range = NULL)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object, typically from a non-robust fit.}

\item{rc}{Robust tuning constant.}

\item{l.grid}{Length of grid.}

\item{tol}{Tolerance}

\item{var.range}{Range of values, min and max, to use in calculations.}


}


\details{ 

Tool for tuning bounds of integral in robust GAMLSS. 

}


\value{

\item{lB, uB}{Lower and upper bounds.}

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gamlss}}  
}


\keyword{robust} 
\keyword{regression} 










