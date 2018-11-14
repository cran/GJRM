\name{rob.const}
\alias{rob.const}
\alias{rob.const}

\title{Bootstrap procedure to help select the robust constant in a GAMLSS}

\description{ 

It helps finding the robust constant for a GAMLSS. 

}

\usage{

rob.const(x, B = 100)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{B}{Number of bootstrap replicates.}

}


\details{ 

It helps finding the robust constant for a GAMLSS based on the mean or median.

}


\value{

\item{rc}{Robust constant used in fitting.}
\item{sw}{Sum of weights for each bootstrap replicate.}
\item{m1}{Mean.}
\item{m2}{Median.}


          
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gamlss}}  
}


\keyword{robust} 
\keyword{regression} 










