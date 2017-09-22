\name{summary.SemiParTRIV}
\alias{summary.SemiParTRIV}
\alias{print.summary.SemiParTRIV}

\title{SemiParTRIV summary}

\description{
It takes a fitted \code{SemiParTRIV} object and produces some summaries from it.

}


\usage{

\method{summary}{SemiParTRIV}(object, n.sim = 100, prob.lev = 0.05, ...)

\method{print}{summary.SemiParTRIV}(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)                      
}


\arguments{ 

\item{object}{A fitted \code{SemiParTRIV} object.}

\item{x}{\code{summary.SemiParTRIV} object produced by \code{summary.SemiParTRIV()}.}


\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals for the association parameter and other measures. It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 
            
\item{digits}{Number of digits printed in output.} 

\item{signif.stars}{By default significance stars are printed alongside output.} 

\item{...}{Other arguments.}            
            

}


\details{ 


\code{print.summary.SemiParTRIV} prints model term summaries. 


}





\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\examples{
## see examples for gjrm
}


\keyword{information criteria}
\keyword{smooth}
\keyword{regression}







