\name{summary.copulaSampleSel}
\alias{summary.copulaSampleSel}
\alias{print.summary.copulaSampleSel}

\title{copulaSampleSel summary}

\description{
It takes a fitted \code{copulaSampleSel} object produced by \code{copulaSampleSel()} and produces some summaries from it.

}


\usage{

\method{summary}{copulaSampleSel}(object, n.sim = 100, prob.lev = 0.05, ...)
  
\method{print}{summary.copulaSampleSel}(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)  
  
}


\arguments{ 

\item{object}{A fitted \code{copulaSampleSel} object as produced by \code{copulaSampleSel()}.}

\item{x}{\code{summary.copulaSampleSel} object produced by \code{summary.copulaSampleSel()}.}


\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals for the association parameter, dispersion coefficient, for instance It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 
                    
\item{digits}{Number of digits printed in output.} 

\item{signif.stars}{By default significance stars are printed alongside output.}           
          
\item{...}{Other arguments.}

}


\details{ 

This function is very similar to \code{summary.SemiParBIV()}.

\code{print.summary.copulaSampleSel} prints model term summaries. 


}



\value{

Very similar to what obtained when using \code{summary.SemiParBIV()}.


}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}





\examples{
## see examples for copulaSampleSel
}


\keyword{information criteria}
\keyword{smooth}
\keyword{regression}
\keyword{non-random sample selection}






