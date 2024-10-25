\name{PE}
\alias{PE}

\title{Partial effect from a binary bivariate model}

\description{ 
\code{PE} can be used to calculate the sample treatment effect from a a binary bivariate model, with 
corresponding interval obtained using posterior simulation.}



\usage{

PE(x1, idx, n.sim = 100, prob.lev = 0.05, 
   plot = FALSE, 
   main = "Histogram of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...)

}


\arguments{ 

\item{x1}{A fitted \code{gjrm} object.}

\item{idx}{This is useful to pick a particular individual and must be provided.}

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the AT distribution used for interval calculations.} 

\item{plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated average effects is produced.}


\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}


}


\details{

PE measures the sample average effect from a binary bivariate model when a binary response (associated 
with a continuous outcome) takes values 0 and 1. Posterior simulation is used to obtain a confidence/credible interval. 

}







\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}





\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}
}





\keyword{copulae} 
\keyword{binary bivariate model} 
\keyword{PE} 
\keyword{average partial effect} 
\keyword{bayesian posterior simulation}









