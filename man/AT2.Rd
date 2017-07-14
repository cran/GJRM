\name{AT2}
\alias{AT2}

\title{Average treatment effect from a two-part model}

\description{ 
\code{AT2} can be used to calculate the sample average treatment effect from a two-part model, with 
corresponding interval obtained using posterior simulation.}



\usage{

AT2(x1, x2, index1, index2, n.sim = 100, prob.lev = 0.05, 
    hd.plot = FALSE, 
    main = "Histogram and Kernel Density of Simulated Average Effects", 
    xlab = "Simulated Average Effects", ...)

}


\arguments{ 

\item{x1}{A fitted \code{SemiParBIV} object as produced by \code{SemiParBIV()}.}

\item{x2}{A fitted \code{SemiParBIV} object as produced by \code{SemiParBIV()}.}

\item{index1}{This is useful to pick a particular individual.}

\item{index2}{As above.}

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the AT distribution used for interval calculations.} 

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated average effects is produced. This can 
               only be produced when \code{delta = FALSE}.}


\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}


}


\details{

AT measures the sample average effect from a two-part model when a binary response (associated 
with a continuous outcome) takes values 0 and 1. Posterior simulation is used to obtain a confidence/credible interval. 

}




\section{WARNINGS}{

This function is only suitable for \code{SemiParBIV()}.

} 






\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}





\seealso{
\code{\link{GJRM-package}}, \code{\link{SemiParBIV}}, \code{\link{summary.SemiParBIV}}  
}

\examples{
## see examples for SemiParBIV
}

\keyword{copulae} 
\keyword{bivariate binary model} 
\keyword{ATE} 
\keyword{average treatment effect} 
\keyword{bayesian posterior simulation}









