\name{res.check}
\alias{res.check}
\alias{int.rescheck}

\title{Diagnostic plots for discrete/continuous response margins}

\description{ 

It produces diagnostic plots based on (randomised) quantile residuals. 

}

\usage{

res.check(x, intervals = FALSE, n.sim = 100, prob.lev = 0.05)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{intervals}{If \code{TRUE} then intervals for the qqplots are produced.}   
  
\item{n.sim}{Number of replicate datasets used to simulate quantiles of the residual distribution.} 

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.} 

}


\details{ 

If the model fits the response well then the plots should look normally distributed.
When fitting models with discrete and/or continuous margins, four plots will be produced. In this case,
the arguments \code{main2} and \code{xlab2} come into play and allow for different
labelling across the plots. 
}


\value{

\item{qr}{It returns the (randomised) quantile residuals for the continuous or discrete margin when 
          fitting a model that involves a binary response.}
\item{qr1}{As above but for first equation (this applies when fitting models with continuous/discrete margins).}
\item{qr2}{As above but for second equation.}
          
}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gjrm}}  
}


\keyword{generalised joint regression modelling} 
\keyword{regression} 
\keyword{Q-Q plot} 
\keyword{density plot}
\keyword{histogram}









