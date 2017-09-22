\name{post.check}
\alias{post.check}
\alias{int.postcheck}

\title{Diagnostic plots for discrete/continuous response margin}

\description{ 

It produces diagnostic plots based on (randomised) quantile residuals. 

}

\usage{

post.check(x, main = "Histogram and Density Estimate of Residuals", 
           main2 = "Histogram and Density Estimate of Residuals",
           xlab = "Quantile Residuals", xlab2 = "Quantile Residuals", 
           intervals = FALSE, n.sim = 100, prob.lev = 0.05, ...)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{main}{Title for the plot.}

\item{main2}{Title for the plot in the second row. This comes into play only when fitting models with two non-binary margins.}

\item{xlab}{Title for the x axis.}

\item{xlab2}{Title for the x axis in the second row. As above.}

\item{intervals}{If \code{TRUE} then intervals for the qqplots are produced.}   
  
\item{n.sim}{Number of replicate datasets used to simulate quantiles of the residual distribution.} 

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.} 

\item{...}{Other graphics parameters to pass on to plotting commands.}

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









