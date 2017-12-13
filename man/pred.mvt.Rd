\name{pred.mvt}
\alias{pred.mvt}

\title{Function to predict mean and variance of marginal distributions, as well as Kendall's tau}


\description{

It takes a fitted \code{gjrm} object produced 
by \code{gjrm()} and 
produces predictions and respective intervals. 

}

\usage{
pred.mvt(x, eq, fun = "mean", n.sim = 100, prob.lev = 0.05, smooth.no = NULL, ...)
}


\arguments{ 
  
\item{x}{A fitted \code{gjrm} object.}

\item{eq}{The equation number.}

\item{fun}{Either mean, variance or tau.}

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals. It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 

\item{smooth.no}{Smooth number if the interest is in a particular smooth and not the additive predictor(s).} 


\item{...}{Other parameters.}


}






\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gjrm}}
}




\keyword{prediction}







