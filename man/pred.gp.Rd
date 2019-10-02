\name{pred.gp}
\alias{pred.gp}

\title{Function to predict quantiles from GP and DGP distributions}


\description{

It takes a fitted \code{gamlss} object produced 
by \code{gamlss()} and 
produces the desired quntities and respective intervals. 

}

\usage{
pred.gp(x, p = 0.5, newdata, n.sim = 100, prob.lev = 0.05)
}




\arguments{ 
  
\item{x}{A fitted \code{gamlss} object.}

\item{p}{Value of p.}

\item{newdata}{A data frame or list containing the values of the model covariates at which predictions are required. 
               If not provided then predictions corresponding to the original data are returned.
               When newdata is provided, it should contain all the variables needed for prediction.}
               
\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals. It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 


}






\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{gamlss}}
}




\keyword{prediction}







