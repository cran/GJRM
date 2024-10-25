\name{marg.mv}
\alias{marg.mv}

\title{Marginal Mean/Variance}

\description{Function \code{marg.mv} can be used to calculate marginal means/variances, with corresponding interval obtained using posterior simulation.}



\usage{

marg.mv(x, eq, newdata, fun = "mean", n.sim = 100, prob.lev = 0.05, bin.model = NULL)

}

\arguments{ 

\item{x}{A fitted \code{marg.mv} object as produced by the respective fitting function.}

\item{eq}{Number of equation of interest.} 

\item{newdata}{A data frame with one row, which must be provided.}

\item{fun}{Either mean or variance.}

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model 
             parameters. } 

\item{prob.lev}{Overall probability of the left and right tails of the simulated distribution used for interval calculations.} 

\item{bin.model}{If a two part or hurdle model is used then this is the object of a binary regression model fitted using gam() from mgcv.} 


}


\details{

marg.mv() calculates the marginal mean or variance. Posterior simulation 
is used to obtain a confidence/credible interval.

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated marginal mean or variance and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.mv}{It returns a vector containing simulated values of the marginal mean or variance. This 
               is used to calculate intervals.}
        
            
}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}} 
}



\keyword{bivariate model} 
\keyword{marginal mean} 
\keyword{marginal variance} 
\keyword{bayesian posterior simulation}








