\name{cond.mv}
\alias{cond.mv}

\title{Conditional Mean/Variance from a Copula Model}

\description{Function \code{cond.mv} can be used to calculate conditional means/variances from a copula model, with corresponding interval obtained using posterior simulation.}



\usage{

cond.mv(x, eq, y1 = NULL, y2 = NULL, newdata, fun = "mean", n.sim = 100, 
        prob.lev = 0.05)

}

\arguments{ 

\item{x}{A fitted \code{cond.mv} object as produced by the respective fitting function.}

\item{eq}{Equation of interest. From this, conditioning is also deduced.} 

\item{y1, y2}{Values for y1 and y2. Depending on the fitted model, one of them may be required.} 

\item{newdata}{A data frame with one row, which must be provided.}

\item{fun}{Either mean or variance.} 

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model 
             parameters. } 

\item{prob.lev}{Overall probability of the left and right tails of the simulated distribution used for interval calculations.} 

}


\details{

cond.mv() calculates the conditional mean or variance of copula models. Posterior simulation 
is used to obtain a confidence/credible interval.

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated conditional mean or variance and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.mv}{It returns a vector containing simulated values of the conditional mean or variance. This 
               is used to calculate intervals.}
        
            
}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}} 
}






\keyword{copulae} 
\keyword{bivariate model} 
\keyword{conditional mean} 
\keyword{conditional variance} 
\keyword{bayesian posterior simulation}








