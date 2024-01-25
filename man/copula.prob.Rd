\name{copula.prob}
\alias{copula.prob}
\alias{jc.probs1}
\alias{jc.probs2}
\alias{jc.probs3}
\alias{jc.probs4}
\alias{jc.probs5}
\alias{jc.probs6}
\alias{jc.probs7}
\alias{jc.probs8}



\title{Copula probabilities (joint and conditional) from a fitted simultaneous model}

\description{ 
\code{copula.prob} can be used to calculate the joint or conditional copula probabilities from a fitted simultaneous model with intervals obtained 
via posterior simulation.}



\usage{

copula.prob(x, y1, y2, y3 = NULL, newdata, joint = TRUE, cond = 0,
            type = "surv", intervals = FALSE, n.sim = 100, prob.lev = 0.05, 
            min.pr = 1e-323, max.pr = 1, cumul = "no")

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object as 
         produced by the respective fitting function.}

\item{y1}{Value of response for first margin.}

\item{y2}{Value of response for second margin.}

\item{y3}{Value of response for third margin if a trivariate model is employed.}

\item{newdata}{A data frame or list containing the values of the model covariates at which predictions are required. 
               If not provided then predictions corresponding to the original data are returned.
               When newdata is provided, it should contain all the variables needed for prediction.}

\item{joint}{If \code{TRUE} then the calculation is done using the fitted joint model. If \code{FALSE} then 
             the calculation is done from univariate fits.}

\item{cond}{There are three possible values: 0 (joint probabilities are delivered), 1 (conditional probabilities are delivered and 
            conditioning is with the respect to the first margin), 2 (as before but conditioning is with the respect to 
            the second margin).}  
  
\item{type}{This argument is only valid for survival copula models. It can take values: "surv", "hazard", "cumhaz".}  
  
\item{intervals}{If \code{TRUE} then intervals for the probabilities are also produced.}   
  
\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used
             for interval calculations.} 

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.} 
  

 \item{min.pr, max.pr}{Allowed minimum and maximum for estimated probabities.}


\item{cumul}{Only used for discrete and continuous margins' case.} 


}


\details{

This function calculates joint or conditional copula probabilities from a fitted simultaneous model or a model assuming independence, with intervals 
obtained via posterior simulation.

}



\value{

\item{res}{It returns several values including: estimated probabilities (\code{p12}), with lower and upper interval limits (\code{CIpr}) 
if \code{intervals = TRUE}, and \code{p1}, \code{p2} and \code{p3} (the marginal probabilities).}
                        
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}  
}




\keyword{generalised joint regression modelling} 
\keyword{copula}
\keyword{marginal distribution}
\keyword{bayesian posterior simulation}





