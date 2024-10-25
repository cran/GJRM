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
            intervals = FALSE, n.sim = 100, prob.lev = 0.05, 
            theta = FALSE, tau = FALSE, min.pr = 1e-323, max.pr = 1)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object as 
         produced by the respective fitting function.}

\item{y1}{Value of response for first margin.}

\item{y2}{Value of response for second margin.}

\item{y3}{Value of response for third margin if a trivariate model is employed.}

\item{newdata}{A data frame with one row, which must be provided.}

\item{joint}{If \code{TRUE} then the calculation is done using the fitted joint model. If \code{FALSE} then 
             the calculation is done from univariate fits.}

\item{cond}{There are three possible values: 0 (joint probabilities are delivered), 1 (conditional probabilities are delivered and 
            conditioning is with the respect to the first margin), 2 (as before but conditioning is with the respect to 
            the second margin).}  
    
\item{intervals}{If \code{TRUE} then intervals for the probabilities are also produced.}   
  
\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used
             for interval calculations.} 

\item{prob.lev}{Overall probability of the left and right tails of the probabilities' distributions used for interval calculations.} 
  
\item{theta}{If \code{TRUE} the theta dependence parameter will be shown. This is especially useful for prediction purposes when theta is specified as 
             a function of covariate effects.}   
  
\item{tau}{If \code{TRUE} the Kendall's tau will also be calculated and provided in output. Note that the calculation adopted here assumes continuous 
           margins. In all other cases, this may provide a rough indication of dependence under certain assumptions. Note that, for the F, PL and J0 (and the related
           rotations), computing times may be longer than for the other cases. This is especially useful for prediction purposes when theta is specified as 
             a function of covariate effects, with an interest in analysing a more interpretable measure of dependence for certain copulae.} 

 \item{min.pr, max.pr}{Allowed minimum and maximum for estimated probabities.}



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





