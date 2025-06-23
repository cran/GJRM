\name{OR}
\alias{OR}

\title{Causal odds ratio of a binary/continuous treatment variable}

\description{ 
\code{OR} can be used to calculate the causal odds ratio of a binary/continuous treatment variable, with 
corresponding interval obtained using posterior simulation.}



\usage{

OR(x, trt, trt.val = NULL, int.var = NULL, joint = TRUE, n.sim = 100, prob.lev = 0.05, 
   length.out = NULL)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}
 
\item{trt}{Name of the treatment variable.}  

\item{trt.val}{Numeric value for the treatment variable. This is only required when the endogenous variable is Gaussian.} 


\item{int.var}{A vector made up of the name of the variable interacted with \code{trt}, and a value for it. It can also be a list.} 

\item{joint}{If \code{FALSE} then the effect is obtained from the univariate model 
            which neglects the presence of unobserved confounders. When \code{TRUE}, the effect is obtained from 
            the simultaneous model which accounts for observed and unobserved confounders.}

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the OR distribution used for interval calculations.} 

\item{length.out}{Ddesired length of the sequence to be used when calculating the effect that a continuous
                  treatment has on a binary outcome.}

}


\details{

OR calculates the causal odds ratio for a binary/continuous Gaussian treatment. Posterior simulation is used to obtain a 
confidence/credible interval. 


}



\value{

\item{prob.lev}{Probability level used.}
\item{sim.OR}{It returns a vector containing simulated values of the average OR. This 
               is used to calculate intervals.}           
\item{Ratios}{For the case of continuous endogenous treatment and binary outcome, it returns a matrix made up of 
               three columns containing the odds ratios for each incremental value in the endogenous variable and respective intervals.}            
          
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}  
}




\keyword{generalised joint regression modelling} 
\keyword{OR} 
\keyword{odds ratio} 
\keyword{bayesian posterior simulation}








