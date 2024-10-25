\name{SATE}
\alias{SATE}

\title{Survival Average Treatment Effects of a binary treatment variable}

\description{ 
\code{SATE} can be used to calculate the survival treatment effects of a binary treatment variable, with 
corresponding interval obtained using posterior simulation.}



\usage{

SATE(x, trt, surv.t = NULL, int.var = NULL, joint = TRUE,
     n.sim = 100, prob.lev = 0.05, ls = 10, plot.type = "survival", ...)

}



\arguments{ 

\item{x}{A fitted \code{gjrm} object as produced by the respective fitting function.}

\item{trt}{Name of the treatment variable.} 

\item{surv.t}{Numeric value for time. If not provided, the function will be calculate the SATE for each time point of a grid of 
              lenght \code{ls}.} 

\item{int.var}{A vector made up of the name of the variable interacted with \code{nm.end}, and a value for it.} 

\item{joint}{If \code{FALSE} then the effects are obtained from the univariate model 
            which neglects the presence of unobserved confounders. When \code{TRUE}, the effects are obtained from 
            the simultaneous model which accounts for observed and unobserved confounders.}

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model 
             parameters. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the SAT distribution used for interval calculations.} 

\item{ls}{Length of sequence to use for time variable. Only used when \code{surv.t = NULL}.}

\item{plot.type}{Used when \code{surv.t = NULL}. Values allowed are: "survival" (for survival plots under treatment = 0 (grey lines)
                 and treatment = 1 (black lines)) and "sate" (SATE evaluated at several time points).}

\item{...}{Other graphics parameters to pass on to plotting commands.}


}


\details{

SATE measures the average survival difference in outcomes under 
treatment (the binary predictor or treatment assumes value 1) and under 
control (the binary treatment assumes value 0). Posterior simulation 
is used to obtain a confidence/credible interval.  

}



\value{

\item{res}{It returns three values: lower interval limit(s), estimated SATE(s) and upper interval limit(s).}
\item{prob.lev}{Probability level used.}
\item{sim.SATE}{It returns a vector containing simulated values of the survival average treatment effect for the case in which a specific time is chosen. This 
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
\keyword{SATE} 
\keyword{survival average treatment effect} 
\keyword{bayesian posterior simulation}








