\name{imputeCounter}
\alias{imputeCounter}

\title{Imputation of Counterfactual}

\description{ 
\code{imputeCounter} imputes counterfactual missing values for a gjrm model object.

}



\usage{

imputeCounter(x, m = 10, nm.end)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{m}{Number of imputed response vectors.} 

\item{nm.end}{Name endogenous variable.}

}


\details{

This function generates m sets of imputed values for the outcome of interest under a fitted joint causal model. The 
algorithm draws parameters from the posterior distribution of the model which are then used to obtain simulated
responses (from the posterior predictive distribution of the missing values) via a rejection algorithm. The 
bound for acceptance/rejection is obtained via a trust region optimisation.
 
The imputed values are used to create m complete imputed datasets and perform complete data 
analysis and inference about the parameters of interest using any summary statistics.

 
}



\references{

Robert C. and  Casella G. (2004). Monte Carlo Statistical Methods. New York: Springer-Verlag.
  
Ripley B. D. (1987) Stochastic Simulation. New York: John Wiley & Sons, Inc.

}


\value{

It returns a list containing m imputed response vectors.

}




\author{


Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}





\seealso{
\code{\link{gjrm}} 
}

\examples{
## see examples for gjrm
}

\keyword{generalised joint regression modelling} 
\keyword{endogeneity}
\keyword{imputation} 









