\name{imputeSS}
\alias{imputeSS}

\title{Missing values' imputation}

\description{ 
\code{imputeSS} imputes missing values for a copulaSampleSel model object.

}



\usage{

imputeSS(x, m)

}

\arguments{ 

\item{x}{A fitted \code{copulaSampleSel} object as produced by the respective fitting function.}

\item{m}{Number of imputed response vectors.} 

}


\details{

This function generates m sets of imputed values for the outcome of interest under a fitted copulaSampleSel model. The 
algorithm draws parameters from the posterior distribution of copulaSampleSel which are then used to obtain simulated
responses (from the posterior predictive distribution of the missing values) via a rejection algorithm. The 
bound for acceptance/rejection is obtained via a trust region optimisation.
 
The imputed values are used to create m complete imputed datasets and perform complete data 
analysis and inference about the parameters of interest using function \code{gamlss()} within this package.

This function is currently suitable for continuous responses only.
 
}



\references{

Robert C. and  Casella G. (2004). Monte Carlo Statistical Methods. New York: Springer-Verlag.
  
Ripley B. D. (1987) Stochastic Simulation. New York: John Wiley & Sons, Inc.

}


\value{

It returns a list containing m imputed response vectors.

}




\author{

Authors: Jose Camarena, Giampiero Marra and Rosalba Radice

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}





\seealso{
\code{\link{copulaSampleSel}} 
}

\examples{
## see examples for copulaSampleSel
}

\keyword{generalised joint regression modelling} 
\keyword{non-random sample selection}
\keyword{imputation} 









