\name{prev}
\alias{prev}

\title{Estimated overall prevalence from a sample selection model}

\description{ 
\code{prev} can be used to calculate the overall estimated prevalence from a sample selection model 
with binay outcome, with corresponding interval
obtained using posterior simulation.}



\usage{

prev(x, sw = NULL, joint = TRUE, n.sim = 100, prob.lev = 0.05)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{sw}{Survey weights.}  

\item{joint}{If \code{FALSE} then the prevalence is obtained from the univariate model 
            which neglects the presence of unobserved confounders. When \code{TRUE}, the prevalence is obtained from 
            the simultaneous model which accounts for observed and unobserved confounders.}
                       
\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. It may 
             be increased if more precision is required.}  
             
\item{prob.lev}{Overall probability of the left and right tails of the prevalence distribution used for interval calculations.}
              

}


\details{

\code{prev} estimates the overall prevalence of a disease (e.g., HIV) when there are missing values that are not at random. 
An interval for the estimated prevalence can be obtained using posterior simulation. 


}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated prevalence and upper confidence interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.prev}{Vector containing simulated values of the prevalence. This 
               is used to calculate an interval.}

}




\author{


Authors: Giampiero Marra, Rosalba Radice, Guy Harling, Mark E McGovern 

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}


}

\references{

Marra G., Radice R., Barnighausen T., Wood S.N. and McGovern M.E. (2017), A Simultaneous Equation Approach to Estimating HIV Prevalence with Non-Ignorable Missing Responses. \emph{Journal of the American Statistical Association}, 112(518), 484-496.



}




\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}}  
}



\keyword{prevalence} 
\keyword{non-random sample selection} 
\keyword{generalised joint regression modelling} 
\keyword{bayesian posterior simulation}










