\name{prev}
\alias{prev}

\title{Estimated overall prevalence from sample selection model}

\description{ 
\code{prev} can be used to calculate the overall estimated prevalence from a sample selection model 
with binay outcome, with corresponding interval
obtained using the delta method or posterior simulation.}



\usage{

prev(x, sw = NULL, type = "joint", ind = NULL, delta = FALSE,  
     n.sim = 100, prob.lev = 0.05, hd.plot = FALSE, 
     main = "Histogram and Kernel Density of Simulated Prevalences", 
     xlab = "Simulated Prevalences", ...)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{sw}{Survey weights.}  

\item{type}{This argument can take three values: \code{"naive"} (the prevalence is calculated ignoring the presence of observed 
            and unobserved confounders), \code{"univariate"} (the prevalence is obtained from the univariate probit/single imputation model 
            which neglects the presence of unobserved confounders) and \code{"joint"} (the prevalence is obtained from the 
            bivariate/trivariate model 
            which accounts for observed and unobserved confounders).}
             
\item{ind}{Binary logical variable. It can be used to calculate the prevalence for a subset of the data.} 

\item{delta}{If \code{TRUE} then the delta method is used for confidence interval calculations, otherwise Bayesian posterior 
             simulation is employed.} 
             
\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.}  
             
\item{prob.lev}{Overall probability of the left and right tails of the prevalence distribution used for interval calculations.}
              

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated prevalences is produced. This can only 
               be produced when \code{delta = FALSE}.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}

}


\details{

\code{prev} estimates the overall prevalence of a disease (e.g., HIV) when there are missing values that are not at random. 
An interval for the estimated prevalence can be obtained using the delta method or posterior simulation. 


}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated prevalence and upper confidence interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.prev}{If \code{delta = FALSE} then it returns a vector containing simulated values of the prevalence. This 
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










