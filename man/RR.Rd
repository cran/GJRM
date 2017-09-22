\name{RR}
\alias{RR}

\title{Causal risk ratio of a binary/continuous/discrete endogenous variable}

\description{ 
\code{RR} can be used to calculate the causal risk ratio of a binary/continuous/discrete endogenous predictor/treatment, with 
corresponding interval obtained using posterior simulation.}



\usage{

RR(x, nm.end, E = TRUE, treat = TRUE, type = "joint", ind = NULL, 
   n.sim = 100, prob.lev = 0.05, length.out = NULL, hd.plot = FALSE,
   rr.plot = FALSE, 
   main = "Histogram and Kernel Density of Simulated Risk Ratios", 
   xlab = "Simulated Risk Ratios", ...)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object.}

\item{nm.end}{Name of the endogenous variable.} 

\item{E}{If \code{TRUE} then \code{RR} calculates the sample RR. If \code{FALSE} then it calculates the sample RR 
         for the treated individuals only.}  

\item{treat}{If \code{TRUE} then \code{RR} calculates the RR using the treated only. If \code{FALSE} then it calculates the ratio using  
the control group. This only makes sense if \code{E = FALSE}.}  

\item{type}{This argument can take three values: \code{"naive"} (the effect is calculated ignoring the presence of observed and 
            unobserved confounders), \code{"univariate"} (the effect is obtained 
            from the univariate model which neglects the presence of unobserved confounders) and \code{"joint"} (the effect is obtained from the bivariate model 
            which accounts for observed and unobserved confounders).}

\item{ind}{Binary logical variable. It can be used to calculate the RR for a subset of the data. Note that it does not make sense to use \code{ind} 
           when some observations are excluded from the RR calculation (e.g., when using \code{E = FALSE}).}  

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the RR distribution used for interval calculations.} 

\item{length.out}{Ddesired length of the sequence to be used when calculating the effect that a continuous/discrete
                  treatment has on a binary outcome.}

\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated risk ratios is produced. This can 
               only be produced when binary responses are used.}

\item{rr.plot}{For the case of continuous/discrete endogenous variable and binary outcome, if \code{TRUE} then a plot (on the log scale)  
                 showing the risk ratios that the binary outcome is equal to 1 for each incremental value of the endogenous variable 
                 and respective intervals is produced.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}


}


\details{

RR calculates the causal risk ratio of the probabilities of positive outcome under 
treatment (the binary predictor or treatment assumes value 1) and under 
control (the binary treatment assumes value 0). Posterior simulation is used to obtain a 
confidence/credible interval. 

RR works also for the case of continuous/discrete endogenous treatment variable.


}



\value{

\item{prob.lev}{Probability level used.}
\item{sim.RR}{It returns a vector containing simulated values of the average RR. This 
               is used to calculate intervals.}           
\item{Ratios}{For the case of continuous/discrete endogenous variable and binary outcome, it returns a matrix made up of 
               three columns containing the risk ratios for each incremental value in the endogenous variable and respective intervals.}            
            
}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}} 
}




\keyword{generalised joint regression modelling} 
\keyword{RR} 
\keyword{risk ratio} 
\keyword{bayesian posterior simulation}








