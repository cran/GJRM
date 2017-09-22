\name{AT}
\alias{AT}

\title{Average treatment effect of a binary/continuous/discrete endogenous variable}

\description{ 
\code{AT} can be used to calculate the treatment effect of a binary/continuous/discrete endogenous predictor/treatment, with 
corresponding interval obtained using posterior simulation.}



\usage{

AT(x, nm.end, eq = NULL, E = TRUE, treat = TRUE, type = "joint", ind = NULL, 
   n.sim = 100, prob.lev = 0.05, length.out = NULL,
   hd.plot = FALSE, te.plot = FALSE, 
   main = "Histogram and Kernel Density of Simulated Average Effects", 
   xlab = "Simulated Average Effects", ...)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object as produced by the respective fitting function.}

\item{nm.end}{Name of the endogenous variable.} 

\item{eq}{Number of equation containing the endogenous variable. This is only used for trivariate models.} 

\item{E}{If \code{TRUE} then \code{AT} calculates the sample ATE. If \code{FALSE} then it calculates the sample AT 
         for the treated individuals only.}  

\item{treat}{If \code{TRUE} then \code{AT} calculates the AT using the treated only. If \code{FALSE} then it calculates the effect on 
the control group. This only makes sense if \code{E = FALSE}.}  

\item{type}{This argument can take three values: \code{"naive"} (the effect is calculated ignoring the presence of observed and unobserved 
            confounders), \code{"univariate"} (the effect is obtained from the univariate model 
            which neglects the presence of unobserved confounders) and \code{"joint"} (the effect is obtained from 
            the simultaneous model which accounts for observed and unobserved confounders).}



\item{ind}{Binary logical variable. It can be used to calculate the AT for a subset of the data. Note that it does not make sense to 
           use \code{ind} when some observations are excluded from the AT calculation (e.g., when using \code{E = FALSE}).}  

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
             when \code{delta = FALSE}. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the AT distribution used for interval calculations.} 

\item{length.out}{Ddesired length of the sequence to be used when calculating the effect that a continuous/discrete
                  treatment has on a binary outcome.}



\item{hd.plot}{If \code{TRUE} then a plot of the histogram and kernel density estimate of the simulated average effects is produced. This can 
               only be produced when when binary responses are used.}

\item{te.plot}{For the case of continuous/discrete endogenous variable and binary outcome, if \code{TRUE} then a plot  
                 showing the treatment effects that the binary outcome is equal to 1 for each incremental value of the endogenous variable 
                 and respective intervals is produced.}

\item{main}{Title for the plot.}

\item{xlab}{Title for the x axis.}

\item{...}{Other graphics parameters to pass on to plotting commands. These are used only when \code{hd.plot = TRUE}.}


}


\details{

AT measures the average difference in outcomes under 
treatment (the binary predictor or treatment assumes value 1) and under 
control (the binary treatment assumes value 0). Posterior simulation 
is used to obtain a confidence/credible interval. See the references below for details. 

AT can also calculate the effect that a continuous/discrete endogenous variable has on a binary outcome. In 
this case the effect will depend on the 
unit increment chosen (as shown by the plot produced). 

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated AT and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.AT}{It returns a vector containing simulated values of the average treatment effect. This 
               is used to calculate intervals.}
\item{Effects}{For the case of continuous/discrete endogenous variable and binary outcome, it returns a matrix made up of 
               three columns containing the effects for each incremental value in the endogenous variable and respective intervals.}            
           
               

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

}


\seealso{
\code{\link{GJRM-package}}, \code{\link{gjrm}} 
}






\keyword{copulae} 
\keyword{bivariate model} 
\keyword{trivariate model} 
\keyword{ATE} 
\keyword{average treatment effect} 
\keyword{bayesian posterior simulation}








