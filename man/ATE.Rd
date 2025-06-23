\name{ATE}
\alias{ATE}

\title{Average Treatment Effect of a binary or continuous treatment variable}

\description{ 
\code{ATE} can be used to calculate the causal average treatment effect of a binary or continuous Gaussian treatment variable, with 
corresponding interval obtained using posterior simulation.}



\usage{

ATE(x, trt, trt.val = NULL, int.var = NULL, eq = NULL, joint = TRUE, n.sim = 100, 
    prob.lev = 0.05, length.out = NULL, percentage = FALSE)

}

\arguments{ 

\item{x}{A fitted \code{gjrm} object as produced by the respective fitting function.}

\item{trt}{Name of the treatment variable.} 

\item{trt.val}{Numeric value for the treatment variable. This is only required when the endogenous variable is Gaussian.} 

\item{int.var}{A vector made up of the name of the variable interacted with \code{trt}, and a value for it. It can also be a list.} 

\item{eq}{Number of equation containing the treatment variable. This is only used for trivariate models.} 


\item{joint}{If \code{FALSE} then the effect is obtained from the univariate model 
            which neglects the presence of unobserved confounders. When \code{TRUE}, the effect is obtained from 
            the simultaneous model which accounts for observed and unobserved confounders.}

\item{n.sim}{Number of simulated coefficient vectors from the posterior distribution of the estimated model 
             parameters. It may be increased if more precision is required.} 

\item{prob.lev}{Overall probability of the left and right tails of the AT distribution used for interval calculations.} 

\item{length.out}{Length of the sequence to be used when calculating the effect that a continuous
                  treatment has on a binary outcome.}

\item{percentage}{Only for the Roy model, when \code{TRUE} it provides results in terms of percentage.}





}


\details{

ATE measures the causal average difference in outcomes under 
treatment (the binary predictor or treatment assumes value 1) and under 
control (the binary treatment assumes value 0). Posterior simulation 
is used to obtain a confidence/credible interval. See the references below for details. 

ATE can also calculate the effect that a continuous Gaussian endogenous variable has on a binary outcome. In 
this case the effect will depend on the 
unit increment chosen (as shown by the plot produced). 

}



\value{

\item{res}{It returns three values: lower confidence interval limit, estimated AT and upper interval limit.}
\item{prob.lev}{Probability level used.}
\item{sim.ATE}{It returns a vector containing simulated values of the average treatment effect. This 
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








