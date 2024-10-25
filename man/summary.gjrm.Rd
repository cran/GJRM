\name{summary.gjrm}
\alias{summary.gjrm}
\alias{print.summary.gjrm}

\title{gjrm summary}

\description{
It takes a fitted \code{gjrm} object and produces some summaries from it.

}


\usage{

\method{summary}{gjrm}(object, n.sim = 100, prob.lev = 0.05, ...)
   


\method{print}{summary.gjrm}(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)
 
}


\arguments{ 

\item{object}{A fitted \code{gjrm} object.}

\item{x}{\code{summary.gjrm} object produced by \code{summary.gjrm()}.}


\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals for the association parameter, dispersion coefficient etc. It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 
    
\item{digits}{Number of digits printed in output.} 

\item{signif.stars}{By default significance stars are printed alongside output.}     
    
\item{...}{Other arguments.}

}


\details{ 

\code{print.summary.gjrm} prints model term summaries. 


}



\value{

\item{tableP1}{Table containing parametric estimates, their standard errors, z-values and p-values for equation 1.}
\item{tableP2,tableP3, ...}{As above but for equation 2 and equations 3 and 4 if present.}
\item{tableNP1}{Table of nonparametric summaries for each smooth component including effective degrees of freedom, estimated rank, 
                approximate Wald statistic for testing the null hypothesis that the smooth term is zero and 
                corresponding p-value, for equation 1.}
\item{tableNP2,tableNP3, ...}{As above but for equation 2 and equations 3 and 4 if present.}
\item{n}{Sample size.}
\item{theta}{Estimated dependence parameter linking the two equations.}
\item{sigma1,sigma2}{Estimated distribution specific parameters for equations 1 and 2.}
\item{nu1,nu2}{Estimated distribution specific parameters for equations 1 and 2.}
\item{formula1,formula2,formula3, ...}{Formulas used for the model equations.}
\item{l.sp1,l.sp2,l.sp3, ...}{Number of smooth components in model equations.}
\item{t.edf}{Total degrees of freedom of the estimated bivariate model.}
\item{CItheta}{Interval(s) for \eqn{\theta}{\theta}.}
\item{CIsig1,CIsig2,CInu1,CInu2}{Intervals for distribution specific parameters} 


}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}




\section{WARNINGS}{

Note that the Kendall's tau (and related interval), as implemented here, is a valid measure of dependence for continuous margins and it will only provide a crude 
indication of dependence in other cases. 

} 



\keyword{information criteria}
\keyword{smooth}
\keyword{regression}







