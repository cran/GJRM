\name{summary.SemiParBIV}
\alias{summary.SemiParBIV}
\alias{print.summary.SemiParBIV}

\title{SemiParBIV summary}

\description{
It takes a fitted \code{SemiParBIV} object and produces some summaries from it.

}

\usage{

\method{summary}{SemiParBIV}(object, n.sim = 100, prob.lev = 0.05, gm = FALSE, ...)
        
\method{print}{summary.SemiParBIV}(x, digits = max(3, getOption("digits") - 3), 
           signif.stars = getOption("show.signif.stars"), ...)    
                      
}


\arguments{ 

\item{object}{A fitted \code{SemiParBIV} object.}

\item{x}{\code{summary.SemiParBIV} object produced by \code{summary.SemiParBIV()}.}

\item{n.sim}{The number of simulated coefficient vectors from the posterior distribution of the estimated model parameters. This is used 
to calculate intervals for the association parameter, dispersion coefficient and other measures (e.g., gamma measure). It may be increased if
            more precision is required.} 

\item{prob.lev}{Probability of the left and right tails of the posterior distribution used for interval calculations.} 
     
\item{gm}{If TRUE then intervals for the gamma measure and odds ratio are calculated.}     
     
\item{digits}{Number of digits printed in output.} 

\item{signif.stars}{By default significance stars are printed alongside output.} 
         
         
\item{...}{Other arguments.}

}


\details{ 

Using some low level functions in \code{mgcv}, based on the results of Marra and Wood (2012), `Bayesian p-values' are returned for the 
smooth terms. These have better frequentist performance than their frequentist counterpart. See the help file of 
\code{summary.gam} in \code{mgcv} for further details. Covariate selection can also be achieved 
using a single penalty shrinkage approach as shown in Marra and Wood (2011). 

Posterior simulation is used to obtain intervals of nonlinear functions of parameters, such as the association and dispersion parameters
as well as the odds ratio and gamma measure discussed by Tajar et al. (2001) if \code{gm = TRUE}. 


\code{print.summary.SemiParBIV} prints model term summaries. 


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
\item{formula1,formula2,formula3, ...}{Formulas used for the model equations.}
\item{l.sp1,l.sp2,l.sp3, ...}{Number of smooth components in model equations.}
\item{t.edf}{Total degrees of freedom of the estimated bivariate model.}
\item{CItheta}{Interval(s) for \eqn{\theta}{\theta}.}
\item{n.sel}{Number of selected observations in the sample selection case.}
\item{OR, CIor}{Odds ratio and related CI. The odds ratio is a measure of association between binary random variables and is defined as 
                p00p11/p10p01. In the case of independence this ratio is equal to 1. It can take values in the range (-Inf, Inf) and 
                it does not depend on the marginal probabilities (Tajar et al., 2001). Interval is calculated using posterior simulation.}
\item{GM, CIgm}{Gamma measure and related CI. This measure of association was proposed by Goodman and Kruskal (1954). It is defined as 
                (\code{OR} - 1)/(\code{OR} + 1), can take values in the range (-1, 1) and does not depend on the marginal probabilities.
                Interval is calculated using posterior simulation.}
\item{tau, CItau}{Kendall's tau and respective intervals.}                

}

\references{

Marra G. and Wood S.N. (2011), Practical Variable Selection for Generalized Additive Models. \emph{Computational Statistics and Data Analysis}, 55(7), 2372-2387.

Marra G. and Wood S.N. (2012), Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, 39(1), 53-74.

Tajar M., Denuit M. and Lambert P. (2001), Copula-Type Representation for Random Couples with Bernoulli Margins. Discussion Papaer 0118, Universite Catholique De Louvain.

}



\section{WARNINGS}{

Note that the summary output will also indeed provide the Kendall's tau and related interval. This is a valid measure of dependence for continuous margins
but it is not for binary margins. However, it is still displayed for the sake of keeping the printed output consistent with that of
other models in the package.

} 



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\seealso{
\code{\link{AT}}, \code{\link{prev}}
}




\keyword{information criteria}
\keyword{smooth}
\keyword{regression}







