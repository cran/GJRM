\name{mb}
\alias{mb}

\title{Nonparametric (worst-case and IV) Manski's bounds}

\description{ 
\code{mb} can be used to calculate the (worst-case and IV) Manski's bounds and confidence interval covering the true effect of interest
with a fixed probability.}



\usage{

mb(treat, outc, IV = NULL, Model, B = 100, sig.lev = 0.05)

}

\arguments{ 

\item{treat}{Binary treatment/selection variable.}
 
\item{outc}{Binary outcome variable.} 

\item{IV}{An instrumental binary variable can be used if available.} 

\item{Model}{Possible values are "B" (model with endogenous variable) and "BSS" (model with non-random sample selection).}

\item{B}{Number of bootstrap replicates. This is used to obtain some components needed for confidence interval calculations.} 

\item{sig.lev}{Significance level.}

}



\value{

\item{LB, UP}{Lower and upper bounds for the true effect of interest.}
\item{CI}{Confidence interval covering the true effect of interest with a fixed probability.}
\item{ate.ra}{Estimated effect of interest assuming random assignment.}



}



\details{

Based on Manski (1990), this function returns the nonparametric lower and upper (worst-case) Manski's bounds for the average 
treatment effect (ATE) when \code{Model = "B"} or prevalence when \code{Model = "BSS"}. When an IV is employed
the function returns IV Manski bounds.

For comparison, it also returns the estimated effect assuming random assignment (i.e., the treatment received or selection relies 
on the assumption of ignorable observed and unobserved selection). Note that this is equivalent to 
what provided by \code{\link{AT}} or \code{\link{prev}} when \code{type = "naive"}, and is different from what obtained
by \code{\link{AT}} or \code{\link{prev}} when \code{type = "univariate"} as observed confounders are accounted for
and the assumption here is of ignorable unobserved selection.

A confidence interval covering the true ATE/prevalence with a fixed probability is also provided. This is based on the approach 
described in Imbens and Manski (2004). NOTE that this interval is typically very close (if not identical) to the lower
and upper bounds.

The ATE can be at most 1 (or 100 in percentage) and the worst-case Manski's bounds have width 1. This means that 
0 is always included within the possibilites of these bounds. Nevertheless, this may be useful to check whether 
the effect from a bivariate recursive model is included within the possibilites of the bounds. 

When estimating a prevalance the worst-case Manski's bounds have width equal to 
the non-response probability,
which provides a measure of the uncertainty about the prevalence caused by non-response. Again, this may be useful to check whether 
the prevalence from a bivariate non-random sample selection model is included within the possibilites of the bounds.

See \code{\link{gjrm}} for some examples.

}




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\references{

Manski C.F. (1990), Nonparametric Bounds on Treatment Effects. \emph{American Economic Review, Papers and Proceedings}, 80(2), 319-323.

Imbens G.W. and Manski C.F (2004), Confidence Intervals for Partially Identified Parameters. \emph{Econometrica}, 72(6), 1845-1857.

}


\seealso{
\code{\link{gjrm}}  
}

\examples{
## see examples for gjrm
}

\keyword{Manski's bounds} 
\keyword{Nonparametric bounds} 
\keyword{Worst-case bounds} 
\keyword{ATE} 
\keyword{confidence interval} 
\keyword{average treatment effect} 
\keyword{prevalence} 










