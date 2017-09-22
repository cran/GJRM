\name{LM.bpm}
\alias{LM.bpm}

\title{Lagrange Multiplier Test (Score Test)}

\description{ 
Before fitting a bivariate probit model, \code{LM.bpm} can be used to test the hypothesis of absence of endogeneity, 
correlated model equations/errors 
or non-random sample selection.}


\usage{

LM.bpm(formula, data = list(), weights = NULL, subset = NULL, Model, 
       hess = TRUE)

}

\arguments{ 


\item{formula}{A list of two formulas, one for equation 1 and the other for equation 2. \code{s} terms are used to specify 
               smooth smooth functions of predictors. Note that if \code{Model = "BSS"} then the first formula MUST refer 
               to the selection equation.} 

\item{data}{An optional data frame, list or environment containing the variables in the model.  If not found in \code{data}, the
            variables are taken from \code{environment(formula)}.} 

\item{weights}{Optional vector of prior weights to be used in fitting.} 

\item{subset}{Optional vector specifying a subset of observations to be used in the fitting process.}

\item{Model}{It indicates the type of model to be used in the analysis. Possible values are "B" (bivariate model) and 
             "BSS" (bivariate model with sample selection). The two marginal equations have probit links.}

\item{hess}{If \code{FALSE} then the expected (rather than observed) information matrix is employed.}


}


\details{

This Lagrange multiplier test (also known as score test) is used here for testing the null 
hypothesis that \eqn{\theta}{\theta} is equal to 0 (i.e. no endogeneity, non-random sample selection or 
correlated model equations/errors, depending 
on  the model being fitted). Its main advantage is that it does 
not require an estimate of the model parameter vector under the alternative hypothesis. Asymptotically, it takes a Chi-squared distribution 
with one degree of freedom. Full details can be found in Marra et al. (2014) and Marra et al. (2017).


}


\value{

It returns a numeric p-value corresponding to the null hypothesis that the correlation, \eqn{\theta}{\theta}, is equal to 0.

}


\section{WARNINGS}{

This test's implementation is ONLY valid for bivariate binary probit models with normal errors.

} 




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Marra G., Radice R. and Filippou P. (2017), Regression Spline Bivariate Probit Models: A Practical Approach to Testing for Exogeneity. \emph{Communications in Statistics - Simulation and Computation}, 46(3), 2283-2298.

Marra G., Radice R. and Missiroli S. (2014), Testing the Hypothesis of Absence of Unobserved Confounding in Semiparametric Bivariate Probit Models. \emph{Computational Statistics}, 29(3-4), 715-741.

}




\seealso{
\code{\link{gjrm}}
}

\examples{
## see examples for gjrm
}

\keyword{joint regression modelling} 
\keyword{lagrange multiplier test} 
\keyword{score test} 
\keyword{endogeneity}
\keyword{non-random sample selection}
\keyword{correlated equations/errors}











