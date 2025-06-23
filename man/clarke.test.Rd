\name{clarke.test}
\alias{clarke.test}

\title{Clarke test}

\description{ 
The Clarke test is a likelihood-ratio-based test that can be used for choosing between two non-nested models.}


\usage{

clarke.test(obj1, obj2, sig.lev = 0.05)

}

\arguments{ 

\item{obj1,obj2}{Objects of the two fitted bivariate non-nested models.} 

\item{sig.lev}{Significance level used for testing.} 

}


\details{

The Clarke (2007) test is a likelihood-ratio-based tests for model selection that use the 
Kullback-Leibler information criterion, and that can be employed for choosing between two bivariate models which 
are non-nested. 

If the two models are statistically equivalent then the log-likelihood ratios of the 
observations should be evenly distributed around zero 
and around half of the ratios should be larger than zero. The test follows asymptotically a binomial distribution with 
parameters \eqn{n} and 0.5. Critical values can be obtained as shown in Clarke (2007). Intuitively, 
model \code{obj1} is preferred over \code{obj2} if the value of the test 
is significantly larger than its expected value under the null hypothesis (\eqn{n/2}), and vice versa. If 
the value is not significantly different from \eqn{n/2} then \code{obj1} can be thought of as equivalent to \code{obj2}.
 

}


\value{

It returns a decision.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Clarke K. (2007), A Simple Distribution-Free Test for Non-Nested Model Selection. \emph{Political Analysis}, 15, 347-363.



}



\examples{
## see examples for gjrm
}

\keyword{Clarke test} 
\keyword{likelihood ratio test} 









