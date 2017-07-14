\name{VuongClarke}
\alias{VuongClarke}

\title{Vuong and Clarke tests}

\description{ 
The Vuong and Clarke tests are likelihood-ratio-based tests that can be used for choosing between two non-nested models.}


\usage{

VuongClarke(obj1, obj2, sig.lev = 0.05)

}

\arguments{ 

\item{obj1,obj2}{Objects of the two fitted bivariate non-nested models.} 

\item{sig.lev}{Significance level used for testing.} 

}


\details{

The Vuong (1989) and Clarke (2007) tests are likelihood-ratio-based tests for model selection that use the 
Kullback-Leibler information criterion. The implemented tests can be used for choosing between two bivariate models which 
are non-nested. 

In the Vuong test, the null hypothesis is that the two models are equally close to the actual model, whereas 
the alternative is that one model is closer. The test follows asymptotically a standard normal 
distribution under the null. Assume that the critical region is \eqn{(-c,c)}, where \eqn{c} is typically set to 1.96. If the value 
of the test is higher than \eqn{c} then we reject the null hypothesis 
that the models are equivalent in favor of model \code{obj1}. Viceversa if the value is smaller than \eqn{c}. If 
the value falls in \eqn{[-c,c]} then we cannot discriminate between the two competing models given the data. 

In the Clarke test, if the two models are statistically equivalent then the log-likelihood ratios of the 
observations should be evenly distributed around zero 
and around half of the ratios should be larger than zero. The test follows asymptotically a binomial distribution with 
parameters \eqn{n} and 0.5. Critical values can be obtained as shown in Clarke (2007). Intuitively, 
model \code{obj1} is preferred over \code{obj2} if the value of the test 
is significantly larger than its expected value under the null hypothesis (\eqn{n/2}), and vice versa. If 
the value is not significantly different from \eqn{n/2} then \code{obj1} can be thought of as equivalent to \code{obj2}.
 

}


\value{

It returns two decisions based on the tests and criteria discussed above.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Clarke K. (2007), A Simple Distribution-Free Test for Non-Nested Model Selection. \emph{Political Analysis}, 15, 347-363.

Vuong Q.H. (1989), Likelihood Ratio Tests for Model Selection and Non-Nested Hypotheses. \emph{Econometrica}, 57(2), 307-333. 



}



\examples{
## see examples for SemiParBIV
}

\keyword{Vuong test} 
\keyword{Clarke test} 
\keyword{likelihood ratio test} 









