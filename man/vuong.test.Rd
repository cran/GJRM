\name{vuong.test}
\alias{vuong.test}

\title{Vuong test}

\description{ 
The Vuong test is likelihood-ratio-based tests that can be used for choosing between two non-nested models.}


\usage{

vuong.test(obj1, obj2, sig.lev = 0.05)

}

\arguments{ 

\item{obj1,obj2}{Objects of the two fitted bivariate non-nested models.} 

\item{sig.lev}{Significance level used for testing.} 

}


\details{

The Vuong test is a likelihood-ratio-based tests for model selection that use the 
Kullback-Leibler information criterion, and that can be employed for choosing between two bivariate models which 
are non-nested. 

The null hypothesis is that the two models are equally close to the actual model, whereas 
the alternative is that one model is closer. The test follows asymptotically a standard normal 
distribution under the null. Assume that the critical region is \eqn{(-c,c)}, where \eqn{c} is typically set to 1.96. If the value 
of the test is higher than \eqn{c} then we reject the null hypothesis 
that the models are equivalent in favor of model \code{obj1}. Viceversa if the value is smaller than \eqn{c}. If 
the value falls in \eqn{[-c,c]} then we cannot discriminate between the two competing models given the data. 

 
}


\value{

It returns a decision.

}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{

Vuong Q.H. (1989), Likelihood Ratio Tests for Model Selection and Non-Nested Hypotheses. \emph{Econometrica}, 57(2), 307-333. 

}



\examples{
## see examples for gjrm
}

\keyword{Vuong test} 
\keyword{likelihood ratio test} 









