\name{gt.bpm}
\alias{gt.bpm}

\title{Gradient test}

\description{ 
\code{gt.bpm} can be used to test the hypothesis of absence of endogeneity, correlated model equations/errors or non-random sample selection
in binary bivariate probit models.}


\usage{

gt.bpm(x)

}

\arguments{ 

\item{x}{A fitted \code{SemiParBIV} object as produced by \code{SemiParBIV()}.} 

}


\details{

The gradient test was first proposed by Terrell (2002) and it is based on classic likelihood 
theory. See Marra et al. (in press) for full details. 



}


\value{

It returns a numeric p-value corresponding to the null hypothesis that the correlation, \eqn{\theta}{\theta}, is equal to 0. 


}



\section{WARNINGS}{

This test's implementation is only valid for bivariate binary probit models with normal errors.

} 




\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}

\references{


Marra G., Radice R. and Filippou P. (2017), Regression Spline Bivariate Probit Models: A Practical Approach to Testing for Exogeneity. \emph{Communications in Statistics - Simulation and Computation}, 46(3), 2283-2298.

Terrell G. (2002), The Gradient Statistic. \emph{Computing Science and Statistics}, 34, 206-215.


}




\seealso{
\code{\link{SemiParBIV}}
}

\examples{
## see examples for SemiParBIV
}

\keyword{generalised joint regression modelling} 
\keyword{gradient test} 
\keyword{endogeneity}
\keyword{non-random sample selection}









