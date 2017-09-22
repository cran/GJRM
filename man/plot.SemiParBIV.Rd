\name{plot.SemiParBIV}
\alias{plot.SemiParBIV}



\title{Plotting function}


\description{It takes a fitted \code{gjrm} object produced 
by \code{gjrm()} and 
plots the estimated smooth functions on the scale of the linear predictors. This function is a 
wrapper of \code{plot.gam()} in \code{mgcv}. Please see 
the documentation of \code{plot.gam()} for full details. 

}

\usage{
\method{plot}{SemiParBIV}(x, eq, ...)

}


\arguments{ 
  
\item{x}{A fitted \code{gjrm} object.}

\item{eq}{The equation from which smooth terms should be considered for printing.}

\item{...}{Other graphics parameters to pass on to plotting commands, as described for \code{plot.gam()} in \code{mgcv}.}


}


\details{

This function produces plots showing the smooth terms of a fitted semiparametric bivariate probit model. In the case of 1-D smooths, the 
x axis of each plot is labelled using the name of the regressor, while the y axis is labelled as \code{s(regr, edf)} 
where \code{regr} is the regressor's name, and \code{edf} the effective degrees of freedom of the smooth. For 2-D smooths, perspective 
plots are produced with the x axes labelled with the first and second variable names and the y axis 
is labelled as \code{s(var1, var2, edf)}, which indicates the variables of which the term is a function and the \code{edf} for the term. 

If \code{seWithMean = TRUE} then the intervals include the uncertainty about the overall mean. Note that the smooths are still shown 
centred. The theoretical arguments 
and simulation study of Marra and Wood (2012) suggest that \code{seWithMean = TRUE} results in intervals with
close to nominal frequentist coverage probabilities. 

}



\value{
The function generates plots.
}



\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}



\references{

Marra G. and Wood S.N. (2012), Coverage Properties of Confidence Intervals for Generalized Additive Model Components. \emph{Scandinavian Journal of Statistics}, 39(1), 53-74.

}


\section{WARNING}{ 

The function can not deal with smooths of more than 2 variables. 


} 



\seealso{
\code{\link{gjrm}}
}



\keyword{smooth} 
\keyword{regression} 
\keyword{hplot}






