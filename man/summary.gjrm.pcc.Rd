\name{summary.gjrm.pcc}
\alias{summary.gjrm.pcc}
\alias{print.summary.gjrm.pcc}

\title{Summary Method for \code{gjrm.pcc} Objects}

\description{
Produces summaries for pair-copula components of a fitted \code{gjrm.pcc} object.
The method extracts and displays information for each estimated pair-copula
in the PCC decomposition, including conditional copulas.
}

\usage{

\method{summary}{gjrm.pcc}(object, ...)

\method{print}{summary.gjrm.pcc}(x, ...)
}

\arguments{

\item{object}{A fitted \code{gjrm.pcc} object.}

\item{x}{An object of class \code{summary.gjrm.pcc} produced by \code{summary.gjrm.pcc()}.}

\item{...}{Additional arguments passed to other methods.}

}

\details{

The function \code{summary.gjrm.pcc} returns the object unchanged but assigns it
the class \code{"summary.gjrm.pcc"} so that the corresponding print method
\code{print.summary.gjrm.pcc} is dispatched.

The print method displays each estimated pair-copula component in a structured
format, using labelled blocks such as \code{C12}, \code{C23} and \code{C13|2}.

This method is intended to provide a clear and interpretable summary of the
dependence structure in the fitted PCC model.
}

\value{

An object of class \code{summary.gjrm.pcc}, which is identical to the input
\code{gjrm.pcc} object but with an additional class attribute to enable
custom printing.

}

\author{
Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}

\seealso{
\code{\link{gjrm.pcc}}, \code{\link{print.gjrm.pcc}}
}

\keyword{regression}
\keyword{copula}
\keyword{summary}
