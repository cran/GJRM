\name{hfunc}
\alias{hfunc}

\title{Compute h-functions (conditional distributions) for a bivariate copula}

\description{
\code{hfunc} computes both h-functions for a given bivariate copula, using a pair of uniform pseudo-observations and the copula parameters.
These h-functions are essential for vine copula constructions and likelihood evaluation in pair-copula decompositions.
}

\usage{
hfunc(p1, p2, theta, copula = "N", dof = NULL)
}

\arguments{
  \item{p1}{A numeric vector of pseudo-observations (marginal uniform variables) from the first variable.}
  \item{p2}{A numeric vector of pseudo-observations from the second variable.}
  \item{theta}{Dependence copula parameter.}
  \item{copula}{A character string indicating the bivariate copula family.}
  \item{dof}{Degrees of freedom parameter, only required for the Student t-copula.}
}

\details{
This function computes the so called \emph{h-functions} which are required when performing pair-copula constructions, simulating from vine copulas or evaluating vine log-likelihoods.

The function returns a named list with two components:
\itemize{
  \item \code{h1}: Partial derivative of \eqn{C(u_1, u_2)} with respect to \eqn{u_1}.
  \item \code{h2}: Partial derivative of \eqn{C(u_1, u_2)} with respect to \eqn{u_2}.
}
}


\author{
Giampiero Marra \email{giampiero.marra@ucl.ac.uk}
}


\keyword{copula}
\keyword{conditional distribution}
\keyword{vine}
