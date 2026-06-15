\name{gjrm.pcc}
\alias{gjrm.pcc}

\title{Three-Dimensional PCC Models with Continuous Margins}

\description{
Fits a 3-dimensional copula model using the pair-copula construction approach by selecting and estimating three bivariate copulas, optionally 
incorporating covariates or smooth effects.
}

\usage{
gjrm.pcc(obj, data, cond.var, eq.theta  = NULL, simplified = FALSE, 
         criterion = "aic", sig.lev = 0.05, dof.grid = 3:12, 
         infl.fac = c(1, 1, 1), copulas = NULL, dofs = NULL)
}



\arguments{

\item{obj}{A list containing three \code{gamlss} objects. Each element has to be fitted using \code{GJRM::gamlss()}.}

  \item{data}{Compulsory \code{data.frame}.}

  \item{cond.var}{Integer (1, 2, or 3) indicating the variable to condition on in the third pair-copula.}
  
  \item{eq.theta}{Optional formula for covariate effects in the copula parameter.}
  
  \item{simplified}{Logical. If \code{TRUE}, then a simplified PCC is fitted.}
  
  \item{criterion}{Model selection criterion. One of: \code{"aic"}, \code{"bic"}.}

  \item{sig.lev}{Significance level for copula selection when using the Vuong and Clarke tests.}

  \item{dof.grid}{Integer vector of candidate degrees of freedom for the Student-t copula.}
  
\item{infl.fac}{Inflation factor for the effective degrees of freedom in the approximate AIC. 
Values greater than 1 increase the penalty on model complexity, leading to smoother and more parsimonious fits. 
Can be a scalar or a vector matching the number of pair-copula components.}
             
  \item{copulas}{Optional character vector specifying the copula family used for each pair-copula in the PCC. 
  Each element must correspond to a pair-copula in the decomposition (in the order defined by the model structure). 
  Allowed values depend on the implementation (e.g. Gaussian, Clayton, Gumbel, Student-t, etc.). 
  If \code{NULL}, the copula families are selected automatically using the specified \code{criterion}. 
  If a Student-t copula is specified (e.g., \code{"T"}), the corresponding degrees of freedom must be provided in \code{dofs}.}

  \item{dofs}{Optional numeric vector of degrees of freedom associated with Student-t copulas in \code{copulas}. 
  This argument must have the same length as \code{copulas}. 
  Entries are only used when the corresponding element of \code{copulas} is a Student-t copula (e.g., \code{"T"}); 
  otherwise they are ignored. If a Student-t copula is included in \code{copulas}, the corresponding element in \code{dofs} must be provided; 
  otherwise the function will stop with an error.}             
             
             
}

\value{
A list containing the fitted pair-copula models, selected copulas, estimated parameters, total AIC/BIC and additional model details. The usual summary and 
plotting functions can be employed for each single fitted model.
}

\keyword{generalised joint regression modelling}
\keyword{smooth} 
\keyword{regression spline}
\keyword{copula}
\keyword{marginal distribution}
\keyword{flexible copula regression modelling} 
\keyword{pair-copula construction}
\keyword{trivariate model}