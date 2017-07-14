\name{GJRM-package}
\alias{GJRM-package}
\docType{package}

\title{Generalised Joint Regression Modelling}

\description{

This package provides functions for fitting various generalised joint regression models with several types of covariate 
effects and distributions. Many modelling options are supported and all parameters of the joint distribution can be specified as 
flexible functions of covariates.

The orginal name of this package was \code{SemiParBIVProbit} which was designed 
to fit flexible bivariate binary response models. However, since then the package has expanded so much that its orginal name 
no longer gave a clue about all modelling options available. The new name should more closely reflect past, current and future developments.

The main fitting functions are listed below.

\code{SemiParBIV()} fits bivariate regression models with binary responses. This is useful to fit bivariate binary models in the presence of 
(i) non-random sample selection or (ii) associated responses/endogeneity or (iii) partial observability. This function includes
the \code{Model} argument which allows the user to fit the model in one of three situations mentioned above.

\code{copulaReg()} fits bivariate models with binary/discrete/continuous/survival margins in the presence of 
associated responses/endogeneity. 

\code{copulaSampleSel()} fits bivariate sample selection models with continuous/discrete response (instead of 
binary as it would be the case when using \code{SemiParBIV()}).  

\code{SemiParTRIV()} fits trivariate binary models (with and without double sample selection).  

\code{gamlss()} fits flexible univariate regression models where the response can be 
binary (only the extreme value distribution is allowed for, continuous, discrete and survival). The 
purpose of this function was only to provide, in some cases, starting values 
for the above functions, but it has now been made available in the form of a proper function should the user wish to fit 
univariate models using the general estimation approach of this package.

We are currently working on several multivariate extensions.


}

\details{
\code{GJRM} provides functions for fitting general joint models in various situations. The estimation approach is
based on a very generic penalized maximum likelihood based framework, where any (parametric) distribution can in principle be employed,
and the smoothers (representing several types of covariate effects) are set up using penalised regression splines.
Several marginal and copula distributions are available and the 
numerical routine carries out function minimization using a trust region algorithm in combination with
an adaptation of an automatic multiple smoothing parameter estimation procedure for GAMs (see \code{mgcv} for more details on this last point). The smoothers 
supported by this package are those available in \code{mgcv}. 

Confidence intervals for smooth components and nonlinear functions of the model
parameters are derived using a Bayesian approach. P-values for testing 
individual smooth terms for equality to the zero function are also provided and based on the approach
implemented in \code{mgcv}. The usual plotting and summary functions are also available. Model/variable 
selection is also possible via the use of shrinakge smoothers and/or information criteria. 


}


\author{

Giampiero Marra (University College London, Department of Statistical Science) and Rosalba Radice (Birkbeck, University of London, Department of Economics, Mathematics and Statistics)

with contributions from Panagiota Filippou.

Thanks to Bear Braumoeller (Department of Political Science, The Ohio State University) for suggesting the implementation of bivariate models with partial observability.

Thanks also to Carmen Cadarso and Francisco Gude for suggesting various modelling extensions and improvements, and importantly for sharing their vision for joint modelling. 

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

Part funded by EPSRC: EP/J006742/1

}

\references{

Key references:

Filippou P., Marra G. and Radice R. (in press), Penalized Likelihood Estimation of a Trivariate Additive Probit Model. \emph{Biostatistics}.

Marra G. and Radice R. (2011), Estimation of a Semiparametric Recursive Bivariate Probit in the Presence of Endogeneity. \emph{Canadian Journal of Statistics}, 39(2), 259-279.

Marra G. and Radice R. (2013), A Penalized Likelihood Estimation Approach to Semiparametric Sample Selection Binary Response Modeling. \emph{Electronic Journal of Statistics}, 7, 1432-1455.

Marra G. and Radice R. (2017), Bivariate Copula Additive Models for Location, Scale and Shape. \emph{Computational Statistics and Data Analysis}, 112, 99-113.

Marra G., Radice R., Barnighausen T., Wood S.N. and McGovern M.E. (2017), A Simultaneous Equation Approach to Estimating HIV Prevalence with Non-Ignorable Missing Responses. \emph{Journal of the American Statistical Association}, 112(518), 484-496.

Marra G. and Wyszynski K. (2016), Semi-Parametric Copula Sample Selection Models for Count Responses. \emph{Computational Statistics and Data Analysis}, 104, 110-129.

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.

Radice R., Marra G. and Wojtys M. (2016), Copula Regression Spline Models for Binary Outcomes. \emph{Statistics and Computing}, 26(5), 981-995. 

Wojtys M. and Marra G. (submitted). Copula-Based Generalized Additive Models with Non-Random Sample Selection.

}


\seealso{\code{\link{SemiParBIV}}, \code{\link{copulaReg}}, \code{\link{copulaSampleSel}}, \code{\link{gamlss}}, \code{\link{SemiParTRIV}}}


\keyword{endogeneity} 
\keyword{package}
\keyword{non-random sample selection} 
\keyword{partial observability} 
\keyword{copula} 
\keyword{penalised regression spline} 
\keyword{smooth} 
\keyword{regression}


