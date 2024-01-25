\name{areds}
\alias{areds}
\docType{data}


\title{AREDS}

\description{Real dataset of bivariate interval and right censored data with 628 subjects 
and three covariates. The dataset is a reshaped version of the AREDS data from the \code{CopulaCenR} package. The dataset 
was selected from the Age-related Eye Disease Study (AREDS) (AREDS Group, 1999). The two events are the 
progression times (in years) to late-AMD in the left and right eyes.}

\usage{
data(areds)
}


\format{\code{war} is a 628 row data frame with the following columns 

\describe{
           
\item{t11, t12}{left and right bounds of the intervals for the left eye. If \code{t12 = NA} then the observation is right-censored.}
\item{t21, t22}{left and right bounds of the intervals for the right eye. If \code{t22 = NA} then the observation is right-censored.}
\item{SevScore1, SevScore2}{baseline AMD severity scores for left and right eyes, respectively. Possible values are: 4, 5, 6, 7, 8.}
\item{age}{age at baseline.}
\item{rs2284665}{a genetic variant covariate highly associated with late-AMD progression. Possible values are: 0, 1, 2.}
\item{cens1, cens2}{type of censoring for left and right eyes.}
\item{cens}{joint censoring indicator for left and right eyes.}


}


} 



\source{

Data are from:

AREDS Group (1999), The Age-Related Eye Disease Study (AREDS): design implications. AREDS report no. 1. \emph{Control Clinical Trials}, 20, 573-600.


}



