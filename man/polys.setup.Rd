\name{polys.setup}
\alias{polys.setup}

\title{Set up geographic polygons}

\usage{

polys.setup(object)



}

\arguments{

  \item{object}{An RDS file object as extracted from http://www.gadm.org.}



}



\value{It produces a list with polygons (\code{polys}), and various names (\code{names0}, \code{names1} - first level of aggregation, 
       \code{names2} - second level of aggregation).}



\description{ 

This function creates geographic polygons in a format suitable for smoothing. 

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

Thanks to Guy Harling for suggesting the implementation of this function.


}


\examples{

?hiv

}


