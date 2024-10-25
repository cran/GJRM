\name{polys.map}
\alias{polys.map}

\title{Geographic map with regions defined as polygons}

\usage{

polys.map(lm, z, scheme = "gray", lab = "", zlim, rev.col = FALSE, ...)



}

\arguments{

  \item{lm}{Named list of matrices where each matrix has two columns. The matrix rows each define the 
   vertex of a boundary polygon.}

  \item{z}{A vector of values associated with each area (item) of \code{lm}. }

  \item{scheme}{Possible values are \code{"heat"}, \code{"terrain"}, \code{"topo"}, \code{"cm"} and \code{"gray"}, indicating how to fill 
               the polygons in accordance with the value of \code{z}.}

  \item{lab}{label for plot.}
  
  \item{zlim}{If missing then the range of z will be chosen using \code{pretty(z)} otherwise the range provided will be used.}  
   
  \item{rev.col}{If \code{TRUE} then coloring scheme is reversed.}  
   
  \item{...}{other arguments to pass to plot.}

}



\value{It produces a plot.}



\description{ 

This function produces a map with geographic regions defined by polygons. It is essentially the same function as
\code{polys.plot()} in \code{mgcv} but with added arguments \code{zlim} and \code{rev.col} and a wider set of choices for 
\code{scheme}. 

}

\details{

See help file of \code{polys.plot} in \code{mgcv}.

}


\author{

Maintainer: Giampiero Marra \email{giampiero.marra@ucl.ac.uk}

}


\keyword{hplot} 
\keyword{smooth} 
\keyword{regression}
