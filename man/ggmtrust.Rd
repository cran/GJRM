\name{ggmtrust}
\alias{ggmtrust}

\title{ggmtrust for penalised network}

\description{ 
penalised network, work in progress.
 
}



\usage{

ggmtrust(s, n, data = NULL, lambda = 1, pen = "lasso", params = NULL, method = "BHHH", 
         w.alasso = NULL, gamma = 1, a = 3.7)
        
        
        
        
}



\arguments{ 

\item{s}{Sample covariance matrix.} 

\item{n}{Sample size.} 

\item{data}{Data.} 

\item{lambda}{Regularisation parameter.} 

\item{pen}{Either "lasso" or "ridge".}

\item{params}{If different from null then these are taken as the starting values.}

\item{method}{Either "H" or "BHHH".}

\item{w.alasso}{weight for alasso.}

\item{gamma}{alasso param.}

\item{a}{scad param.}


}


\details{

penalised network, work in progress.


}


\value{
The function returns an object of class \code{ggmtrust}.
}


