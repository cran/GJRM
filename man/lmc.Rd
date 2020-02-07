\name{lmc}
\alias{lmc}

\title{Linear Model Fitting with Constraints}

\description{ 

Linear model fitting with positivity and sum-to-one constraints on the model's coefficients. 

}



\usage{

lmc(y, X)
        
        
        
        
}



\arguments{ 

\item{y}{Response vector.} 

\item{X}{Design matrix.} 

}


\details{

Linear model fitting with positivity and sum-to-one constraints on the model's coefficients.


}


\value{
The function returns an object of class \code{lmc}.
}



\examples{

\dontrun{ 

set.seed(1)

n    <- 1000
beta <- c(0.07, 0.08, 0.21, 0.12, 0.15, 0.17, 0.2)
l    <- length(beta)
X    <- matrix(runif(n*l), n, l)

y    <- X\%*\%beta + rnorm(n)

out <- lmc(y, X)
conv.check(out)


coef(out)                    # estimated   coefficients
round(out$c.coefficients, 3) # constrained coefficients


}

}

\keyword{linear model} 
\keyword{sum-to-one constraint} 
\keyword{positivity constraint} 


