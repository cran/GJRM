\name{lmc}
\alias{lmc}

\title{Linear Model Fitting with Constraints}

\description{ 

Linear model fitting with positivity and sum-to-one constraints on the model's coefficients. 

}



\usage{

lmc(y, X, start.v = NULL, lambda = 1, pen = "none", gamma = 1, a = 3.7) 
        
        
}



\arguments{ 

\item{y}{Response vector.} 

\item{X}{Design matrix.} 

\item{start.v}{Starting values.}

\item{lambda}{Tuning parameter.} 

\item{pen}{Type of penalty. Choices are: none, ridge, lasso, alasso, scad.} 

\item{gamma}{Power parameter of adaptive lasso.} 

\item{a}{Scad parameter.} 

}


\details{

Linear model fitting with positivity and sum-to-one constraints on the model's coefficients.


}


\value{
The function returns an object of class \code{lmc}.
}



\examples{

\dontrun{ 

library(GJRM)

set.seed(1)

n    <- 1000
beta <- c(0.07, 0.08, 0.21, 0.12, 0.15, 0.17, 0.2)
l    <- length(beta)
X    <- matrix(runif(n*l), n, l)

y    <- X\%*\%beta + rnorm(n)

out <- lmc(y, X)
conv.check(out)

out1 <- lmc(y, X, start.v = beta)
conv.check(out1)


coef(out)                    # estimated   coefficients
round(out$c.coefficients, 3) # constrained coefficients
sum(out$c.coefficients)

round(out1$c.coefficients, 3) 
sum(out1$c.coefficients)


# penalised estimation

out1 <- lmc(y, X, pen = "alasso", lambda = 0.02)
conv.check(out1)

coef(out1)                    
round(out1$c.coefficients, 3)
sum(out1$c.coefficients)


AIC(out, out1)
BIC(out, out1)

round(cbind(out$c.coefficients, out1$c.coefficients), 3)

# scad

n    <- 10000
beta <- c(0.2, 0, 0, 0.02, 0.01, 0.01, 0.01, 0.08, 0.21, 0.12, 0.15, 0.17, 0.02)
l    <- length(beta)
X    <- matrix(runif(n*l), n, l)

y    <- X\%*\%beta + rnorm(n)

out1 <- lmc(y, X, pen = "scad", lambda = 0.01)
conv.check(out1)

coef(out1)  
sum(out1$c.coefficients)
                  
round(cbind(beta, out1$c.coefficients), 2)

}

}

\keyword{linear model} 
\keyword{sum-to-one constraint} 
\keyword{positivity constraint} 


