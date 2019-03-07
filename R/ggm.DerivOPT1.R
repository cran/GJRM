ggm.DerivOPT1 <- function(params, s, n, idx, lambda, pen, VC, w.alasso, gamma, a){


params[idx] <- esp.tr(params[idx], "N")$vrb # this does the exp of params[idx]
                                            # and avoids the 0 problem

p     <- ncol(s)
p1    <- length(params)
omega <- matrix(0, p, p)

omega[lower.tri(omega, diag = TRUE)] <- params
omega <- t(omega) + omega - diag(diag(omega))

#### check PD ####

res.omega <- PDef(omega)

omega <- res.omega$res
sigma <- res.omega$res.inv

countPD <- VC$my.env$countPD

if(res.omega$check.eigen == TRUE){ params <- omega[lower.tri(omega, diag = TRUE)]; countPD <- countPD + 1; VC$my.env$countPD <- countPD}

##################

sc.f <- n*0.5

ll <- ( determinant(omega)$modulus[1]-sum(omega*s) )*sc.f

res <- -ll


S <- sc.f*Dpens(params, type = pen, lambda, w.alasso, gamma, a)
diag(S)[idx] <- 0
       
S1 <- 0.5*crossprod(params, S)%*%params 


Sres <- res
res  <- Sres + S1

res


}

