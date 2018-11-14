ggm.Deriv <- function(params, s, n, idx, lambda = 1, pen = "ridge", VC){

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

a.gr           <- sigma - s
a.gr           <- ( a.gr + t(a.gr) - diag(diag(a.gr)) )*sc.f
a.gr           <- a.gr[lower.tri(a.gr, diag = TRUE)]
a.gr.star      <- a.gr
a.gr.star[idx] <- a.gr[idx]*params[idx] 



a.kr  <- -kronecker(sigma, sigma)
a.dup <- duplication.matrix(p) 
a.h   <- t(a.dup)%*%a.kr%*%a.dup
a.h   <- matrix(a.h, p1, p1)*sc.f



# new a.h allowing for omega*

a.h2            <- a.h
a.h2[idx,]      <- a.h[idx,]*params[idx]
a.h2[-idx, idx] <- t(a.h2[idx, -idx])
a.h2[idx,idx]   <- t(t(a.h2[idx,idx])*params[idx])

#

a.hdiag         <- diag(a.h)[idx]
a.gr.sigma      <- a.gr[idx]
dsig.dsigstar.2 <- params[idx]^2
d2sig.d2sigstar <- params[idx]

a.hdiag <- a.hdiag*dsig.dsigstar.2 + a.gr.sigma*d2sig.d2sigstar 
diag(a.h2)[idx] <- a.hdiag
a.h <- a.h2

a.h


#####################

res <- -ll
G   <- -a.gr.star
H   <- -a.h

# since we will be minimising

#####################
# penalty setup
#####################
                            
if(pen == "ridge"){

diag.el <- rep(1, p1)

diag.el[idx] <- 0  # since we do not want to penalise the variances/precisions

S  <- lambda*diag(diag.el) 
S1 <- 0.5*crossprod(params, S)%*%params # this is to fix as well with sc.f
S2 <- S%*%params   

}


if(pen == "lasso"){

#diag.el <- 1/sqrt(params^2 + 1e-08)  

pL <- penAprL(params)
pen0 <- pL$pen0
pen1 <- pL$pen1
pen2 <- pL$pen2

pen0[idx] <- pen1[idx] <- pen2[idx] <- 0

S1 <- sc.f*lambda*sum(pen0)
S2 <- sc.f*lambda*pen1
S  <- sc.f*lambda*diag(pen2) # 0.5 removed...

}
         

#########################################
# final quantities to use in optimisation
#########################################

Sres <- res
res  <- Sres + S1
G    <- G + S2
H    <- H + S


list(value = res, gradient = G, hessian = H, S = S, l = Sres, p = p, countPD = countPD)

}

