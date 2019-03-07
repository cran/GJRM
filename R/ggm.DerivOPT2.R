ggm.DerivOPT2 <- function(params, s, n, idx, lambda, pen, VC, w.alasso, gamma, a){

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

a.gr           <- sigma - s
a.gr           <- ( a.gr + t(a.gr) - diag(diag(a.gr)) )*sc.f
a.gr           <- a.gr[lower.tri(a.gr, diag = TRUE)]
a.gr.star      <- a.gr
a.gr.star[idx] <- a.gr[idx]*params[idx] 

G   <- -a.gr.star

S <- sc.f*Dpens(params, type = pen, lambda, w.alasso, gamma, a)
diag(S)[idx] <- 0
       
S2 <- S%*%params   




G <- G + S2

G


}


