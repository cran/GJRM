llpsi <- function (mlk, tc){

psi    <- log1p(exp(mlk + tc)) - log1p(exp(tc))                                       # log( (1 + exp(mlk + tc))/(1 + exp(tc))  )  
d.psi  <- exp(mlk + tc)/(1 + exp(mlk + tc))                                           # exp(mlk + tc)/(1 + exp(tc))/((1 + exp(mlk + tc))/(1 + exp(tc)))
d2.psi <- (1 - exp(mlk + tc)/(1 + exp(mlk + tc))) * exp(mlk + tc)/(1 + exp(mlk + tc)) # exp(mlk + tc)/(1 + exp(mlk + tc)) - exp(mlk + tc) * exp(mlk + tc)/(1 + exp(mlk + tc))^2
 
list(psi = psi, d.psi = d.psi, d2.psi = d2.psi) 
 
}


