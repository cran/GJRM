rob.const <- function(x, B = 100, left.trunc = 0){

cont <- c("tN","N", "N2", "GU", "rGU", "LO", "LN", "WEI","IG", "GA", "DAGUM", "SM", "TW","BE", "FISK","GP","GPII","GPo")
disc <- c("tNBI", "tNBII", "tPIG", "NBI", "NBII", "PIG", "P", "tP","DGP","DGPII","DGP0") 
margin <- x$margin[1]

n <- x$n

min.dn <- 1e-160 # x$VC$min.dn


eta1   <- x$eta1
sigma2 <- x$sigma2 
nu     <- x$nu


sw <- NA

for(i in 1:B){

y2s <- sim.resp(margin, n, eta1, sigma2, nu, setseed = FALSE, left.trunc = left.trunc)

if(margin %in% cont) lpdf <- log(distrHsAT(y2s, eta1, sigma2, nu, margin, min.dn = min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr, left.trunc = left.trunc)$pdf2) 
if(margin %in% disc) lpdf <- log(distrHsATDiscr2(y2s, eta1, sigma2, nu, margin, min.dn = min.dn, left.trunc = left.trunc)$pdf2) 


sw[i] <- sum(llpsi(lpdf, x$VC$rc)$d.psi)
                 
              }

v1 <- 1/B*sum(sw/n)
v2 <- median(sw/n)
                
list(rc = x$VC$rc, sw = sw, m1 = v1, m2 = v2)


}
