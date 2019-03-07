rob.const <- function(x, B = 100){

cont <- c("N", "N2", "GU", "rGU", "LO", "LN", "WEI","iG", "GA", "DAGUM", "SM", "BE", "FISK","GP")
disc <- c("NBI", "NBII", "PIG", "PO", "ZTP","DGP") 
margin <- x$margin[1]

n <- x$n

eta1   <- x$eta1
sigma2 <- x$sigma2 
nu     <- x$nu


sw <- NA

for(i in 1:B){

y2s <- sim.resp(margin, n, eta1, sigma2, nu, setseed = FALSE)

if(margin %in% cont) lpdf <- log(distrHsAT(y2s, eta1, sigma2, nu, margin)$pdf2) 
if(margin %in% disc) lpdf <- log(distrHsATDiscr2(y2s, eta1, sigma2, nu, margin)$pdf2) 

sw[i] <- sum(llpsi(lpdf, x$VC$rc)$d.psi)
                 
              }

v1 <- 1/B*sum(sw/n)
v2 <- median(sw/n)
                
list(rc = x$VC$rc, sw = sw, m1 = v1, m2 = v2)


}
