rob.const <- function(x, B = 100){

cont <- c("N", "N2", "GU", "rGU", "LO", "LN", "WEI","iG", "GA", "DAGUM", "SM", "BE", "FISK")
disc <- c("NBI", "NBII", "PIG", "PO", "ZTP") 
margin <- x$margin[1]

if(margin == "ZTP") rZTP <- function(n, mu) qpois(runif(n, dpois(0, mu), 1), mu)

n <- x$n

eta1   <- x$eta1
sigma2 <- x$sigma2 
nu     <- x$nu


sw <- NA

for(i in 1:B){

if(margin == "N")     y2s <- rNO(   n,    mu = eta1,         sigma = sqrt(sigma2)) 
if(margin == "N2")    y2s <- rNO(   n,    mu = eta1,         sigma = sigma2) 
if(margin == "GU")    y2s <- rGU(   n,    mu = eta1,         sigma = sqrt(sigma2)) 
if(margin == "rGU")   y2s <- rRG(   n,    mu = eta1,         sigma = sqrt(sigma2)) 
if(margin == "LO")    y2s <- rLO(   n,    mu = eta1,         sigma = sqrt(sigma2)) 
if(margin == "LN")    y2s <- rLOGNO(n,    mu = eta1,         sigma = sqrt(sigma2)) 
if(margin == "WEI")   y2s <- rWEI(  n,    mu = exp(eta1),    sigma = sqrt(sigma2)) 
if(margin == "iG")    y2s <- rIG(   n,    mu = exp(eta1),    sigma = sqrt(sigma2)) 
if(margin == "GA")    y2s <- rGA(   n,    mu = exp(eta1),    sigma = sqrt(sigma2)) 
if(margin == "GAi")   y2s <- rGA(   n,    mu = eta1,         sigma = sqrt(sigma2)) 
if(margin == "DAGUM") y2s <- rGB2(  n,    mu = exp(eta1),    sigma = sqrt(sigma2), nu = nu, tau = 1) 
if(margin == "SM")    y2s <- rGB2(  n,    mu = exp(eta1),    sigma = sqrt(sigma2), nu = 1,  tau = nu) 
if(margin == "BE")    y2s <- rBE(   n,    mu = plogis(eta1), sigma = sqrt(sigma2))
if(margin == "FISK")  y2s <- rGB2(  n,    mu = exp(eta1),    sigma = sqrt(sigma2), nu = 1,  tau = 1)
if(margin == "NBI")   y2s <- rNBI(  n,    mu = exp(eta1),    sigma = sqrt(sigma2)) 
if(margin == "NBII")  y2s <- rNBII( n,    mu = exp(eta1),    sigma = sqrt(sigma2)) 
if(margin == "PIG")   y2s <- rPIG(  n,    mu = exp(eta1),    sigma = sqrt(sigma2)) 
if(margin == "PO")    y2s <- rPO(   n,    mu = exp(eta1)) 
if(margin == "ZTP")   y2s <- rZTP(  n,    mu = exp(eta1))   



if(margin %in% cont) lpdf <- log(distrHsAT(y2s, eta1, sigma2, nu, margin)$pdf2) 
if(margin %in% disc) lpdf <- log(distrHsATDiscr2(y2s, eta1, sigma2, nu, margin)$pdf2) 


sw[i] <- sum(llpsi(lpdf, x$VC$rc)$d.psi)

                 
}

v1 <- 1/B*sum(sw/n)
v2 <- median(sw/n)
                
list(rc = x$VC$rc, sw = sw, m1 = v1, m2 = v2)


}
