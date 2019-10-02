pred.gp <- function(x, p = 0.5, newdata, n.sim = 100, prob.lev = 0.05){

mar <- x$margins[1]
bs  <- rMVN(n.sim, mean = x$coefficients, sigma = x$Vb)


if(!missing(newdata)){ 

if(mar %in% c("GPII","GPo") ) xi <- exp( predict(x, eq = 1, newdata = newdata) ) - 0.5
if(mar %in% c("GP","DGP") )   xi <-      predict(x, eq = 1, newdata = newdata)
if(mar %in% c("DGPII") )      xi <-      predict(x, eq = 1, newdata = newdata)^2




if(mar %in% c("GP", "DGP", "GPII","DGPII") ){
   if( !is.null(x$X2) ) sigma <- exp( predict(x, eq = 2, newdata = newdata) ) else sigma <- x$sigma
                                            }


if(mar %in% c("GPo") ){
   if( !is.null(x$X2) ) sigma <- exp( predict(x, eq = 2, newdata = newdata) )/(1 + xi) else sigma <- x$sigma/(1 + xi)
                      }


}



if(missing(newdata)){ 

if(mar %in% c("GPII","GPo") ) xi <- exp( x$eta1 ) - 0.5
if(mar %in% c("GP","DGP") )   xi <-      x$eta1
if(mar %in% c("DGPII") )      xi <-      x$eta1^2



if(mar %in% c("GP", "DGP", "GPII","DGPII") ){
          if( !is.null(x$X2) ) sigma <- exp( x$eta2 ) else sigma <- x$sigma
                                            }

if(mar %in% c("GPo") ){
          if( !is.null(x$X2) ) sigma <- exp( x$eta2 )/(1 + xi) else sigma <- x$sigma/(1 + xi)  
                      }



           
}


# the above can be made more efficient by creating a function that carries out
# the transformations


q.p <- (sigma/xi)*((1-p)^(-xi)-1)


###########
# INTERVALS
###########


ind1 <- (1:x$X1.d2)

if(!missing(newdata)) X1m  <- predict(x, eq = 1, newdata = newdata, type = "lpmatrix") 
if( missing(newdata)) X1m  <- predict(x, eq = 1,                    type = "lpmatrix") 



if( !is.null(x$X2) ){
   ind2 <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2)
   if(!missing(newdata)) X2m  <- predict(x, eq = 2, newdata = newdata, type = "lpmatrix")
   if( missing(newdata)) X2m  <- predict(x, eq = 2,                    type = "lpmatrix")   
                    } 


if( is.null(x$X2) ){
   ind2 <- x$X1.d2 + 1
   X2m  <- matrix(1, dim(X1m), 1) 
}


eta1S <- X1m%*%t(bs[, ind1])
eta2S <- X2m%*%t(bs[, ind2])



if(mar %in% c("GPII","GPo") ) xiS <- exp( eta1S ) - 0.5
if(mar %in% c("GP","DGP") )   xiS <-      eta1S   
if(mar %in% c("DGPII") )      xiS <-      eta1S^2


if(mar %in% c("GP", "DGP", "GPII","DGPII") ) sigmaS <- exp( eta2S )                                          
if(mar %in% c("GPo") )                       sigmaS <- exp( eta2S )/(1 + xiS)


q.pS <- (sigmaS/xiS)*((1-p)^(-xiS)-1)


CIqp  <- rowQuantiles(q.pS,   probs = c(prob.lev/2, 1-prob.lev/2), na.rm = TRUE)
CIxi  <- rowQuantiles(xiS,    probs = c(prob.lev/2, 1-prob.lev/2), na.rm = TRUE)
CIsig <- rowQuantiles(sigmaS, probs = c(prob.lev/2, 1-prob.lev/2), na.rm = TRUE)


list(qp = q.p, CIqp = CIqp, xi = xi , sigma = sigma, CIxi = CIxi, CIsig = CIsig)

}


