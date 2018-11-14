approx.CLM <- function(num, den, epsilon){

w.den <- which(den == epsilon)
w.num <- num[w.den]
res   <- ifelse(10 * abs(w.num) > epsilon, w.num * epsilon, w.num)
num[w.den] <- res

return(num)
   
}



