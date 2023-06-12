mice.impute.copulaSS <- function(y, ry, x, mice.formula, margins, copula, ...) {

x <- as.data.frame(x) 

x[, "ry"] <- ry
x[, "y"]  <- y

out <- gjrm(mice.formula, data = x, model = "BSS", margins = margins, copula = copula)
return(imputeSS(out, 1)[[1]])

}



  