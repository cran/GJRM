mice.impute.copulaSS <- function(y, ry, x, mice.formula, margins, BivD, ...) {

x <- as.data.frame(x) 

x[, "ry"] <- ry
x[, "y"]  <- y

out <- gjrm(mice.formula, data = x, Model = "BSS", margins = margins, BivD = BivD)
return(imputeSS(out, 1)[[1]])

}



  