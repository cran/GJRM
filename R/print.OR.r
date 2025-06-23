print.OR <- function(x, ...){



es <- format(x$res, digits = 3, trim=TRUE)

cat("\nOdds ratio with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")




invisible(x)

}

