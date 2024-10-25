print.cond.mv <- function(x, ...){

es <- format(x$res, digits = 3, trim = TRUE)

if(x$fun == "mean") cat("\nConditional mean with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")
if(x$fun == "variance") cat("\nConditional variance with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")


invisible(x)

}
