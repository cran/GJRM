print.SATE <- function(x, ...){

es <- format(x$res, digits = 3, trim = TRUE) 

if(dim(x$res)[1] > 1) cat("\nSurvival average treatment effects with ",(1-x$prob.lev)*100,"% intervals:\n\n",sep="")
if(dim(x$res)[1] == 1) cat("\nSurvival average treatment effect with ",(1-x$prob.lev)*100,"% intervals:\n\n",sep="")

print(es)

invisible(x)

}
