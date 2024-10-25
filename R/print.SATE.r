print.SATE <- function(x, ...){

es <- format(x$res, digits = 3, trim = TRUE) 

if(x$full == FALSE ) cat("\nSurvival average treatment effect with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")
if(x$full == TRUE  ) cat("\nSurvival average treatment effects with ",(1-x$prob.lev)*100,"% intervals:\n\n",sep="")


if(x$full == FALSE ) cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

if(x$full == TRUE ) print(es)

invisible(x)

}
