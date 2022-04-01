print.AT <- function(x, ...){





if(x$triv == FALSE){


if(x$type != "naive" && !(x$mar2 %in% x$bl) && x$eq == 1){

es <- format(x$Effects, digits = 3, trim=TRUE)

cat("\n")
print(es)
cat("\n")


}else{


if(x$mar2 %in% x$bl || ( !(x$mar2 %in% x$bl) && x$eq == 1) ) es <- format(x$res*100, digits = 3, trim=TRUE)

if(!(x$mar2 %in% x$bl) && x$eq == 2) es <- format(x$res, digits = 3, trim=TRUE)

if(x$mar2 %in% x$bl || ( !(x$mar2 %in% x$bl) && x$eq == 1) ) cat("\nTreatment effect (%) with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

if( !(x$mar2 %in% x$bl) && x$eq == 2) cat("\nTreatment effect with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")

if(x$Model != "ROY") cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")

}





}


if(x$triv == TRUE){


es <- format(x$res*100, digits = 3, trim=TRUE)

cat("\nTreatment effect (%) with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")
cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")



}



if(x$Model == "ROY"){


es <- format(x$res, digits = 3, trim = TRUE)

cat("\nAverage treatment effect with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")
cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")



}





invisible(x)

}
