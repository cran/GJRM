print.ATE <- function(x, ...){


if(x$type != "naive" && x$mar2 == "N" && x$eq == 1){

    es <- format(x$Effects, digits = 3, trim = TRUE)

    cat("\n")
    print(es)
    cat("\n")
    
}else{



es <- format(x$res, digits = 3, trim = TRUE)

cat("\nAverage treatment effect with ",(1-x$prob.lev)*100,"% interval:\n\n",sep="")
cat(es[2]," (",es[1],",",es[3],")\n\n",sep="")


     }





invisible(x)

}
