vis.gjrm <- function(x, eq, fun = NULL, ...){


#maybe only eq 1 and d 2 make sense
#we need to say distribution false true or something?
#check se fit
#will onlt work for two parameters distributons in both equations for now

 mar <- xx1 <- NULL

 if(missing(eq)) stop("You must provide the equation number.")
 if(!is.null(fun) && !(fun %in% c("mean", "variance")) ) stop("The value for fun can either be mean or variance.")
 if(!is.null(fun) && !(eq %in% c(1,2)) )                 stop("The value for eq can be either 1 or 2.")

 if(!is.null(fun)){ if(eq == 1) mar <- x$margins[1]; if(eq == 2) mar <- x$margins[2]  } 
 
 
 if(!is.null(fun) && !(x$margins[1] %in% c(x$VC$m2, x$VC$m2d)) && !(x$margins[2] %in% c(x$VC$m2, x$VC$m2d)) ) stop("This is currently implemented for margins with two parameters.\nGet in touch to check progress on the other cases.") 

 if(eq > x$l.flist) stop("The fitted model has a smaller number of equations.") 

 if(eq==1 && x$l.sp1==0) stop("There is no model component to plot.")   
 if(eq==2 && x$l.sp2==0) stop("There is no model component to plot.")   
 if(eq==3 && x$l.sp3==0) stop("There is no model component to plot.")    
 if(eq==4 && x$l.sp4==0) stop("There is no model component to plot.")   
 if(eq==5 && x$l.sp5==0) stop("There is no model component to plot.")   
 if(eq==6 && x$l.sp6==0) stop("There is no model component to plot.")   
 if(eq==7 && x$l.sp7==0) stop("There is no model component to plot.") 
 if(eq==8 && x$l.sp8==0) stop("There is no model component to plot.")   

 
 if(eq==1){ ss.plot <- x$gam1
            ind <- 1:x$X1.d2 
            } 
 if(eq==2){ ss.plot <- x$gam2
            ind <- (x$X1.d2 + 1):(x$X1.d2 + x$X2.d2) 
            }
 if(eq==3){ ss.plot <- x$gam3
            ind <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) }
            
 if(eq==4){ ss.plot <- x$gam4
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2) }
            
 if(eq==5){ ss.plot <- x$gam5
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2) }
            
 if(eq==6){ ss.plot <- x$gam6
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2) }   
            
 if(eq==7){ ss.plot <- x$gam7
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2) }             
               
 if(eq==8){ ss.plot <- x$gam8
            ind <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2 + x$X5.d2 + x$X6.d2 + x$X7.d2 + + x$X8.d2) }             
                            
           

  ss.plot$coefficients <- x$coefficients[ind]
  ss.plot$Vp <- x$Vb[ind,ind]
  ss.plot$Vp.t <- x$Vb.t[ind,ind]
  ss.plot$sig2 <- 1
  ss.plot$edf <- diag(x$F)[ind]
  ss.plot$scale.estimated <- FALSE 
  ss.plot$call$data <- x$call$data



# this works at the moment only for 2 parameter distributions 
# for both margins


if(!is.null(fun) && eq == 1 && mar %in% x$VC$m2){

 ss.plot1 <- x$gam3
 ind1     <- (x$X1.d2 + x$X2.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2) 

           ss.plot1$coefficients <- x$coefficients[ind1]
           ss.plot1$Vp <- x$Vb[ind1,ind1]
           ss.plot1$Vp.t <- x$Vb.t[ind1,ind1]
           ss.plot1$sig2 <- 1
           ss.plot1$edf <- diag(x$F)[ind1]
           ss.plot1$scale.estimated <- FALSE 
           ss.plot1$call$data <- x$call$data

}


if(!is.null(fun) && eq == 2 && mar %in% x$VC$m2){

 ss.plot1 <- x$gam4
 ind1     <- (x$X1.d2 + x$X2.d2 + x$X3.d2 + 1):(x$X1.d2 + x$X2.d2 + x$X3.d2 + x$X4.d2) 

           ss.plot1$coefficients <- x$coefficients[ind1]
           ss.plot1$Vp <- x$Vb[ind1,ind1]
           ss.plot1$Vp.t <- x$Vb.t[ind1,ind1]
           ss.plot1$sig2 <- 1
           ss.plot1$edf <- diag(x$F)[ind1]
           ss.plot1$scale.estimated <- FALSE 
           ss.plot1$call$data <- x$call$data

}   


if(is.null(fun))  suppressWarnings( vis.gam(ss.plot, ...) )
if(!is.null(fun)) suppressWarnings( vis.gam2(x = ss.plot, xx1 = ss.plot1, fun = fun, mar = mar, ...) ) 
  
  
}










































