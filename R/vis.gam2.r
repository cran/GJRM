vis.gam2 <- function (x, view = NULL, cond = list(), n.grid = 30, too.far = 0, 
    col = NA, color = "heat", contour.col = NULL, se = -1, 
    plot.type = "persp", zlim = NULL, nCol = 50, eq = eq, fun = fun, mar = mar, xx1 = NULL, 
    xxx1 = NULL, ...){
    
type <- "response"
    fac.seq <- function(fac, n.grid) {
        fn <- length(levels(fac))
        gn <- n.grid
        if (fn > gn) 
            mf <- factor(levels(fac))[1:gn]
        else {
            ln <- floor(gn/fn)
            mf <- rep(levels(fac)[fn], gn)
            mf[1:(ln * fn)] <- rep(levels(fac), rep(ln, fn))
            mf <- factor(mf, levels = levels(fac))
        }
        mf
    }
    dnm <- names(list(...))
    v.names <- names(x$var.summary)
    if (is.null(view)) {
        k <- 0
        view <- rep("", 2)
        for (i in 1:length(v.names)) {
            ok <- TRUE
            if (is.matrix(x$var.summary[[i]])) 
                ok <- FALSE
            else if (is.factor(x$var.summary[[i]])) {
                if (length(levels(x$var.summary[[i]])) <= 1) 
                  ok <- FALSE
            }
            else {
                if (length(unique(x$var.summary[[i]])) == 1) 
                  ok <- FALSE
            }
            if (ok) {
                k <- k + 1
                view[k] <- v.names[i]
            }
            if (k == 2) 
                break
        }
        if (k < 2) 
            stop("Model does not seem to have enough terms to do anything useful")
    }
    else {
        if (sum(view %in% v.names) != 2) 
            stop(gettextf("view variables must be one of %s", 
                paste(v.names, collapse = ", ")))
        for (i in 1:2) if (!inherits(x$var.summary[[view[i]]], 
            c("numeric", "factor"))) 
            stop("Don't know what to do with parametric terms that are not simple numeric or factor variables")
    }
    ok <- TRUE
    for (i in 1:2) if (is.factor(x$var.summary[[view[i]]])) {
        if (length(levels(x$var.summary[[view[i]]])) <= 1) 
            ok <- FALSE
    }
    else {
        if (length(unique(x$var.summary[[view[i]]])) <= 1) 
            ok <- FALSE
    }
    if (!ok) 
        stop(gettextf("View variables must contain more than one value. view = c(%s,%s).", 
            view[1], view[2]))
    if (is.factor(x$var.summary[[view[1]]])) 
        m1 <- fac.seq(x$var.summary[[view[1]]], n.grid)
    else {
        r1 <- range(x$var.summary[[view[1]]])
        m1 <- seq(r1[1], r1[2], length = n.grid)
    }
    if (is.factor(x$var.summary[[view[2]]])) 
        m2 <- fac.seq(x$var.summary[[view[2]]], n.grid)
    else {
        r2 <- range(x$var.summary[[view[2]]])
        m2 <- seq(r2[1], r2[2], length = n.grid)
    }
    v1 <- rep(m1, n.grid)
    v2 <- rep(m2, rep(n.grid, n.grid))
    newd <- data.frame(matrix(0, n.grid * n.grid, 0))
    for (i in 1:length(x$var.summary)) {
        ma <- cond[[v.names[i]]]
        if (is.null(ma)) {
            ma <- x$var.summary[[i]]
            if (is.numeric(ma)) 
                ma <- ma[2]
        }
        if (is.matrix(x$var.summary[[i]])) 
            newd[[i]] <- matrix(ma, n.grid * n.grid, ncol(x$var.summary[[i]]), 
                byrow = TRUE)
        else newd[[i]] <- rep(ma, n.grid * n.grid)
    }
    names(newd) <- v.names
    newd[[view[1]]] <- v1
    newd[[view[2]]] <- v2
    if (type == "link") 
        zlab <- paste("linear predictor")
    else if (type == "response") 
        zlab <- type
    else stop("type must be \"link\" or \"response\"")



###################
#####################
# standar errors not correct at the moment



           fv  <- predict.gam(x,    newdata = newd, se.fit = TRUE, type = "link")
           
           if(!is.null(xx1))  fv2 <- predict.gam(xx1,  newdata = newd, se.fit = TRUE, type = "link")
           if(!is.null(xxx1)) fv3 <- predict.gam(xxx1, newdata = newd, se.fit = TRUE, type = "link")
           

    #  if(mar %in% c("LO")){
    #          if(fun == "mean")     fv$fit <- fv$fit
    #          if(fun == "variance") fv$fit <- pi^2*exp(fv2$fit)/3
    #  } 
    #
    #
    #  if(mar %in% c("LN")){
    #       if(fun == "mean")     fv$fit <- exp(fv$fit)*sqrt(exp(exp(fv2$fit)))
    #       if(fun == "variance") fv$fit <- exp(exp(fv2$fit))*( exp(exp(fv2$fit)) - 1 )*exp(2*fv$fit)
    #                      }      




 if(mar %in% c("SM")){
  
     if(fun == "mean")     fv$fit <- exp(fv$fit)/gamma(exp(fv3$fit))*gamma( 1+1/sqrt(exp(fv2$fit)) )*gamma( -1/sqrt(exp(fv2$fit))+exp(fv3$fit) )             
     if(fun == "variance") fv$fit <- exp(fv$fit)^2*( gamma(1+2/sqrt(exp(fv2$fit)))*gamma(exp(fv3$fit))*gamma(-2/sqrt(exp(fv2$fit))+exp(fv3$fit))-gamma(1+1/sqrt(exp(fv2$fit)))^2*gamma(-1/sqrt(exp(fv2$fit))+exp(fv3$fit))^2 )
                          
                                                  
 } 



 if(mar %in% c("BE")){
 
    if(fun == "mean")     fv$fit <- exp(fv$fit)                 
    if(fun == "variance") fv$fit <- exp(fv$fit)*(1-exp(fv$fit))*exp(fv2$fit)
                                     
 } 
 
 
 if(mar %in% c("FISK")){
 
    if(fun == "mean")     fv$fit <- exp(fv$fit)*pi/sqrt(exp(fv2$fit))/sin(pi/sqrt(exp(fv2$fit))) 
    if(fun == "variance") fv$fit <- exp(fv$fit)^2*( 2*pi/sqrt(exp(fv2$fit))/sin(2*pi/sqrt(exp(fv2$fit)))-(pi/sqrt(exp(fv2$fit)))^2/sin(pi/sqrt(exp(fv2$fit)))^2 )
                                     
 }  
 




 if(mar %in% c("GU")){
 
    if(fun == "mean")     fv$fit <- fv$fit - 0.57722*sqrt(exp(fv2$fit)) 
    if(fun == "variance") fv$fit <- pi^2*exp(fv2$fit)/6
                                     
 } 
 
 
 
 if(mar %in% c("rGU")){
 
    if(fun == "mean")     fv$fit <- fv$fit + 0.57722*sqrt(exp(fv2$fit)) 
    if(fun == "variance") fv$fit <- pi^2*exp(fv2$fit)/6
                                     
 }  




 if(mar %in% c("LO")){
 
    if(fun == "mean")     fv$fit <- fv$fit 
    if(fun == "variance") fv$fit <- pi^2*exp(fv2$fit)/3
                                     
 } 
 
 
 if(mar %in% c("N")){
 
    if(fun == "mean")     fv$fit <- fv$fit 
    if(fun == "variance") fv$fit <- exp(fv2$fit)
                                     
 } 
 
 if(mar %in% c("N2")){
 
    if(fun == "mean")     fv$fit <- fv$fit 
    if(fun == "variance") fv$fit <- sqrt(exp(fv2$fit))
                                     
 }  
 
 
 

 if(mar %in% c("LN")){
 
    if(fun == "mean")     fv$fit <- exp(fv$fit)*sqrt(exp(exp(fv2$fit))) 
    if(fun == "variance") fv$fit <- exp(exp(fv2$fit))*( exp(exp(fv2$fit)) - 1 )*exp(2*fv$fit)
                                    
 } 
 

 if(mar %in% c("iG")){
 
    if(fun == "mean")     fv$fit <- exp(fv$fit) 
    if(fun == "variance") fv$fit <- exp(fv$fit)^3*exp(fv2$fit)
                                    
 }  
 
 
 
  if(mar %in% c("GA")){
  
     if(fun == "mean")     fv$fit <- exp(fv$fit) 
     if(fun == "variance") fv$fit <- exp(fv$fit)^2*exp(fv2$fit)
                                     
 } 
 
 
  if(mar %in% c("WEI")){
  
     if(fun == "mean")     fv$fit <- exp(fv$fit)*gamma(1+1/sqrt(exp(fv2$fit))) 
     if(fun == "variance") fv$fit <- exp(fv$fit)^2*( gamma(1+2/sqrt(exp(fv2$fit))) - gamma( 1+1/sqrt(exp(fv2$fit)) )^2  )
                                     
 }  
 

  if(mar %in% c("DAGUM")){
  
     if(fun == "mean")      fv$fit <- -(exp(fv$fit)/sqrt(exp(fv2$fit)))*gamma(-1/sqrt(exp(fv2$fit)))*gamma(1/sqrt(exp(fv2$fit))+exp(fv3$fit))/gamma(exp(fv3$fit))             
     if(fun == "variance")  fv$fit <- -(exp(fv$fit)/sqrt(exp(fv2$fit)))^2*( 2*sqrt(exp(fv2$fit))*gamma(-2/sqrt(exp(fv2$fit)))*gamma(2/sqrt(exp(fv2$fit)) + exp(fv3$fit))/gamma(exp(fv3$fit)) + ( gamma(-1/sqrt(exp(fv2$fit)))*gamma(1/sqrt(exp(fv2$fit)) + exp(fv3$fit))/gamma(exp(fv3$fit))  )^2   )
                                                                           
 } 
 
 
 
 if(mar %in% c("PO")){ 
 
    if(fun == "mean" || fun == "variance")  fv$fit <- exp(fv$fit) 
                                     
 } 
 
 
 
 if(mar %in% c("NBI", "PIG")){ 
 
    if(fun == "mean")     fv$fit <- exp(fv$fit) 
    if(fun == "variance") fv$fit <- exp(fv$fit) + sqrt(exp(fv2$fit))*exp(fv$fit)^2
                                     
 } 
 
 
 if(mar %in% c("NBII")){ 
 
    if(fun == "mean")     fv$fit <- exp(fv$fit) 
    if(fun == "variance") fv$fit <- ( 1 + sqrt(exp(fv2$fit)) )*exp(fv$fit)
                                     
 }  
 
 if(mar %in% c("ZTP")){ 
 
    if(fun == "mean")     fv$fit <- exp(fv$fit)/( 1 - exp(-exp(fv$fit)) ) 
    if(fun == "variance") fv$fit <- ( exp(fv$fit)*( 1 - exp(-exp(fv$fit))*(exp(fv$fit) + 1)) )/( 1 - exp(-exp(fv$fit)) )^2 
                                     
 }   
  
 





#####################
#####################

    z <- fv$fit
    if (too.far > 0) {
        ex.tf <- exclude.too.far(v1, v2, x$model[, view[1]], 
            x$model[, view[2]], dist = too.far)
        fv$se.fit[ex.tf] <- fv$fit[ex.tf] <- NA
    }
    if (is.factor(m1)) {
        m1 <- as.numeric(m1)
        m1 <- seq(min(m1) - 0.5, max(m1) + 0.5, length = n.grid)
    }
    if (is.factor(m2)) {
        m2 <- as.numeric(m2)
        m2 <- seq(min(m1) - 0.5, max(m2) + 0.5, length = n.grid)
    }
    if (se <= 0) {
        old.warn <- options(warn = -1)
        av <- matrix(c(0.5, 0.5, rep(0, n.grid - 1)), n.grid, 
            n.grid - 1)
        options(old.warn)
        max.z <- max(z, na.rm = TRUE)
        z[is.na(z)] <- max.z * 10000
        z <- matrix(z, n.grid, n.grid)
        surf.col <- t(av) %*% z %*% av
        surf.col[surf.col > max.z * 2] <- NA
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        }
        else {
            min.z <- min(fv$fit, na.rm = TRUE)
            max.z <- max(fv$fit, na.rm = TRUE)
        }
        surf.col <- surf.col - min.z
        surf.col <- surf.col/(max.z - min.z)
        surf.col <- round(surf.col * nCol)
        con.col <- 1
        if (color == "heat") {
            pal <- heat.colors(nCol)
            con.col <- 3
        }
        else if (color == "topo") {
            pal <- topo.colors(nCol)
            con.col <- 2
        }
        else if (color == "cm") {
            pal <- cm.colors(nCol)
            con.col <- 1
        }
        else if (color == "terrain") {
            pal <- terrain.colors(nCol)
            con.col <- 2
        }
        else if (color == "gray" || color == "bw") {
            pal <- gray(seq(0.1, 0.9, length = nCol))
            con.col <- 1
        }
        else stop("color scheme not recognised")
        if (is.null(contour.col)) 
            contour.col <- con.col
        surf.col[surf.col < 1] <- 1
        surf.col[surf.col > nCol] <- nCol
        if (is.na(col)) 
            col <- pal[as.array(surf.col)]
        z <- matrix(fv$fit, n.grid, n.grid)
        if (plot.type == "contour") {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                ifelse("main" %in% dnm, "", ",main=zlab"), ",...)", 
                sep = "")
            if (color != "bw") {
                txt <- paste("image(m1,m2,z,col=pal,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
                txt <- paste("contour(m1,m2,z,col=contour.col,zlim=c(min.z,max.z)", 
                  ifelse("add" %in% dnm, "", ",add=TRUE"), ",...)", 
                  sep = "")
                eval(parse(text = txt))
            }
            else {
                txt <- paste("contour(m1,m2,z,col=1,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
            }
        }
        else {
            stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
                ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), 
                ifelse("zlab" %in% dnm, "", ",zlab=zlab"), ",...)", 
                sep = "")
            if (color == "bw") {
                op <- par(bg = "white")
                txt <- paste("persp(m1,m2,z,col=\"white\",zlim=c(min.z,max.z) ", 
                  stub, sep = "")
                eval(parse(text = txt))
                par(op)
            }
            else {
                txt <- paste("persp(m1,m2,z,col=col,zlim=c(min.z,max.z)", 
                  stub, sep = "")
                eval(parse(text = txt))
            }
        }
    }
    else {
        if (color == "bw" || color == "gray") {
            subs <- paste("grey are +/-", se, "s.e.")
            lo.col <- "gray"
            hi.col <- "gray"
        }
        else {
            subs <- paste("red/green are +/-", se, "s.e.")
            lo.col <- "green"
            hi.col <- "red"
        }
        if (!is.null(zlim)) {
            if (length(zlim) != 2 || zlim[1] >= zlim[2]) 
                stop("Something wrong with zlim")
            min.z <- zlim[1]
            max.z <- zlim[2]
        }
        else {
            max.z <- max(fv$fit + fv$se.fit * se, na.rm = TRUE)
            min.z <- min(fv$fit - fv$se.fit * se, na.rm = TRUE)
            zlim <- c(min.z, max.z)
        }
        z <- fv$fit - fv$se.fit * se
        z <- matrix(z, n.grid, n.grid)
        if (plot.type == "contour") 
            warning("sorry no option for contouring with errors: try plot.gam")
        stub <- paste(ifelse("xlab" %in% dnm, "", ",xlab=view[1]"), 
            ifelse("ylab" %in% dnm, "", ",ylab=view[2]"), ifelse("zlab" %in% 
                dnm, "", ",zlab=zlab"), ifelse("sub" %in% dnm, 
                "", ",sub=subs"), ",...)", sep = "")
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
            dnm, "", ",border=lo.col"), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- fv$fit
        z <- matrix(z, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
            dnm, "", ",border=\"black\""), stub, sep = "")
        eval(parse(text = txt))
        par(new = TRUE)
        z <- fv$fit + se * fv$se.fit
        z <- matrix(z, n.grid, n.grid)
        txt <- paste("persp(m1,m2,z,col=col,zlim=zlim", ifelse("border" %in% 
            dnm, "", ",border=hi.col"), stub, sep = "")
        eval(parse(text = txt))
    }
}





























