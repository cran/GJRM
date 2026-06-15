gjrm.pcc <- function(obj, data, cond.var, eq.theta = NULL, simplified = FALSE, criterion = "aic", 
                     sig.lev = 0.05, dof.grid = 3:12, infl.fac = c(1, 1, 1), copulas = NULL, dofs = NULL){ 

#if (!is.list(obj) || length(obj) != 3 || !all(sapply(obj, inherits, "gamlss"))) stop("The object must be a list of three 'gamlss' objects.")

if (!is.list(obj) || length(obj) != 3 || !all(sapply(obj, function(x) inherits(x, "gamlss") && inherits(x, "gjrm")))) {
  stop("Currently, the object must be a list of three 'gamlss' objects obtained using the GJRM package.")
}


obj1 <- obj[[1]]; obj2 <- obj[[2]]; obj3 <- obj[[3]]

b12 <- b13 <- b23.1 <- b23 <- b13.2 <- b12.3 <- NULL

if(missing(data)) stop("You must provide a data frame.")
if(!missing(data)){ 
                   if(!is.data.frame(data)) stop("You must provide a data frame.")
                  }

tmp_p1 <- ".gjrm_p1"
tmp_p2 <- ".gjrm_p2"
tmp_p3 <- ".gjrm_p3"


# safety check
reserved_names <- c(tmp_p1, tmp_p2, tmp_p3)
if (any(reserved_names %in% names(data))) {
  stop("The data contains reserved variable names: ",
       paste(reserved_names, collapse = ", "),
       ". Please rename them before calling this function.")
}


data[[tmp_p1]] <- obj1$p
data[[tmp_p2]] <- obj2$p
data[[tmp_p3]] <- obj3$p

pf1 <- obj1$fit$pdf2 # to use for count case
pf2 <- obj2$fit$pdf2
pf3 <- obj3$fit$pdf2


if(!(cond.var %in% c(1:3))) stop("The cond.var argument can be set to 1, 2 or 3.")
if (!(criterion %in% c("vuong", "clarke", "aic", "bic"))) stop("Invalid criterion argument. Must be one of: 'aic', 'bic', 'vuong', 'clarke'.")


if (!is.null(copulas) && !is.null(dofs)) {
  if (length(dofs) != length(copulas)) {
    stop("copulas and dofs must have the same length.")
  }
}

if (!is.null(copulas)) {

  # positions where t-copula is used
  t_pos <- which(copulas == "T")

  if (length(t_pos) > 0) {

    # dofs must be provided
    if (is.null(dofs)) {
      stop("dofs must be provided when copulas include 'T'.")
    }

    # check missing or NA dofs at T positions
    if (any(is.na(dofs[t_pos]))) {
      stop("Missing dofs for copulas equal to 'T' at positions: ",
           paste(t_pos[is.na(dofs[t_pos])], collapse = ", "))
    }
  }
}

if (length(infl.fac) == 1) {
  infl.fac <- rep(infl.fac, 3)
} else if (length(infl.fac) != 3) {
  stop("infl.fac must be either length 1 or length 3.")
}


############################################################################################################################
## EQUATIONS' SET UP
############################################################################################################################

if( is.null(eq.theta) ){

	
eq12 <- list(as.formula(paste(tmp_p1, "~ 1")),
             as.formula(paste(tmp_p2, "~ 1")))
eq13 <- list(as.formula(paste(tmp_p1, "~ 1")),
             as.formula(paste(tmp_p3, "~ 1")))
             
eq23 <- list(as.formula(paste(tmp_p2, "~ 1")),
             as.formula(paste(tmp_p3, "~ 1")))	
	
	
			}else{

	
eq12 <- list(as.formula(paste(tmp_p1, "~ 1")),
             as.formula(paste(tmp_p2, "~ 1")),
             eq.theta)
eq13 <- list(as.formula(paste(tmp_p1, "~ 1")),
             as.formula(paste(tmp_p3, "~ 1")),
             eq.theta)
eq23 <- list(as.formula(paste(tmp_p2, "~ 1")),
             as.formula(paste(tmp_p3, "~ 1")),
             eq.theta)	
	
	
	
			}



if(simplified == TRUE){ 

  if(is.null(eq.theta)) {
    eqhs <- list(h.fm ~ 1, h.sm ~ 1)
  } else {
    eqhs <- list(h.fm ~ 1, h.sm ~ 1, eq.theta)
  }
                   
} else {

  if(cond.var == 1){
    if(is.null(eq.theta)){ 
      eqhs <- list(h.fm ~ 1, h.sm ~ 1,
                   as.formula(paste("~ s(", tmp_p1, ")", sep = ""))) 
    } else {
      eq.theta <- update(eq.theta,
                         as.formula(paste("~ . + s(", tmp_p1, ")", sep = "")))
      eqhs     <- list(h.fm ~ 1, h.sm ~ 1, eq.theta)
    } 
  }
		
  if(cond.var == 2){
    if(is.null(eq.theta)){ 
      eqhs <- list(h.fm ~ 1, h.sm ~ 1,
                   as.formula(paste("~ s(", tmp_p2, ")", sep = ""))) 
    } else {
      eq.theta <- update(eq.theta,
                         as.formula(paste("~ . + s(", tmp_p2, ")", sep = "")))
      eqhs     <- list(h.fm ~ 1, h.sm ~ 1, eq.theta)
    } 
  }

  if(cond.var == 3){
    if(is.null(eq.theta)){ 
      eqhs <- list(h.fm ~ 1, h.sm ~ 1,
                   as.formula(paste("~ s(", tmp_p3, ")", sep = ""))) 
    } else {
      eq.theta <- update(eq.theta,
                         as.formula(paste("~ . + s(", tmp_p3, ")", sep = "")))
      eqhs     <- list(h.fm ~ 1, h.sm ~ 1, eq.theta)
    } 
  }		
		
}






############################################################################################################################
############################################################################################################################





if(cond.var == 1){  
                 

              if(is.null(copulas))  b12 <- fix.fenv(cop.mod.sel(eq12, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[1]))
 
 
              if(!is.null(copulas)) b12 <- gjrm(eq12, fixed.margins = TRUE, data = data, copula = copulas[1], dof = dofs[1], infl.fac = infl.fac[1])
 
 
 
              if(b12$copula %in% b12$BivD2){
               
                h.fm <- NA
                
                h.fm[b12$fit$teta.ind1] <- hfunc(data[[tmp_p1]][b12$fit$teta.ind1], data[[tmp_p2]][b12$fit$teta.ind1], b12$theta[b12$fit$teta.ind1], b12$fit$Cop1, b12$dof)$h1
                h.fm[b12$fit$teta.ind2] <- hfunc(data[[tmp_p1]][b12$fit$teta.ind2], data[[tmp_p2]][b12$fit$teta.ind2], b12$theta[b12$fit$teta.ind2], b12$fit$Cop2, b12$dof)$h1
              
                data$h.fm <- h.fm 

              }else data$h.fm <- hfunc(data[[tmp_p1]], data[[tmp_p2]], b12$theta, b12$copula, b12$dof)$h1 


              if(is.null(copulas))  b13 <- fix.fenv(cop.mod.sel(eq13, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[2]))   
              
              if(!is.null(copulas)) b13 <- gjrm(eq13, fixed.margins = TRUE, data = data, copula = copulas[2], dof = dofs[2], infl.fac = infl.fac[2])   
              
              
              if(b13$copula %in% b13$BivD2){
 
               h.sm <- NA
 
               h.sm[b13$fit$teta.ind1] <- hfunc(data[[tmp_p1]][b13$fit$teta.ind1], data[[tmp_p3]][b13$fit$teta.ind1], b13$theta[b13$fit$teta.ind1], b13$fit$Cop1, b13$dof)$h1
               h.sm[b13$fit$teta.ind2] <- hfunc(data[[tmp_p1]][b13$fit$teta.ind2], data[[tmp_p3]][b13$fit$teta.ind2], b13$theta[b13$fit$teta.ind2], b13$fit$Cop2, b13$dof)$h1
              
                data$h.sm <- h.sm 

              }else data$h.sm <- hfunc(data[[tmp_p1]], data[[tmp_p3]], b13$theta, b13$copula, b13$dof)$h1
 
 
 
 
              if(is.null(copulas))  b23.1 <- fix.fenv(cop.mod.sel(eqhs, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[3]))  
 
              if(!is.null(copulas)) b23.1 <- gjrm(eqhs, fixed.margins = TRUE, data = data, copula = copulas[3], dof = dofs[3], infl.fac = infl.fac[3])  
 
 
              # this is needed for Rosenblatt residuals
              
              if(b23.1$copula %in% b23.1$BivD2){
 
               third.res <- NA
 
               third.res[b23.1$fit$teta.ind1] <- hfunc(data$h.fm[b23.1$fit$teta.ind1], data$h.sm[b23.1$fit$teta.ind1], b23.1$theta[b23.1$fit$teta.ind1], b23.1$fit$Cop1, b23.1$dof)$h1
               third.res[b23.1$fit$teta.ind2] <- hfunc(data$h.fm[b23.1$fit$teta.ind2], data$h.sm[b23.1$fit$teta.ind2], b23.1$theta[b23.1$fit$teta.ind2], b23.1$fit$Cop2, b23.1$dof)$h1
               

              }else third.res <- hfunc(data$h.fm, data$h.sm, b23.1$theta, b23.1$copula, b23.1$dof)$h1
  
 

              aicT      <- AIC(b12) + AIC(b13) + AIC(b23.1)                    
              bicT      <- BIC(b12) + BIC(b13) + BIC(b23.1)

              cop.sel   <- c(b12$copula, b13$copula, b23.1$copula)   
              dof.vec   <- c(b12$dof, b13$dof, b23.1$dof)
              
              ros.resids <- cbind(obj1$p, data$h.fm, third.res)

                            
                 } 




if(cond.var == 2){  
                           
              if(is.null(copulas))  b12 <- fix.fenv(cop.mod.sel(eq12, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[1]))
              if(!is.null(copulas)) b12 <- gjrm(eq12, fixed.margins = TRUE, data = data, copula = copulas[1], dof = dofs[1], infl.fac = infl.fac[1])
            
 
               if(b12$copula %in% b12$BivD2){
                
                 h.fm <- NA
                 
                 h.fm[b12$fit$teta.ind1] <- hfunc(data[[tmp_p1]][b12$fit$teta.ind1], data[[tmp_p2]][b12$fit$teta.ind1], b12$theta[b12$fit$teta.ind1], b12$fit$Cop1, b12$dof)$h2
                 h.fm[b12$fit$teta.ind2] <- hfunc(data[[tmp_p1]][b12$fit$teta.ind2], data[[tmp_p2]][b12$fit$teta.ind2], b12$theta[b12$fit$teta.ind2], b12$fit$Cop2, b12$dof)$h2
               
                 data$h.fm <- h.fm 
 
              }else data$h.fm <- hfunc(data[[tmp_p1]], data[[tmp_p2]], b12$theta, b12$copula, b12$dof)$h2




              if(is.null(copulas))  b23 <- fix.fenv(cop.mod.sel(eq23, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[2])) 
              if(!is.null(copulas)) b23 <- gjrm(eq23, fixed.margins = TRUE, data = data, copula = copulas[2], dof = dofs[2], infl.fac = infl.fac[2]) 
              
               
                if(b23$copula %in% b23$BivD2){
                 
                  h.sm <- NA
                  
                  h.sm[b23$fit$teta.ind1] <- hfunc(data[[tmp_p2]][b23$fit$teta.ind1], data[[tmp_p3]][b23$fit$teta.ind1], b23$theta[b23$fit$teta.ind1], b23$fit$Cop1, b23$dof)$h1
                  h.sm[b23$fit$teta.ind2] <- hfunc(data[[tmp_p2]][b23$fit$teta.ind2], data[[tmp_p3]][b23$fit$teta.ind2], b23$theta[b23$fit$teta.ind2], b23$fit$Cop2, b23$dof)$h1
                
                  data$h.sm <- h.sm 
  
              }else data$h.sm <- hfunc(data[[tmp_p2]], data[[tmp_p3]], b23$theta, b23$copula, b23$dof)$h1
 
 
 
 
              if(is.null(copulas))  b13.2 <- fix.fenv(cop.mod.sel(eqhs, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[3]))
              if(!is.null(copulas)) b13.2 <- gjrm(eqhs, fixed.margins = TRUE, data = data, copula = copulas[3], dof = dofs[3], infl.fac = infl.fac[3])
  
  
              # this is needed for Rosenblatt residuals
              
              if(b13.2$copula %in% b13.2$BivD2){
 
               third.res <- NA
 
               third.res[b13.2$fit$teta.ind1] <- hfunc(data$h.fm[b13.2$fit$teta.ind1], data$h.sm[b13.2$fit$teta.ind1], b13.2$theta[b13.2$fit$teta.ind1], b13.2$fit$Cop1, b13.2$dof)$h1
               third.res[b13.2$fit$teta.ind2] <- hfunc(data$h.fm[b13.2$fit$teta.ind2], data$h.sm[b13.2$fit$teta.ind2], b13.2$theta[b13.2$fit$teta.ind2], b13.2$fit$Cop2, b13.2$dof)$h1
               

              }else third.res <- hfunc(data$h.fm, data$h.sm, b13.2$theta, b13.2$copula, b13.2$dof)$h1  
  
  
  
              aicT      <- AIC(b12) + AIC(b23) + AIC(b13.2)    
              bicT      <- BIC(b12) + BIC(b23) + BIC(b13.2)
              
              cop.sel   <- c(b12$copula, b23$copula, b13.2$copula)
              dof.vec   <- c(b12$dof, b23$dof, b13.2$dof)
              
              
              
              ros.resids <- cbind(obj2$p, data$h.fm, third.res)
 
              
  
                 }  
                 
                 
                 
if(cond.var == 3){  
                           
              if(is.null(copulas))  b13 <- fix.fenv(cop.mod.sel(eq13, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[1])) 
              if(!is.null(copulas)) b13 <- gjrm(eq13, fixed.margins = TRUE, data = data, copula = copulas[1], dof = dofs[1], infl.fac = infl.fac[1]) 
 
 
              if(b13$copula %in% b13$BivD2){
 
               h.fm <- NA
 
               h.fm[b13$fit$teta.ind1] <- hfunc(data[[tmp_p1]][b13$fit$teta.ind1], data[[tmp_p3]][b13$fit$teta.ind1], b13$theta[b13$fit$teta.ind1], b13$fit$Cop1, b13$dof)$h2
               h.fm[b13$fit$teta.ind2] <- hfunc(data[[tmp_p1]][b13$fit$teta.ind2], data[[tmp_p3]][b13$fit$teta.ind2], b13$theta[b13$fit$teta.ind2], b13$fit$Cop2, b13$dof)$h2
              
                data$h.fm <- h.fm 

              }else data$h.fm <- hfunc(data[[tmp_p1]], data[[tmp_p3]], b13$theta, b13$copula, b13$dof)$h2 
 
 
              if(is.null(copulas))  b23 <- fix.fenv(cop.mod.sel(eq23, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[2])) 
              if(!is.null(copulas)) b23 <- gjrm(eq23, fixed.margins = TRUE, data = data, copula = copulas[2], dof = dofs[2], infl.fac = infl.fac[2]) 
             
              
              
                if(b23$copula %in% b23$BivD2){
                 
                  h.sm <- NA
                  
                  h.sm[b23$fit$teta.ind1] <- hfunc(data[[tmp_p2]][b23$fit$teta.ind1], data[[tmp_p3]][b23$fit$teta.ind1], b23$theta[b23$fit$teta.ind1], b23$fit$Cop1, b23$dof)$h2
                  h.sm[b23$fit$teta.ind2] <- hfunc(data[[tmp_p2]][b23$fit$teta.ind2], data[[tmp_p3]][b23$fit$teta.ind2], b23$theta[b23$fit$teta.ind2], b23$fit$Cop2, b23$dof)$h2
                
                  data$h.sm <- h.sm 
  
              }else data$h.sm <- hfunc(data[[tmp_p2]], data[[tmp_p3]], b23$theta, b23$copula, b23$dof)$h2



 
              if(is.null(copulas))  b12.3 <- fix.fenv(cop.mod.sel(eqhs, data = data, sig.lev = sig.lev, crit = criterion, dof.grid = dof.grid, infl.fac = infl.fac[3])) 
              if(!is.null(copulas)) b12.3 <- gjrm(eqhs, fixed.margins = TRUE, data = data, copula = copulas[3], dof = dofs[3], infl.fac = infl.fac[3]) 
   
             # this is needed for Rosenblatt residuals
              
              if(b12.3$copula %in% b12.3$BivD2){
 
               third.res <- NA
 
               third.res[b12.3$fit$teta.ind1] <- hfunc(data$h.fm[b12.3$fit$teta.ind1], data$h.sm[b12.3$fit$teta.ind1], b12.3$theta[b12.3$fit$teta.ind1], b12.3$fit$Cop1, b12.3$dof)$h1
               third.res[b12.3$fit$teta.ind2] <- hfunc(data$h.fm[b12.3$fit$teta.ind2], data$h.sm[b12.3$fit$teta.ind2], b12.3$theta[b12.3$fit$teta.ind2], b12.3$fit$Cop2, b12.3$dof)$h1
               

              }else third.res <- hfunc(data$h.fm, data$h.sm, b12.3$theta, b12.3$copula, b12.3$dof)$h1  
              
   
   
              aicT      <- AIC(b13) + AIC(b23) + AIC(b12.3)             
              bicT      <- BIC(b13) + BIC(b23) + BIC(b12.3)
              
              cop.sel   <- c(b13$copula, b23$copula, b12.3$copula)
              dof.vec   <- c(b13$dof, b23$dof, b12.3$dof)
 
              ros.resids <- cbind(obj3$p, data$h.fm, third.res)
 
 
                 }                   
  
  
L <- list(margins = c(obj1$margins[1], obj2$margins[1], obj3$margins[1]), infl.fac = infl.fac,
     m1 = obj1, m2 = obj2, m3 = obj3, 
     b12 = b12, b13 = b13, b23.1 = b23.1, b23 = b23, b13.2 = b13.2, b12.3 = b12.3, 
     AIC.pcc = aicT, BIC.pcc = bicT, cond.var = cond.var, copula.sel = cop.sel, 
     dof.vec = dof.vec, simplified = simplified, eq.theta = eq.theta,
     left.trunc.v = c(obj1$VC$left.trunc, obj2$VC$left.trunc, obj3$VC$left.trunc), ros.resids = ros.resids)


class(L) <- c("gjrm.pcc")


L 
 
}
