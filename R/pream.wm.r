pream.wm <- function(formula, margins, M, l.flist, type = "copR"){
  
  #scs <- c("WEI", "FISK", "LN", "LO", "N", "N2") # for survival models
  scs <- c("WEI", "FISK") # for survival models

  
if(type == "biv"){


  if(M$BivD == "T" && (M$dof <=2 || M$dof > 249)) stop("dof must be a number greater than 2 and smaller than 249.")

  if(M$intf == FALSE && margins[2] %in% c(M$m1d,M$m2d,M$m2,M$m3)) stop("Please use copulaReg() for models involving binary and continuous/discrete margins.")
  
  if(!is.null(M$theta.fx) && M$BivD != "N") stop("This approach is not currently implemented for non-Gaussian bivariate distributions.")
  if(!is.null(M$theta.fx)) { if(M$Model != "B" && !(margins[2] %in% M$bl)) stop("Only bivariate Gaussian binary models with fixed theta are currenlty allowed for.")}
  
  if(!is.null(M$theta.fx) ) { if( abs(M$theta.fx) > 0.999 ) stop("The theta value must be in the interval [-0.999,0.999].") }

  if(M$Model == "BPO" && M$BivD != "N") stop("This model is not defined for copulae.")
  if(M$Model == "BPO" && margins[1] != "probit" && margins[2] != "probit") stop("This model is not defined for the chosen margins.")

  if(!(M$Model %in% M$mb)) stop("Error in parameter Model value. It should be one of:\nB, BSS, BPO, BPO0.")
  
  if( M$Model == "BSS" && M$BivD %in% M$BivD2 ) stop("Mixed copulae can not be implemented for selection models.")
  
  
  if(!(M$BivD %in% c(M$opc,M$BivD2))) stop("Error in parameter BivD value. It should be one of:\nN, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F, AMH, FGM, T, PL, HO,\nC0C90, C0C270, C180C90, C180C270, J0J90, J0J270, J180J90, J180J270,\nG0G90, G0G270, G180G90, G180G270.")
  
  
  if(!(M$extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of:\nt, pC or sED.")
  
  if(!(margins[1] %in% M$bl) ) stop("Error in first margin value. It should be one of:\nprobit, logit, cloglog.")
  if(!(margins[2] %in% c(M$bl,M$m2,M$m3)) && M$intf == FALSE ) stop("Error in second margin value. It should be one of:\nprobit, logit, cloglog.")  
  if(!(margins[2] %in% c(M$bl,M$m2,M$m3,M$m1d,M$m2d)) && M$intf == TRUE ) stop("Error in second margin value. It should be one of:\nN, N2, GU, rGU, LO, LN, WEI, iG, GA, DAGUM, SM, BE, FISK, PO, ZTP, NBI, NBII, PIG.")  

  if(margins[2] %in% c(M$m2,M$m3,M$m1d,M$m2d) && (M$Model == "BPO" || M$Model == "BPO0") ) stop("For continuous/discrete responses, partial observability models\nare not allowed for when using this function.")   
  
  if(margins[2] %in% c(M$m2,M$m3,M$m1d,M$m2d) && M$Model == "BSS" ) stop("Please use function copulaSampleSel().")   
    
  if(l.flist > 2 && margins[2] %in% c(M$bl,M$m1d)){ if(l.flist!=3) stop("You need to specify three equations.") } 
  if(l.flist > 2 && margins[2] %in% c(M$m2,M$m2d)){ if(l.flist!=4) stop("You need to specify four equations.") } 
  if(l.flist > 2 && margins[2] %in% M$m3         ){ if(l.flist!=5) stop("You need to specify five equations.") }  
  
  if( l.flist > 2  && M$Model == "BPO0")                 stop("You only need to specify two equations.\nThe chosen model does not have a correlation parameter.")
  if( l.flist > 2  && M$Model == "B" && !is.null(M$theta.fx)) stop("You only need to specify two equations.\nThe chosen model is not allowed to estimate the theta parameter.")
  


#####
  if(margins[2] %in% c(M$m1d,M$m2d)) stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")






}
  
  
  
  
  
  
  
if(type == "triv"){  
  

    mmar <- c("probit", "logit", "cloglog")

    if(!(M$Model %in% M$mb)) stop("Error in parameter Model value. It should be one of: T, TSS, TESS.")
      
      
    if(!(margins[1] %in% mmar) ) stop("Error in first margin value. It should be one of:\nprobit, logit, cloglog.")  
    if(!(margins[2] %in% mmar) ) stop("Error in second margin value. It should be one of:\nprobit, logit, cloglog.") 
    if(!(margins[3] %in% mmar) ) stop("Error in third margin value. It should be one of:\nprobit, logit, cloglog.")     
      
      
    if(!(M$penCor %in% c("unpen", "ridge", "lasso", "alasso"))) stop("Error in parameter penCor value. It should be one of:\nunpen, ridge, lasso, alasso.")
  
    if(is.null(M$w.alasso) && M$penCor == "alasso") stop("You must provide a vector of three weights when using alasso.")
  
    if(length(M$w.alasso)!=3 && M$penCor == "alasso") stop("You must supply a vector made up of three weights.")
  
    if(!(M$extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of:\nt, pC or sED.")
    
    #if(l.flist > 3) stop("You are only allowed to specify three equations.") 
    
    if(l.flist > 3 && M$penCor %in% c("ridge", "lasso", "alasso")) stop("This option is not currently allowed.") 

    if(l.flist > 3 && M$Chol == FALSE) stop("You must set Chol = TRUE.")
    
    if(l.flist > 3 && M$Model %in% c("TSS", "TESS")) stop("You can not currently use more than three equations with the model chosen.") 


    
    
    #####
  if(M$Model == "TSS") stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")
  
  
  
}  
  
  
if(type == "copR"){


if(M$surv == TRUE){


  if( !(margins[1] %in% c(M$m2,M$m3,M$bl)) ) stop("The first marginal distribution must be either continuous or probit, PO or PH.")
  if( !(margins[2] %in% c(M$bl,scs)) )       stop("The second marginal distribution must be probit, PO or PH.")

  #if(!(margins[1] %in% c(scs,M$bl)) || !(margins[2] %in% c(scs,M$bl))  ) stop("The marginal distributions must be WEI, FISK, probit, PO or PH.")
  #if( (margins[1] %in% scs && margins[2] %in% M$bl) || (margins[1] %in% M$bl && margins[2] %in% scs) ) stop("You can not mix semiparametric and parametric margins.")


  if(!is.null(M$c1) && length(table( M$c1 %in% c(0,1) ) ) > 1 ) stop("Your first censoring indicator is not binary. Please fix.")
  if(!is.null(M$c2) && length(table( M$c2 %in% c(0,1) ) ) > 1 ) stop("Your second censoring indicator is not binary. Please fix.")
  
  if(length(M$c1) != length(M$c2)) stop("The two censoring indicators must have the same length.")

  if(l.flist > 2 && margins[1] %in% c(M$bl) && margins[2] %in% c(M$bl)){ if(l.flist!=3) stop("You need to specify three equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m2) && margins[2] %in% c(M$bl)){ if(l.flist!=4) stop("You need to specify four equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m3) && margins[2] %in% c(M$bl)){ if(l.flist!=5) stop("You need to specify five equations.") } 



   
}  
  
  if(M$BivD == "T" && (M$dof <=2 || M$dof > 249)) stop("dof must be a number greater than 2 and smaller than 249.")

  if( margins[1] %in% c(M$m2d) && margins[2] %in% c(M$m1d) ) stop("Please swap the two equations (and hence margins' specification).\nThe second instead of the first margin has to be a two-parameter discrete distribution.")
  if( margins[1] %in% c(M$m2,M$m3) && margins[2] %in% c(M$m1d,M$m2d) ) stop("Please swap the two equations (and hence margins' specification).\nThe first instead of the second margin has to be discrete.")

  if(!(M$BivD %in% c(M$opc,M$BivD2))) stop("Error in parameter BivD value. It should be one of:\nN, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F, AMH, FGM, T, PL, HO, \nC0C90, C0C270, C180C90, C180C270, J0J90, J0J270, J180J90, J180J270,\nG0G90, G0G270, G180G90, G180G270.")
 
  if(!(M$extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of:\nt, pC or sED.")
  
  if(!(margins[1] %in% c(M$m2,M$m3,M$m1d,M$m2d,M$bl)) ) stop("Error in first margin value. It should be one of:\nN, N2, GU, rGU, LO, LN, WEI, iG, GA, DAGUM, SM, BE, FISK, NBI, NBII, PIG, PO, ZTP\nor probit, logit, cloglog for survival models.")  
  if(!(margins[2] %in% c(M$m2,M$m3,M$m1d,M$m2d,M$bl)) ) stop("Error in second margin value. It should be one of:\nN, N2, GU, rGU, LO, LN, WEI, iG, GA, DAGUM, SM, BE, FISK, NBI, NBII, PIG, PO, ZTP\nor probit, logit, cloglog for survival models.")  
  

  
  

 
  if(M$BivD == "T" && margins[1] %in% c(M$m2,M$m3) && margins[2] %in% c(M$m2,M$m3)){   
  
  if(l.flist > 2 && margins[1] %in% c(M$m2) && margins[2] %in% c(M$m2)){ if(l.flist!=6) stop("You need to specify six equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m2) && margins[2] %in% c(M$m3)){ if(l.flist!=7) stop("You need to specify seven equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m3) && margins[2] %in% c(M$m2)){ if(l.flist!=7) stop("You need to specify seven equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m3) && margins[2] %in% c(M$m3)){ if(l.flist!=8) stop("You need to specify eight equations.") } 
  
  }else{
  
  
  if(l.flist > 2 && margins[1] %in% c(M$m1d,M$bl)    && margins[2] %in% c(M$m1d,M$bl))   { if(l.flist!=3) stop("You need to specify three equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m2d)    && margins[2] %in% c(M$m1d))             { if(l.flist!=4) stop("You need to specify four equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m1d)    && margins[2] %in% c(M$m2,M$m2d))        { if(l.flist!=4) stop("You need to specify four equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m1d)    && margins[2] %in% c(M$m3,M$m3d))        { if(l.flist!=5) stop("You need to specify five equations.") }   
  
  if(l.flist > 2 && margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d)){ if(l.flist!=5) stop("You need to specify five equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m3,M$m3d)){ if(l.flist!=6) stop("You need to specify six equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m3,M$m3d) && margins[2] %in% c(M$m2,M$m2d)){ if(l.flist!=6) stop("You need to specify six equations.") } 
  if(l.flist > 2 && margins[1] %in% c(M$m3,M$m3d) && margins[2] %in% c(M$m3,M$m3d)){ if(l.flist!=7) stop("You need to specify seven equations.") } 
  
  }
 
 
 
 ######
 if(margins[1] %in% c(M$m1d,M$m2d) && margins[2] %in% c(M$m1d,M$m2d,M$m2,M$m3) ) stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")
 if( M$surv == TRUE ) stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")
 if( M$surv == TRUE && margins[1] %in% c(M$m2,M$m3) && margins[2] %in% c(M$bl) ) stop("Survival models with first continuous margin and second survival one not ready yet.\nGet in touch to check progress.")


}


if(type == "copSS"){

  if(M$BivD == "T" && (M$dof <=2 || M$dof > 249)) stop("dof must be a number greater than 2 and smaller than 249.")

    
  if(!(M$BivD %in% M$opc)) stop("Error in parameter BivD value. It should be one of: N, C0, C90, C180, C270, J0, J90, J180, J270, G0, G90, G180, G270, F, AMH, FGM, T, PL, HO.")
  if(!(M$extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of: t, pC or sED.")
  
  if(!(M$margins[1] %in% M$bl) ) stop("Error in first margin value. It should be one of:\nprobit, logit, cloglog.")
  if(!(M$margins[2] %in% c(M$m2,M$m3,M$m1d,M$m2d)) ) stop("Error in second margin value. It should be one of:\nN, N2, GU, rGU, LO, LN, WEI, iG, GA, DAGUM, SM, BE, FISK, PO, ZTP, NBI, NBII, PIG.")  
  
  if(l.flist > 2 && M$margins[2] %in% c(M$m1d))   { if(l.flist!=3) stop("You need to specify three equations.") } 
  if(l.flist > 2 && M$margins[2] %in% c(M$m2,M$m2d)){ if(l.flist!=4) stop("You need to specify four equations.") } 
  if(l.flist > 2 && M$margins[2] %in% M$m3       ){ if(l.flist!=5) stop("You need to specify five equations.") }  
 
 
 #####
 if(margins[2] %in% c("GU", "rGU", "LO", "LN", "WEI","iG", "DAGUM", "SM", "BE", "FISK") ) stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")
 

}







if(type == "gamls"){


if(M$surv == TRUE){

  if( !(M$margin %in% c(scs,M$bl)) ) stop("The marginal distribution can be WEI, FISK, probit, PO or PH.")
  if(is.null(M$cens) ) stop("You must provide the binary censoring indicator.")
  
  if(!is.null(M$cens) && length(table( M$cens %in% c(0,1) ) ) > 1 ) stop("Your censoring indicator is not binary. Please fix.")
     
}  

  if(M$surv == TRUE && M$robust == TRUE) stop("It is not currently possible to fit robust survival models.")


  
  if(!(M$extra.regI %in% c("t","pC","sED"))) stop("Error in parameter extra.regI value. It should be one of:\nt, pC or sED.")
  if(!(M$margin %in% c(M$m2,M$m3,M$m1d,M$m2d, M$bl)) ) stop("Error in margin value. It should be one of:\nN, N2, GU, rGU, LO, LN, WEI, iG, GA, DAGUM, SM, BE, FISK, NBI, NBII, PIG, PO, ZTP, GEVlink.")  
  if(l.flist > 1 && M$margin %in% c(M$m1d)    )                   stop("You need to specify one equation.")  
  if(l.flist > 1 && M$margin %in% c(M$m2,M$m2d) ){ if(l.flist!=2) stop("You need to specify two equations.")   } 
  if(l.flist > 1 && M$margin %in% c(M$m3,M$m3d) ){ if(l.flist!=3) stop("You need to specify three equations.") } 
  
 if(M$robust == TRUE) stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")
 if( M$surv == TRUE ) stop("Check the next release for the final tested version of this model\nor get in touch to check progress.")

}
  
  
  
  
  
  
  
 
}




















