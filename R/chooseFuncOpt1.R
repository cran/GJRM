chooseFuncOpt1 <- function(margin, m1d, m2d, m2, m3, bl, informative, type.cens, hrate, truncation.time){
    
  if (margin %in% c(m1d, m2d, m2)) func.opt1 <- bprobgHsContUniv
  if (margin %in% c(m3)) func.opt1 <- bprobgHsContUniv3
  if (margin %in% c(bl) && informative == "no" && type.cens == "R" && is.null(hrate)) func.opt1 <- bcontSurvGuniv
  if (margin %in% c(bl) && informative == "no" && type.cens == "L" && is.null(hrate)) func.opt1 <- bcontSurvGunivL
  if (margin %in% c(bl) && informative == "no" && type.cens == "I" && is.null(hrate)) func.opt1 <- bcontSurvGunivI
  if (margin %in% c(bl) && informative == "no" && type.cens == "mixed" && is.null(hrate)) func.opt1 <- bcontSurvGunivMIXED
  if (margin %in% c(bl) && informative == "yes") func.opt1 <- bcontSurvGunivInform 
  if (margin %in% c(bl) && informative == "no" && type.cens == "R" && !is.null(hrate)) func.opt1 <- bcontSurvGuniv_ExcessHazard
  if (margin %in% c(bl) && informative == "no" && type.cens == "L" && !is.null(hrate)) func.opt1 <- bcontSurvGunivL_ExcessHazard
  if (margin %in% c(bl) && informative == "no" && type.cens == "I" && !is.null(hrate)) func.opt1 <- bcontSurvGunivI_ExcessHazard
  if (margin %in% c(bl) && informative == "no" && type.cens == "mixed" && !is.null(hrate)) func.opt1 <- bcontSurvGunivMIXED_ExcessHazard
  if (margin %in% c(bl) && informative == "no" && type.cens == "mixed" && is.null(hrate) && !is.null(truncation.time)) func.opt1 <- bcontSurvGunivMIXED_LeftTruncation
  if (margin %in% c(bl) && informative == "no" && type.cens == "mixed" && !is.null(hrate) && !is.null(truncation.time)) func.opt1 <- bcontSurvGunivMIXED_ExcessHazard_LeftTruncation
  
  return(func.opt1)
}