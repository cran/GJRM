func.OPT <- function(margins, M, type = "copR"){
  
if(type == "ROY"){

  if(margins[2] %in% M$m1d && margins[3] %in% M$m1d) func.opt <- bprobgHsDiscr1ROY
  if(margins[2] %in% M$m2d && margins[3] %in% M$m2d) func.opt <- bprobgHsDiscr2ROY
  
  if(margins[2] %in% M$m2  && margins[3] %in% M$m2 ) func.opt <- bprobgHsCont2ROY
  if(margins[2] %in% M$m3  && margins[3] %in% M$m3 ) func.opt <- bprobgHsCont3ROY
  
  if(margins[2] %in% M$bl  && margins[3] %in% M$bl ) func.opt <- bprobgHsBinROY  

  
  # mixed cases like m1 and m2 not done at the moment as they may not be plausible empirically
    
  
}  
  
  

if(type == "biv"){

  if(M$Model=="B" && margins[2] %in% M$bl  &&  is.null(M$theta.fx)) func.opt <- bprobgHs   
  if(M$Model=="B" && margins[2] %in% M$bl  && !is.null(M$theta.fx)) func.opt <- bprobgHsFixTheta

  if(M$Model=="B" && margins[2] %in% M$m1d  ) {func.opt <- bprobgHsDiscr1 ; func.optUniv <- bprobgHsContUniv}   
  if(M$Model=="B" && margins[2] %in% M$m2d  ) {func.opt <- bprobgHsDiscr2 ; func.optUniv <- bprobgHsContUniv} 
  
  if(M$Model=="BPO")                        func.opt <- bprobgHsPO 
  if(M$Model=="BPO0")                       func.opt <- bprobgHsPO0   
  if(M$Model=="BSS")                        func.opt <- bprobgHsSS   
  
  if(M$Model=="B" && margins[2] %in% M$m2 &&  is.null(M$K1) ) {func.opt <- bprobgHsCont  ; func.optUniv <- bprobgHsContUniv }   
  if(M$Model=="B" && margins[2] %in% M$m3 &&  is.null(M$K1) ) {func.opt <- bprobgHsCont3 ; func.optUniv <- bprobgHsContUniv3}  
  
  if(M$Model=="B" && margins[2] == "TW"   &&  is.null(M$K1) ) {func.optUniv <- bprobgHsContUniv3}  

  if(M$Model=="B" && margins[2] %in% M$m2 && !is.null(M$K1) ) {func.opt <- bCopulaCLMgHsCont}
  
}
  
  
  
if(type == "copR"){


  
  if(margins[1] %in% M$m1d && margins[2] %in% M$m1d) func.opt  <- bdiscrdiscr11 
  if(margins[1] %in% M$m1d && margins[2] %in% M$m2d) func.opt  <- bdiscrdiscr12
  if(margins[1] %in% M$m2d && margins[2] %in% M$m2d) func.opt  <- bdiscrdiscr  
  
  
  if(margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$BivD != "T") {if(M$robust == FALSE) func.opt <- bcont else func.opt <- bcontROB } 
  if(margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$BivD == "T") func.opt <- bconttwoParC  

  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && M$BivD != "T") func.opt <- bcont3       
  if(margins[1] %in% M$m3 && margins[2] %in% M$m3 && M$BivD == "T") func.opt <- bcont3twoParC   
  
  if(margins[1] %in% M$m2 && margins[2] %in% M$m3 && M$BivD != "T") func.opt <- bcont23
  if(margins[1] %in% M$m2 && margins[2] %in% M$m3 && M$BivD == "T") func.opt <- bcont23twoParC
  
  if(margins[1] %in% M$m3 && margins[2] %in% M$m2 && M$BivD != "T") func.opt <- bcont32 
  if(margins[1] %in% M$m3 && margins[2] %in% M$m2 && M$BivD == "T") func.opt <- bcont32twoParC 
  
  
  
  
  
  
  if(margins[1] %in% M$m2 && margins[2] %in% M$m2 && M$surv == TRUE)  func.opt <- bcontSurv   # not really used
  
  if(margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == TRUE && M$dep.cens == FALSE)  func.opt <- bcontSurvG
  
  
  # *** NEW ***
  # New function - includes all types of censoring ***

  if(margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == TRUE && M$dep.cens == FALSE 
     && M$type.cens1 == 'mixed' && M$type.cens2 == 'mixed')  func.opt <- bcontSurvG_extended
  
  # **************************************************  
  
  
  
  
  

  
  if(margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == TRUE && M$dep.cens == TRUE &&  is.null(M$c3)) func.opt <- bcontSurvGDep
  if(margins[1] %in% M$bl && margins[2] %in% M$bl && M$surv == TRUE && M$dep.cens == TRUE && !is.null(M$c3)) func.opt <- bcontSurvGDepA
  

}


if(type == "copSS"){


  if(margins[2] %in% M$m2    ) func.opt <- bprobgHsContSS    
  if(margins[2] %in% M$m3    ) func.opt <- bprobgHsCont3SS  
  
  if(margins[2] %in% M$m1d   ) func.opt <- bprobgHsDiscr1SS
  if(margins[2] %in% M$m2d   ) func.opt <- bprobgHsDiscr2SS  

  

}
  

func.opt


}

