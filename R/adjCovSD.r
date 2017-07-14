adjCovSD <- function(x, design){

scores <- aCov(x)

    if (inherits(design,"survey.design2"))
      covsan <- svyrecvar(scores%*%x$Vb,design$cluster,design$strata,design$fpc,postStrata=design$postStrata)
    else if (inherits(design, "twophase"))
      covsan <- twophasevar(scores%*%x$Vb, design)
    else if (inherits(design, "twophase2"))
      covsan <- twophase2var(scores%*%x$Vb, design)
    else
      covsan <- svyCprod(scores%*%x$Vb,design$strata,design$cluster[[1]],design$fpc, design$nPSU,
                  design$certainty,design$postStrata)
                 
  x$Vb <- covsan
  
  rm(covsan, scores)
  
  x              
                                             
  }
  
  