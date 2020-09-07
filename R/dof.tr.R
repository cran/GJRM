dof.tr <- function(var.st){
  
   var.st <- ifelse( var.st > 5.51,   5.51, var.st ) # can't b over 250  
   var.st <- ifelse( var.st < -10, -10, var.st ) 
 
   vao <- exp(var.st) + 2
    
 list(var.st = var.st, vao = vao )  
 
}    