pp <- function(x){

   cont1par <- c("PO","ZTP","probit","logit","cloglog","GEVlink")  
   cont2par <- c("N","N2","GU","rGU","LO","LN","WEI","GO","iG","GA","GA2","GAi","BE","FISK","NBI", "NBII","NBIa", "NBIIa","PIG","GP","GPII","GPo","DGP","DGPII")  
   cont3par <- c("DAGUM","SM","DEL","SICHEL","GGA","TW")   

if(x$univar.gamlss == TRUE) x$BivD <- "N" 

  if(x$BivD=="FGM")  {cop <- "FGM"                ;lind <- "atanh"} 
  if(x$BivD=="T")    {cop <- paste("Student-t (dof = ",format(x$dof.a, digits=3),")",sep="")   ;lind <- "atanh"} 
  if(x$BivD=="AMH")  {cop <- "AMH"                ;lind <- "atanh"} 
  if(x$BivD=="N")    {cop <- "Gaussian"           ;lind <- "atanh"} 
  if(x$BivD=="F")    {cop <- "Frank"              ;lind <- "identity"}       
  
  if(x$BivD=="PL")   {cop <- "Plackett"           ;lind <- "log"}
  if(x$BivD=="HO")   {cop <- "Hougaard"           ;lind <- "qlogis"} 
  
  if(x$BivD=="C0")   {cop <- "Clayton"            ;lind <- "log"}   
  if(x$BivD=="C90")  {cop <- "90\u00B0 Clayton"   ;lind <- "log(- \u00B7)"}                 
  if(x$BivD=="C180") {cop <- "180\u00B0 Clayton"  ;lind <- "log"}                    
  if(x$BivD=="C270") {cop <- "270\u00B0 Clayton"  ;lind <- "log(- \u00B7)"}    
  if(x$BivD=="J0")   {cop <- "Joe"                ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="J90")  {cop <- "90\u00B0 Joe"       ;lind <- "log(- \u00B7 - 1)"}
  if(x$BivD=="J180") {cop <- "180\u00B0 Joe"      ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="J270") {cop <- "270\u00B0 Joe"      ;lind <- "log(- \u00B7 - 1)"}
  if(x$BivD=="G0")   {cop <- "Gumbel"             ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="G90")  {cop <- "90\u00B0 Gumbel"    ;lind <- "log(- \u00B7 - 1)"}
  if(x$BivD=="G180") {cop <- "180\u00B0 Gumbel"   ;lind <- "log(\u00B7 - 1)"} 
  if(x$BivD=="G270") {cop <- "270\u00B0 Gumbel"   ;lind <- "log(- \u00B7 - 1)"} 
  
  if(x$BivD=="C0C90")    {cop <- "Clayton & 90\u00B0 Clayton"            ;lind <- "log & log(- \u00B7)"}   
  if(x$BivD=="C0C270")   {cop <- "Clayton & 270\u00B0 Clayton"           ;lind <- "log & log(- \u00B7)"}                 
  if(x$BivD=="C180C90")  {cop <- "180\u00B0 Clayton & 90\u00B0 Clayton"  ;lind <- "log & log(- \u00B7)"}                    
  if(x$BivD=="C180C270") {cop <- "180\u00B0 Clayton & 270\u00B0 Clayton" ;lind <- "log & log(- \u00B7)"}  
 
  if(x$BivD=="J0J90")    {cop <- "Joe & 90\u00B0 Joe"            ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}   
  if(x$BivD=="J0J270")   {cop <- "Joe & 270\u00B0 Joe"           ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                 
  if(x$BivD=="J180J90")  {cop <- "180\u00B0 Joe & 90\u00B0 Joe"  ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                    
  if(x$BivD=="J180J270") {cop <- "180\u00B0 Joe & 270\u00B0 Joe" ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}  
  
  if(x$BivD=="G0G90")    {cop <- "Gumbel & 90\u00B0 Gumbel"            ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}   
  if(x$BivD=="G0G270")   {cop <- "Gumbel & 270\u00B0 Gumbel"           ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                 
  if(x$BivD=="G180G90")  {cop <- "180\u00B0 Gumbel & 90\u00B0 Gumbel"  ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                    
  if(x$BivD=="G180G270") {cop <- "180\u00B0 Gumbel & 270\u00B0 Gumbel" ;lind <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}    
 

  
  # if(x$BivD=="N" && x$Model=="BPO0") cop <- "Independent"   
    
  mml <- c("LN","WEI","GO","iG","GA","GA2","GGA","DAGUM","SM","FISK","NBI","NBII","NBIa","NBIIa","PIG","PO","ZTP","GP","GPII","GPo","DGP","DGPII","TW")  
  
  
  
    
  if(x$margins[1] %in% c("N","N2","GU","rGU","LO","GAi") )  m1l <- "identity"
  if(x$margins[2] %in% c("N","N2","GU","rGU","LO","GAi") )  m2l <- "identity"
  
  if(x$margins[1] %in% mml )                                m1l <- "log" 
  if(x$margins[2] %in% mml )                                m2l <- "log" 
  
  if(x$margins[1] %in% c("GP","DGP") )                    m1l <- "identity"
  if(x$margins[2] %in% c("GP","DGP") )                    m2l <- "identity"  
  
  if(x$margins[1] %in% c("DGPII") )                       m1l <- "sqrt"
  if(x$margins[2] %in% c("DGPII") )                       m2l <- "sqrt"  
  
  if(x$margins[1] %in% c("GPII","GPo") )                       m1l <- "log(\u00B7 + 0.5)"
  if(x$margins[2] %in% c("GPII","GPo") )                       m2l <- "log(\u00B7 + 0.5)"     
  
  if(x$margins[1] %in% c("BE") )                            m1l <- "qlogis" 
  if(x$margins[2] %in% c("BE") )                            m2l <- "qlogis"   
 
  if(x$margins[1]=="probit")  m1l <- "probit"
  if(x$margins[1]=="logit")   m1l <- "logit"
  if(x$margins[1]=="cloglog") m1l <- "cloglog"
  if(x$margins[1]=="cauchit") m1l <- "cauchit"
  if(x$margins[1]=="GEVlink") m1l <- "GEVlink" 
  
  
  if(x$margins[2]=="probit")  m2l <- "probit"
  if(x$margins[2]=="logit")   m2l <- "logit"
  if(x$margins[2]=="cloglog") m2l <- "cloglog"
  if(x$margins[2]=="cauchit") m2l <- "cauchit"  
  if(x$margins[2]=="GEVlink") m2l <- "GEVlink" 
  
 
 
list(cont1par = cont1par, cont2par = cont2par, cont3par = cont3par, cop = cop, lind = lind, m1l = m1l, m2l = m2l)
  

}

