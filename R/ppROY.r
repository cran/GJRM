ppROY <- function(x){

  cont1par <- c("P","tP","DGP0","probit","logit","cloglog","GEVlink")  
  cont2par <- c("N","N2","GU","rGU","LO","LN","WEI","GO","IG","GA","GA2","GAi","BE","FISK","tNBI", "tNBII","tPIG","NBI", "NBII","PIG","GP","GPII","GPo","DGP","DGPII")  
  cont3par <- c("DAGUM","SM","DEL","SICHEL","GGA","TW")   

  if(x$BivD1=="FGM")          {cop1 <- "FGM"                                     ;lind1 <- "atanh"} 
  if(x$BivD1=="T")            {cop1 <- "Student-t"                               ;lind <- "atanh"} # {cop1 <- paste("Student-t (dof = ",format(x$dof12.a, digits=3),")",sep=""); lind1 <- "atanh"} 
  if(x$BivD1=="AMH")          {cop1 <- "AMH"                                     ;lind1 <- "atanh"} 
  if(x$BivD1=="N")            {cop1 <- "Gaussian"                                ;lind1 <- "atanh"} 
  if(x$BivD1=="F")            {cop1 <- "Frank"                                   ;lind1 <- "identity"}       
  if(x$BivD1=="PL")           {cop1 <- "Plackett"                                ;lind1 <- "log"}
  if(x$BivD1=="HO")           {cop1 <- "Hougaard"                                ;lind1 <- "qlogis"} 
  if(x$BivD1=="C0")           {cop1 <- "Clayton"                                 ;lind1 <- "log"}   
  if(x$BivD1=="C90")          {cop1 <- "90\u00B0 Clayton"                        ;lind1 <- "log(- \u00B7)"}                 
  if(x$BivD1=="C180")         {cop1 <- "180\u00B0 Clayton"                       ;lind1 <- "log"}                    
  if(x$BivD1=="C270")         {cop1 <- "270\u00B0 Clayton"                       ;lind1 <- "log(- \u00B7)"}    
  if(x$BivD1=="GAL0")         {cop1 <- "Galambos"                                ;lind1 <- "log"}   
  if(x$BivD1=="GAL90")        {cop1 <- "90\u00B0 Galambos"                       ;lind1 <- "log(- \u00B7)"}                 
  if(x$BivD1=="GAL180")       {cop1 <- "180\u00B0 Galambos"                      ;lind1 <- "log"}                    
  if(x$BivD1=="GAL270")       {cop1 <- "270\u00B0 Galambos"                      ;lind1 <- "log(- \u00B7)"}    
  if(x$BivD1=="J0")           {cop1 <- "Joe"                                     ;lind1 <- "log(\u00B7 - 1)"} 
  if(x$BivD1=="J90")          {cop1 <- "90\u00B0 Joe"                            ;lind1 <- "log(- \u00B7 - 1)"}
  if(x$BivD1=="J180")         {cop1 <- "180\u00B0 Joe"                           ;lind1 <- "log(\u00B7 - 1)"} 
  if(x$BivD1=="J270")         {cop1 <- "270\u00B0 Joe"                           ;lind1 <- "log(- \u00B7 - 1)"}
  if(x$BivD1=="G0")           {cop1 <- "Gumbel"                                  ;lind1 <- "log(\u00B7 - 1)"} 
  if(x$BivD1=="G90")          {cop1 <- "90\u00B0 Gumbel"                         ;lind1 <- "log(- \u00B7 - 1)"}
  if(x$BivD1=="G180")         {cop1 <- "180\u00B0 Gumbel"                        ;lind1 <- "log(\u00B7 - 1)"} 
  if(x$BivD1=="G270")         {cop1 <- "270\u00B0 Gumbel"                        ;lind1 <- "log(- \u00B7 - 1)"} 
  if(x$BivD1=="C0C90")        {cop1 <- "Clayton & 90\u00B0 Clayton"              ;lind1 <- "log & log(- \u00B7)"}   
  if(x$BivD1=="C0C270")       {cop1 <- "Clayton & 270\u00B0 Clayton"             ;lind1 <- "log & log(- \u00B7)"}                 
  if(x$BivD1=="C180C90")      {cop1 <- "180\u00B0 Clayton & 90\u00B0 Clayton"    ;lind1 <- "log & log(- \u00B7)"}                    
  if(x$BivD1=="C180C270")     {cop1 <- "180\u00B0 Clayton & 270\u00B0 Clayton"   ;lind1 <- "log & log(- \u00B7)"} 
  if(x$BivD1=="GAL0GAL90")    {cop1 <- "Galambos & 90\u00B0 Galambos"            ;lind1 <- "log & log(- \u00B7)"}   
  if(x$BivD1=="GAL0GAL270")   {cop1 <- "Galambos & 270\u00B0 Galambos"           ;lind1 <- "log & log(- \u00B7)"}                 
  if(x$BivD1=="GAL180GAL90")  {cop1 <- "180\u00B0 Galambos & 90\u00B0 Galambos"  ;lind1 <- "log & log(- \u00B7)"}                    
  if(x$BivD1=="GAL180GAL270") {cop1 <- "180\u00B0 Galambos & 270\u00B0 Galambos" ;lind1 <- "log & log(- \u00B7)"}    
  if(x$BivD1=="J0J90")        {cop1 <- "Joe & 90\u00B0 Joe"                      ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}   
  if(x$BivD1=="J0J270")       {cop1 <- "Joe & 270\u00B0 Joe"                     ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                 
  if(x$BivD1=="J180J90")      {cop1 <- "180\u00B0 Joe & 90\u00B0 Joe"            ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                    
  if(x$BivD1=="J180J270")     {cop1 <- "180\u00B0 Joe & 270\u00B0 Joe"           ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}  
  if(x$BivD1=="G0G90")        {cop1 <- "Gumbel & 90\u00B0 Gumbel"                ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}   
  if(x$BivD1=="G0G270")       {cop1 <- "Gumbel & 270\u00B0 Gumbel"               ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                 
  if(x$BivD1=="G180G90")      {cop1 <- "180\u00B0 Gumbel & 90\u00B0 Gumbel"      ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                    
  if(x$BivD1=="G180G270")     {cop1 <- "180\u00B0 Gumbel & 270\u00B0 Gumbel"     ;lind1 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}    
 
  if(x$BivD2=="FGM")          {cop2 <- "FGM"                                     ;lind2 <- "atanh"} 
  if(x$BivD2=="T")            {cop2 <- "Student-t"                               ;lind <- "atanh"} # {cop2 <- paste("Student-t (dof = ",format(x$dof13.a, digits=3),")",sep=""); lind2 <- "atanh"} 
  if(x$BivD2=="AMH")          {cop2 <- "AMH"                                     ;lind2 <- "atanh"} 
  if(x$BivD2=="N")            {cop2 <- "Gaussian"                                ;lind2 <- "atanh"} 
  if(x$BivD2=="F")            {cop2 <- "Frank"                                   ;lind2 <- "identity"}       
  if(x$BivD2=="PL")           {cop2 <- "Plackett"                                ;lind2 <- "log"}
  if(x$BivD2=="HO")           {cop2 <- "Hougaard"                                ;lind2 <- "qlogis"} 
  if(x$BivD2=="C0")           {cop2 <- "Clayton"                                 ;lind2 <- "log"}   
  if(x$BivD2=="C90")          {cop2 <- "90\u00B0 Clayton"                        ;lind2 <- "log(- \u00B7)"}                 
  if(x$BivD2=="C180")         {cop2 <- "180\u00B0 Clayton"                       ;lind2 <- "log"}                    
  if(x$BivD2=="C270")         {cop2 <- "270\u00B0 Clayton"                       ;lind2 <- "log(- \u00B7)"}    
  if(x$BivD2=="GAL0")         {cop2 <- "Galambos"                                ;lind2 <- "log"}   
  if(x$BivD2=="GAL90")        {cop2 <- "90\u00B0 Galambos"                       ;lind2 <- "log(- \u00B7)"}                 
  if(x$BivD2=="GAL180")       {cop2 <- "180\u00B0 Galambos"                      ;lind2 <- "log"}                    
  if(x$BivD2=="GAL270")       {cop2 <- "270\u00B0 Galambos"                      ;lind2 <- "log(- \u00B7)"}    
  if(x$BivD2=="J0")           {cop2 <- "Joe"                                     ;lind2 <- "log(\u00B7 - 1)"} 
  if(x$BivD2=="J90")          {cop2 <- "90\u00B0 Joe"                            ;lind2 <- "log(- \u00B7 - 1)"}
  if(x$BivD2=="J180")         {cop2 <- "180\u00B0 Joe"                           ;lind2 <- "log(\u00B7 - 1)"} 
  if(x$BivD2=="J270")         {cop2 <- "270\u00B0 Joe"                           ;lind2 <- "log(- \u00B7 - 1)"}
  if(x$BivD2=="G0")           {cop2 <- "Gumbel"                                  ;lind2 <- "log(\u00B7 - 1)"} 
  if(x$BivD2=="G90")          {cop2 <- "90\u00B0 Gumbel"                         ;lind2 <- "log(- \u00B7 - 1)"}
  if(x$BivD2=="G180")         {cop2 <- "180\u00B0 Gumbel"                        ;lind2 <- "log(\u00B7 - 1)"} 
  if(x$BivD2=="G270")         {cop2 <- "270\u00B0 Gumbel"                        ;lind2 <- "log(- \u00B7 - 1)"} 
  if(x$BivD2=="C0C90")        {cop2 <- "Clayton & 90\u00B0 Clayton"              ;lind2 <- "log & log(- \u00B7)"}   
  if(x$BivD2=="C0C270")       {cop2 <- "Clayton & 270\u00B0 Clayton"             ;lind2 <- "log & log(- \u00B7)"}                 
  if(x$BivD2=="C180C90")      {cop2 <- "180\u00B0 Clayton & 90\u00B0 Clayton"    ;lind2 <- "log & log(- \u00B7)"}                    
  if(x$BivD2=="C180C270")     {cop2 <- "180\u00B0 Clayton & 270\u00B0 Clayton"   ;lind2 <- "log & log(- \u00B7)"} 
  if(x$BivD2=="GAL0GAL90")    {cop2 <- "Galambos & 90\u00B0 Galambos"            ;lind2 <- "log & log(- \u00B7)"}   
  if(x$BivD2=="GAL0GAL270")   {cop2 <- "Galambos & 270\u00B0 Galambos"           ;lind2 <- "log & log(- \u00B7)"}                 
  if(x$BivD2=="GAL180GAL90")  {cop2 <- "180\u00B0 Galambos & 90\u00B0 Galambos"  ;lind2 <- "log & log(- \u00B7)"}                    
  if(x$BivD2=="GAL180GAL270") {cop2 <- "180\u00B0 Galambos & 270\u00B0 Galambos" ;lind2 <- "log & log(- \u00B7)"}    
  if(x$BivD2=="J0J90")        {cop2 <- "Joe & 90\u00B0 Joe"                      ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}   
  if(x$BivD2=="J0J270")       {cop2 <- "Joe & 270\u00B0 Joe"                     ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                 
  if(x$BivD2=="J180J90")      {cop2 <- "180\u00B0 Joe & 90\u00B0 Joe"            ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                    
  if(x$BivD2=="J180J270")     {cop2 <- "180\u00B0 Joe & 270\u00B0 Joe"           ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}  
  if(x$BivD2=="G0G90")        {cop2 <- "Gumbel & 90\u00B0 Gumbel"                ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}   
  if(x$BivD2=="G0G270")       {cop2 <- "Gumbel & 270\u00B0 Gumbel"               ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                 
  if(x$BivD2=="G180G90")      {cop2 <- "180\u00B0 Gumbel & 90\u00B0 Gumbel"      ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}                    
  if(x$BivD2=="G180G270")     {cop2 <- "180\u00B0 Gumbel & 270\u00B0 Gumbel"     ;lind2 <- "log(\u00B7 - 1) & log(- \u00B7 - 1)"}      
    
    
    
    
    
    
    
  mml <- c("LN","WEI","GO","IG","GA","GA2","GGA","DAGUM","SM","FISK","tNBI","tNBII","tPIG","NBI","NBII","PIG","P","DGP0","tP","GP","GPII","GPo","DGP","DGPII","TW")  
  
  
  
    
  if(x$margins[1] %in% c("N","N2","GU","rGU","LO","GAi") )  m1l <- "identity"
  if(x$margins[2] %in% c("N","N2","GU","rGU","LO","GAi") )  m2l <- "identity"
  if(x$margins[3] %in% c("N","N2","GU","rGU","LO","GAi") )  m3l <- "identity"

  if(x$margins[1] %in% mml )                                m1l <- "log" 
  if(x$margins[2] %in% mml )                                m2l <- "log" 
  if(x$margins[3] %in% mml )                                m3l <- "log" 
  
  if(x$margins[1] %in% c("GP","DGP") )                      m1l <- "identity"
  if(x$margins[2] %in% c("GP","DGP") )                      m2l <- "identity"  
  if(x$margins[3] %in% c("GP","DGP") )                      m3l <- "identity"  
 
  if(x$margins[1] %in% c("DGPII") )                         m1l <- "log" # "sqrt"
  if(x$margins[2] %in% c("DGPII") )                         m2l <- "log" # "sqrt"  
  if(x$margins[3] %in% c("DGPII") )                         m3l <- "log" # "sqrt"  
  
  if(x$margins[1] %in% c("GPII","GPo") )                    m1l <- "log(\u00B7 + 0.5)"
  if(x$margins[2] %in% c("GPII","GPo") )                    m2l <- "log(\u00B7 + 0.5)"    
  if(x$margins[3] %in% c("GPII","GPo") )                    m3l <- "log(\u00B7 + 0.5)"     
  
  if(x$margins[1] %in% c("BE") )                            m1l <- "qlogis" 
  if(x$margins[2] %in% c("BE") )                            m2l <- "qlogis"   
  if(x$margins[3] %in% c("BE") )                            m3l <- "qlogis"   
  

  if(x$margins[1] == "probit")  m1l <- "probit"
  if(x$margins[1] == "logit")   m1l <- "logit"
  if(x$margins[1] == "cloglog") m1l <- "cloglog"
  if(x$margins[1] == "cauchit") m1l <- "cauchit"
  if(x$margins[1] == "GEVlink") m1l <- "GEVlink" 
  
  if(x$margins[2] == "probit")  m2l <- "probit"
  if(x$margins[2] == "logit")   m2l <- "logit"
  if(x$margins[2] == "cloglog") m2l <- "cloglog"
  if(x$margins[2] == "cauchit") m2l <- "cauchit"  
  if(x$margins[2] == "GEVlink") m2l <- "GEVlink" 
  
  if(x$margins[3] == "probit")  m3l <- "probit"
  if(x$margins[3] == "logit")   m3l <- "logit"
  if(x$margins[3] == "cloglog") m3l <- "cloglog"
  if(x$margins[3] == "cauchit") m3l <- "cauchit"  
  if(x$margins[3] == "GEVlink") m3l <- "GEVlink"   
  
 
 
list(cont1par = cont1par, cont2par = cont2par, cont3par = cont3par, cop1 = cop1, cop2 = cop2, lind1 = lind1, lind2 = lind2, m1l = m1l, m2l = m2l, m3l = m3l)
  

}

