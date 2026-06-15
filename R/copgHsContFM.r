copgHsContFM <- function(p1, p2, teta, teta.st, BivD, par2 = NULL, nu.st = NULL, T.BivD = NULL){
     
    
########################################################################################   
# Transformations on teta parameter
########################################################################################   

BivD2mixed <- c("C0C90","C0C270","C180C90","C180C270", "J0J90","J0J270","J180J90","J180J270", "G0G90","G0G270","G180G90","G180G270", "GAL0GAL90", "GAL0GAL270", "GAL180GAL90", "GAL180GAL270") 


#if(BivD %in% c("T") ) {
#
#dernu.dernu.st <- der2nu.dernu.stnu.st <- exp(nu.st)  not used for fixed margins
#
#}
 
 
if(BivD %in% c("N","AMH","FGM","T") ) {

derteta.derteta.st <- 1/cosh(teta.st)^2
der2teta.derteta.stteta.st <- -(2 * (sinh(teta.st) * cosh(teta.st))/(cosh(teta.st)^2)^2)
       
}  

if(BivD %in% c("F") ) {

derteta.derteta.st <- 1
der2teta.derteta.stteta.st <- 0
       
} 
   


if(BivD %in% c("HO") ) {

derteta.derteta.st         <- exp(-teta.st)/(1 + exp(-teta.st))^2 
der2teta.derteta.stteta.st <- -((1 - 2 * (exp(-teta.st)/(1 + exp(-teta.st)))) * exp(-teta.st)/(1 + exp(-teta.st))^2) 
       
} 
   
   
   

if(BivD %in% c("C0", "C180","GAL0","GAL180","J0", "J180","G0", "G180","PL") )   derteta.derteta.st <- der2teta.derteta.stteta.st <-  exp(teta.st) 
if(BivD %in% c("C90","C270","GAL90","GAL270","J90","J270","G90","G270") )       derteta.derteta.st <- der2teta.derteta.stteta.st <- -exp(teta.st)  



if(BivD %in% c("C0", "C180","GAL0","GAL180","J0", "J180","G0", "G180")  && T.BivD %in% BivD2mixed){ derteta.derteta.st         <- 1 # der teta.st wrt teta.st  
                                                                                                    der2teta.derteta.stteta.st <- 0} 

if(BivD %in% c("C90","C270","GAL90","GAL270","J90","J270","G90","G270") && T.BivD %in% BivD2mixed){ derteta.derteta.st         <- 1 # see comments in main file
                                                                                                    der2teta.derteta.stteta.st <- 0}


  

########################################################################################
########################################################################################


########################################################################################
# Rotations
########################################################################################

if(BivD %in% c("C90","J90","G90","GAL90") ) {
p1 <- 1 - p1 
teta <- -teta
}  

if(BivD %in% c("C180","J180","G180","GAL180") ) {
p1 <- 1 - p1
p2 <- 1 - p2
}  

if(BivD %in% c("C270","J270","G270","GAL270") ) {
p2 <- 1 - p2 
teta <- -teta 
}   
   
########################################################################################   
########################################################################################
     
         
if(BivD %in% c("GAL0","GAL90","GAL180","GAL270")){
     

der2h.derp1teta <- (((-log(p1))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p1))-(-log(p1))^(1+teta)*((1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2+(((((-log(p2))^(1+
teta)+teta*(-log(p2))^(1+teta)*log(-log(p2)))*(1+
1/teta)-(-log(p2))^(1+teta)/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)-(-log(p2))^(1+teta)*((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)*(1+1/teta))/(-log(p1))^(1+
teta)-((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*((1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))-1))/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2+((log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)+log(1/(-log(p1))^teta+1/(-log(p2))^teta)/(teta*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)))*(1+
teta*(-log(p2))^(1+teta)*(1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+teta)*
((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2)-((1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))/((-log(p1))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))+1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))))/teta-
teta*(-log(p2))^(1+teta)*((-log(p1))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2*log(-log(p1))+2*((-log(p1))^(1+teta)*
(-log(p2))^(1+teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p2))-(-log(p2))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)))*
(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2)^2)*exp((1/(-log(p1))^teta+
1/(-log(p2))^teta)^-(1/teta))


der2c.derrho.derrho <- ((((-log(p1))^(1+teta)*((((((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta+2*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))/teta+(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta))/teta+(((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)*(1+1/teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))/teta+
(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))^2/(-log(p1))^teta+log(-log(p2))^2/(-log(p2))^teta))+
((-log(p1))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p1))-2*((-log(p1))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2)))*
log(-log(p1)))*(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))+((-log(p1))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*log(-log(p1))-
(-log(p1))^(1+teta)*((1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2))*(((-log(p2))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*log(-log(p2))-
(-log(p2))^(1+teta)*((1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2))/((-log(p2))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2-
2*((-log(p1))^(1+teta)*((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p1))-(-log(p1))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)/((-log(p1))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2)))/((-log(p1))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2+
(((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*((1+
1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(((-log(p1))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p1))-(-log(p1))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))/((-log(p1))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2-2*((-log(p2))^(1+teta)*(((((-log(p2))^(1+
teta)+teta*(-log(p2))^(1+teta)*log(-log(p2)))*
(1+1/teta)-(-log(p2))^(1+teta)/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)-(-log(p2))^(1+teta)*
((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)*
(1+1/teta))/(-log(p1))^(1+teta)-((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))-1))*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2))+(((1+1/teta)*(2*(-log(p2))^(1+
teta)+teta*(-log(p2))^(1+teta)*log(-log(p2)))*
log(-log(p2))-((-log(p2))^(1+teta)*log(-log(p2))+
((-log(p2))^(1+teta)+teta*(-log(p2))^(1+teta)*
log(-log(p2))-(-log(p2))^(1+teta))/teta)/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)+((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta)*(2*((-log(p2))^(1+teta)/teta)-
((-log(p2))^(1+teta)+teta*(-log(p2))^(1+teta)*
log(-log(p2)))*(1+1/teta))/teta-((-log(p2))^(1+
teta)*((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta)*log(-log(p2))-(-log(p2))^(1+
teta)*((((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-1)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))/teta+
((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-1)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
2)*(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta))*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))^2/(-log(p1))^teta+
log(-log(p2))^2/(-log(p2))^teta)))*(1+1/teta))/(-log(p1))^(1+
teta)-((-log(p1))^(1+teta-2*(1+teta))*((((-log(p2))^(1+
teta)+teta*(-log(p2))^(1+teta)*log(-log(p2)))*
(1+1/teta)-(-log(p2))^(1+teta)/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)-(-log(p2))^(1+teta)*
((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-1)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)*
(1+1/teta))*log(-log(p1))+((-log(p2))^(1+teta)*
((((((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta+
2*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))/teta+(1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta))*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))/teta+
(((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)*
(1+1/teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))/teta+
(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))^2/(-log(p1))^teta+log(-log(p2))^2/(-log(p2))^teta))+
((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-2*((-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2)))*
log(-log(p2)))*(1/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))-1)))/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2+((((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta))-
1)*(1+teta*(-log(p2))^(1+teta)*(1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2)-
((1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))+1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))))/teta+2*(((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p1))-(-log(p1))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2))*(1-1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)))/((-log(p1))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2)+
2*((((((-log(p2))^(1+teta)+teta*(-log(p2))^(1+teta)*
log(-log(p2)))*(1+1/teta)-(-log(p2))^(1+teta)/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)-(-log(p2))^(1+
teta)*((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta)*(1+1/teta))/(-log(p1))^(1+
teta)-((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p2))-(-log(p2))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))-1))/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2)-2*(teta*(-log(p2))^(1+
teta)*((-log(p1))^(1+teta)*((-log(p2))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2*
log(-log(p1))+2*((-log(p1))^(1+teta)*(-log(p2))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*((1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)))*
(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2)^2))*((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)))-
(((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)-((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta))*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/(teta*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta))^2+
(1/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))-((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(2*(1+1/teta)))*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(log(-log(p1))^2/(-log(p1))^teta+
log(-log(p2))^2/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))*(1+teta*(-log(p2))^(1+
teta)*(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2)-((1-1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))/((-log(p1))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))+1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))))/teta-(((((-log(p1))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2*
log(-log(p1))+2*((-log(p1))^(1+teta)*(-log(p2))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*((1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)))*
((-log(p2))^(1+teta)+teta*(-log(p2))^(1+teta)*log(-log(p2)))+
teta*(-log(p2))^(1+teta)*(((-log(p1))^(1+teta)*
((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2*log(-log(p1))+2*((-log(p1))^(1+
teta)*(-log(p2))^(1+teta)*((-log(p2))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(-log(p2))-(-log(p2))^(1+teta)*((1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2))*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))*log(-log(p1))+2*(((-log(p1))^(1+teta)*
(-log(p2))^(1+teta)*((-log(p2))^(1+teta)*((((((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta+2*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))/teta+(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta))/teta+(((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta)*(1+1/teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta))/teta+(1+1/teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))^2/(-log(p1))^teta+
log(-log(p2))^2/(-log(p2))^teta))+((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-2*((-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2)))*
log(-log(p2)))+((-log(p1))^(1+teta)*(-log(p2))^(1+
teta)*log(-log(p1))+(-log(p1))^(1+teta)*(-log(p2))^(1+
teta)*log(-log(p2)))*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p2))-(-log(p2))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2)))*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)-
(-log(p1))^(1+teta)*(-log(p2))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))))*
(1+1/teta)-(-log(p2))^(1+teta)*((-log(p1))^(1+teta)*
((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2*log(-log(p1))+2*((-log(p1))^(1+teta)*
(-log(p2))^(1+teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(-log(p2))-(-log(p2))^(1+
teta)*((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2))*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)-
((-log(p1))^(1+teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2*log(-log(p1))+
2*((-log(p1))^(1+teta)*(-log(p2))^(1+teta)*
((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)))*
((-log(p2))^(1+teta)*((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)+
2*(teta*(-log(p1))^(1+teta)*(-log(p2))^(1+
teta)*((-log(p1))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2*log(-log(p1))+2*((-log(p1))^(1+
teta)*(-log(p2))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)))*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)/((-log(p1))^(1+
teta)*((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))^2)^2))*(1+
1/teta))/((-log(p1))^(1+teta)*((-log(p2))^(1+teta)*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2)^2)*
exp((1/(-log(p1))^teta+1/(-log(p2))^teta)^-(1/teta))




der2h.derteta.teta.st <- p1*(((-log(p2))^(1+teta)*((((((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta+2*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))/teta+(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta))/teta+(((1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta-1)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta)*(1+1/teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta))/teta+
(1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))^2/(-log(p1))^teta+log(-log(p2))^2/(-log(p2))^teta))+
((-log(p2))^(1+teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-2*((-log(p2))^(1+teta)*
((1+1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2)))*
log(-log(p2))-2*((-log(p2))^(1+teta)*((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(-log(p2))-(-log(p2))^(1+teta)*((1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))^2*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2))/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta))^2+
((((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta))-
1)*(1-1/((-log(p2))^(1+teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)))/teta+2*(((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta)*log(-log(p2))-(-log(p2))^(1+teta)*((1+
1/teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta^2))/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))^2))*((log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)+log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)))-
(((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)-
((1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta-
1)*(log(-log(p1))/(-log(p1))^teta+log(-log(p2))/(-log(p2))^teta)+
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta)*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/teta))*
log(1/(-log(p1))^teta+1/(-log(p2))^teta)/(teta*
(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1/teta))^2+
(1/(teta*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))-((1+1/teta)*(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1/teta)*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta)*log(1/(-log(p1))^teta+
1/(-log(p2))^teta)/teta^2)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(2*(1+1/teta)))*(log(-log(p1))/(-log(p1))^teta+
log(-log(p2))/(-log(p2))^teta)+(log(-log(p1))^2/(-log(p1))^teta+
log(-log(p2))^2/(-log(p2))^teta)/(1/(-log(p1))^teta+
1/(-log(p2))^teta)^(1+1/teta))*(1-1/((-log(p2))^(1+
teta)*(1/(-log(p1))^teta+1/(-log(p2))^teta)^(1+
1/teta))))/teta)*exp((1/(-log(p1))^teta+1/(-log(p2))^teta)^-(1/teta))



   
}     
     
     

if(BivD %in% c("HO")){


der2h.derteta.teta.st <- ((-log(p2))^(1/teta - 1) * ((((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - (((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2)) - (((-log(p1))^(1/teta) * 
    log(-log(p1))/teta + 2 * (-log(p1))^(1/teta)) * log(-log(p1)) + 
    ((-log(p2))^(1/teta) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta)) * 
        log(-log(p2))) * (teta - 1)/teta) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) + ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) + (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 3) * (teta - 2)/teta^2) * (teta - 
    1)))/teta^2) + (2 * (((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 
    1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta)) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2)) - ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((2 * (-log(p2))^(1/teta - 
    1) - ((-log(p2))^(1/teta - 1) * log(-log(p2))/teta + 2 * 
    (-log(p2))^(1/teta - 1) + 2 * (-log(p2))^(1/teta - 1))) * 
    log(-log(p2))/teta + 2 * (-log(p2))^(1/teta - 1) - 2 * (-log(p2))^(1/teta - 
    1))/teta)/teta - ((-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1)/teta + (1 - ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta))) * 
    (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
        (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 1)/teta) - (((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) + (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) - ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
    log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * (teta - 1)/teta)/teta) + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * (((-log(p1))^(1/teta) - ((-log(p1))^(1/teta) * log(-log(p1))/teta + 
    2 * (-log(p1))^(1/teta))) * log(-log(p1)) + ((-log(p2))^(1/teta) - 
    ((-log(p2))^(1/teta) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta))) * 
    log(-log(p2)))/teta)/teta) + 2 * (((-log(p2))^(1/teta - 1) * 
    (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 1) * 
        log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
        log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * 
        ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 2) * 
        (teta - 1)/teta^2) + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * ((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 1) + 
    (-log(p2))^(1/teta - 1) * log(-log(p2))/teta))/teta) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
    ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
        log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1)/teta)))) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/p2 



der2h.derp1teta <- ((-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * ((-log(p2))^(1/teta - 1) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) + ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((-log(p2))^(1/teta - 1) - 
    ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta))/teta) + 
    (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta - 1) * (((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
        (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) + 
        ((-log(p1))^(1/teta - 1) - ((-log(p1))^(1/teta - 1) + 
            (-log(p1))^(1/teta - 1) * log(-log(p1))/teta)) * 
            ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
                1)/teta - (-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 1) * (((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
        (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 1)/teta)) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 1) - ((-log(p2))^(1/teta - 
    1) * ((-log(p1))^(1/teta - 1) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 3) * (teta - 2)/teta^2) * (teta - 
    1) + ((-log(p1))^(1/teta - 1) - ((-log(p1))^(1/teta - 1) + 
    (-log(p1))^(1/teta - 1) * log(-log(p1))/teta) * (teta - 1)/teta) * 
    ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 2)) + 
    ((-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        2) * ((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 
        1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta))/teta - 
        (-log(p1))^(1/teta - 1) * (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 2) * (((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
            (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 1)/teta)) * (teta - 
        1))/teta) * exp(-((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta)/(p1 * 
    p2) 



der2c.derrho.derrho <- ((-log(p1))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * ((-log(p2))^(1/teta - 1) * ((((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - (((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2)) - (((-log(p1))^(1/teta) * 
    log(-log(p1))/teta + 2 * (-log(p1))^(1/teta)) * log(-log(p1)) + 
    ((-log(p2))^(1/teta) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta)) * 
        log(-log(p2))) * (teta - 1)/teta) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) + ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) + (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 3) * (teta - 2)/teta^2) * (teta - 
    1)))/teta^2) + (2 * (((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 
    1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta)) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2)) - ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((2 * (-log(p2))^(1/teta - 
    1) - ((-log(p2))^(1/teta - 1) * log(-log(p2))/teta + 2 * 
    (-log(p2))^(1/teta - 1) + 2 * (-log(p2))^(1/teta - 1))) * 
    log(-log(p2))/teta + 2 * (-log(p2))^(1/teta - 1) - 2 * (-log(p2))^(1/teta - 
    1))/teta)/teta - ((-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1)/teta + (1 - ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta))) * 
    (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
        (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 1)/teta) - (((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) + (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) - ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
    log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    2) * (teta - 1)/teta)/teta) + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * (((-log(p1))^(1/teta) - ((-log(p1))^(1/teta) * log(-log(p1))/teta + 
    2 * (-log(p1))^(1/teta))) * log(-log(p1)) + ((-log(p2))^(1/teta) - 
    ((-log(p2))^(1/teta) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta))) * 
    log(-log(p2)))/teta)/teta) + 2 * (((-log(p2))^(1/teta - 1) * 
    (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 1) * 
        log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
        log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * 
        ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 2) * 
        (teta - 1)/teta^2) + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * ((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 1) + 
    (-log(p2))^(1/teta - 1) * log(-log(p2))/teta))/teta) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
    ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
        log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1)/teta)))) + (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta - 
    1) * ((((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - (((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2)) - (((-log(p1))^(1/teta) * 
    log(-log(p1))/teta + 2 * (-log(p1))^(1/teta)) * log(-log(p1)) + 
    ((-log(p2))^(1/teta) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta)) * 
        log(-log(p2))) * (teta - 1)/teta) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) + ((-log(p1))^(1/teta) * 
    log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) + (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 3) * (teta - 2)/teta^2) * (teta - 
    1)))/teta^2) + (2 * (((-log(p1))^(1/teta - 1) - ((-log(p1))^(1/teta - 
    1) + (-log(p1))^(1/teta - 1) * log(-log(p1))/teta)) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2)) - ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((2 * (-log(p1))^(1/teta - 
    1) - ((-log(p1))^(1/teta - 1) * log(-log(p1))/teta + 2 * 
    (-log(p1))^(1/teta - 1) + 2 * (-log(p1))^(1/teta - 1))) * 
    log(-log(p1))/teta + 2 * (-log(p1))^(1/teta - 1) - 2 * (-log(p1))^(1/teta - 
    1))/teta)/teta) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1) + 2 * (((-log(p1))^(1/teta - 1) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) + ((-log(p1))^(1/teta - 
    1) - ((-log(p1))^(1/teta - 1) + (-log(p1))^(1/teta - 1) * 
    log(-log(p1))/teta)) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
    1)/teta) * ((-log(p2))^(1/teta - 1) * (((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
    (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta^2) + ((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^(teta - 1) * ((-log(p2))^(1/teta - 1) - 
    ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta))/teta - 
    (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        1) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * 
        log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
        log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * 
        ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 1)/teta))) - 
    ((-log(p1))^(1/teta - 1) * ((-log(p2))^(1/teta - 1) * ((((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 2) * log((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
        (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 3) * (teta - 2)/teta^2) * 
        log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - (((-log(p1))^(1/teta) * 
        log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2)) - 
        (((-log(p1))^(1/teta) * log(-log(p1))/teta + 2 * (-log(p1))^(1/teta)) * 
            log(-log(p1)) + ((-log(p2))^(1/teta) * log(-log(p2))/teta + 
            2 * (-log(p2))^(1/teta)) * log(-log(p2))) * (teta - 
            2)/teta) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        3) + ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
        log(-log(p2))) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        3) + (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        3) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
        ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
            log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
            4) * (teta - 3)/teta^2) * (teta - 2)))/teta^2) * 
        (teta - 1) + 2 * (((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 
        1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta) * 
        (teta - 1)/teta) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
        2) * log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
        ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
            log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
            3) * (teta - 2)/teta^2)) - ((-log(p1))^(1/teta) + 
        (-log(p2))^(1/teta))^(teta - 2) * ((2 * ((-log(p2))^(1/teta - 
        1) * log(-log(p2))) - (((-log(p2))^(1/teta - 1) * log(-log(p2))/teta + 
        2 * (-log(p2))^(1/teta - 1) + 2 * (-log(p2))^(1/teta - 
        1)) * log(-log(p2))/teta + 2 * (-log(p2))^(1/teta - 1)) * 
        (teta - 1))/teta + 2 * (-log(p2))^(1/teta - 1))/teta) + 
        (2 * (((-log(p1))^(1/teta - 1) - ((-log(p1))^(1/teta - 
            1) + (-log(p1))^(1/teta - 1) * log(-log(p1))/teta)) * 
            ((-log(p2))^(1/teta - 1) * (((-log(p1))^(1/teta) + 
                (-log(p2))^(1/teta))^(teta - 2) * log((-log(p1))^(1/teta) + 
                (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * 
                log(-log(p1)) + (-log(p2))^(1/teta) * log(-log(p2))) * 
                ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
                  3) * (teta - 2)/teta^2) * (teta - 1) + ((-log(p1))^(1/teta) + 
                (-log(p2))^(1/teta))^(teta - 2) * ((-log(p2))^(1/teta - 
                1) - ((-log(p2))^(1/teta - 1) + (-log(p2))^(1/teta - 
                1) * log(-log(p2))/teta) * (teta - 1)/teta))) - 
            (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta) + 
                (-log(p2))^(1/teta))^(teta - 2) * ((2 * (-log(p1))^(1/teta - 
                1) - ((-log(p1))^(1/teta - 1) * log(-log(p1))/teta + 
                2 * (-log(p1))^(1/teta - 1) + 2 * (-log(p1))^(1/teta - 
                1))) * log(-log(p1))/teta + 2 * (-log(p1))^(1/teta - 
                1) - 2 * (-log(p1))^(1/teta - 1)) * (teta - 1)/teta)/teta - 
        ((-log(p1))^(1/teta - 1) * (-log(p2))^(1/teta - 1) * 
            ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
                2) * ((((-log(p1))^(1/teta) * log(-log(p1)) + 
            (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 1)/teta + (1 - ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^teta) * log((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))) * (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * 
            log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
            ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
                log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
                1)/teta) - (((-log(p1))^(1/teta) * log(-log(p1)) + 
            (-log(p2))^(1/teta) * log(-log(p2))) * (((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 1) * log((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta)) + (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
            1) - ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
            log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
            2) * (teta - 1)/teta)/teta) + ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 1) * (((-log(p1))^(1/teta) - 
            ((-log(p1))^(1/teta) * log(-log(p1))/teta + 2 * (-log(p1))^(1/teta))) * 
            log(-log(p1)) + ((-log(p2))^(1/teta) - ((-log(p2))^(1/teta) * 
            log(-log(p2))/teta + 2 * (-log(p2))^(1/teta))) * 
            log(-log(p2)))/teta)/teta) * (teta - 1) + 2 * (((-log(p1))^(1/teta - 
            1) * ((-log(p2))^(1/teta - 1) * (((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 2) * log((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta)) - ((-log(p1))^(1/teta) * log(-log(p1)) + 
            (-log(p2))^(1/teta) * log(-log(p2))) * ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 3) * (teta - 2)/teta^2) * 
            (teta - 1) + ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
            2) * ((-log(p2))^(1/teta - 1) - ((-log(p2))^(1/teta - 
            1) + (-log(p2))^(1/teta - 1) * log(-log(p2))/teta) * 
            (teta - 1)/teta)) + (-log(p2))^(1/teta - 1) * ((-log(p1))^(1/teta - 
            1) - ((-log(p1))^(1/teta - 1) + (-log(p1))^(1/teta - 
            1) * log(-log(p1))/teta)) * ((-log(p1))^(1/teta) + 
            (-log(p2))^(1/teta))^(teta - 2) * (teta - 1)/teta) * 
            (((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^teta * 
                log((-log(p1))^(1/teta) + (-log(p2))^(1/teta)) - 
                ((-log(p1))^(1/teta) * log(-log(p1)) + (-log(p2))^(1/teta) * 
                  log(-log(p2))) * ((-log(p1))^(1/teta) + (-log(p2))^(1/teta))^(teta - 
                  1)/teta))))/teta) * exp(-((-log(p1))^(1/teta) + 
    (-log(p2))^(1/teta))^teta)/(p1 * p2) 


}



if(BivD %in% c("PL")){



der2h.derteta.teta.st <- -(((4 * (p1 + p2) - 8 * p1)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1))) - ((2 * (((p1 + p2) * 
    (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 1) + 
    4 * teta)) * (2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - p1 * 
    (4 * (teta - 1) + 4 * teta)) + 0.5 * ((2 * ((p1 + p2) * (teta - 
    1) + 1) - 4 * (p1 * teta)) * (2 * (p1 + p2)^2 - (1.5 * ((2 * 
    (((p1 + p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * 
    (teta - 1) + 4 * teta))^2/(((p1 + p2) * (teta - 1) + 1)^2 - 
    4 * (p1 * p2 * teta * (teta - 1)))) + 8 * (p1 * p2))) * (teta - 
    1)))/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1)))^1.5)/(4 * (teta - 1)) + (4 * (1 - 0.5 * ((2 * 
    (1 + 2 * ((p1 + p2) * (teta - 1))) - p1 * (4 * (teta - 1) + 
    4 * teta))/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * 
    p2 * teta * (teta - 1))) - 0.5 * ((2 * (((p1 + p2) * (teta - 
    1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta)) * 
    (2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (teta - 
    1)/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1)))^1.5))) - 16 * ((teta - (0.5 * ((2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p1 * teta)) * (teta - 1)/sqrt(((p1 + 
    p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))) + 
    1)) * (teta - 1)/(2 * (teta - 1))^2))/(2 * (teta - 1))^2) 



der2h.derp1teta <- ((0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * 
    (2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - (1.5 * ((2 * (((p1 + 
        p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 
        1) + 4 * teta)) * (2 * ((p1 + p2) * (teta - 1) + 1) - 
        4 * (p2 * teta)) * (teta - 1)/(((p1 + p2) * (teta - 1) + 
        1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))) + p2 * (4 * 
        (teta - 1) + 4 * teta)))) + 0.5 * ((2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p2 * teta)) * (2 * (1 + 2 * ((p1 + 
    p2) * (teta - 1))) - p1 * (4 * (teta - 1) + 4 * teta)))) * 
    (teta - 1)/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * 
    teta * (teta - 1)))^1.5 + (0.5 * ((2 * (((p1 + p2) * (teta - 
    1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta)) * 
    (2 * (teta - 1) - 4 * teta) * (teta - 1)/(((p1 + p2) * (teta - 
    1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5) + 4 * 
    (teta/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * 
        teta * (teta - 1))))))/(4 * (teta - 1)) + ((2 * (teta - 
    1) - 4 * teta)/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * 
    (p1 * p2 * teta * (teta - 1))) - 0.5 * ((2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p1 * teta)) * (2 * ((p1 + p2) * (teta - 
    1) + 1) - 4 * (p2 * teta)) * (teta - 1)/(((p1 + p2) * (teta - 
    1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5)) * (teta - 
    1)/(2 * (teta - 1))^2  


der2c.derrho.derrho <- -((((0.5 * ((2 * (((p1 + p2) * (teta - 1) + 1) * (p1 + p2)) - 
    p1 * p2 * (4 * (teta - 1) + 4 * teta)) * (2 * (teta - 1) - 
    4 * teta)) + 0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * 
    (p1 * teta)) * (2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - 
    p2 * (4 * (teta - 1) + 4 * teta))) + 0.5 * ((2 * ((p1 + p2) * 
    (teta - 1) + 1) - 4 * (p2 * teta)) * (2 * (1 + 2 * ((p1 + 
    p2) * (teta - 1))) - (1.5 * ((2 * (((p1 + p2) * (teta - 1) + 
    1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta)) * 
    (2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (teta - 
    1)/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1)))) + p1 * (4 * (teta - 1) + 4 * teta)))))/(((p1 + 
    p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5 + 
    8 * (((2 * (teta - 1) - 4 * teta)/sqrt(((p1 + p2) * (teta - 
        1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1))) - 0.5 * 
        ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * 
            (2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p2 * teta)) * 
            (teta - 1)/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * 
            (p1 * p2 * teta * (teta - 1)))^1.5)) * (teta - 1)/(2 * 
        (teta - 1))^2)) * (teta - 1) + (0.5 * ((2 * (((p1 + p2) * 
    (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 1) + 
    4 * teta)) * (2 * (teta - 1) - 4 * teta)) + 0.5 * ((2 * ((p1 + 
    p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (2 * (1 + 2 * 
    ((p1 + p2) * (teta - 1))) - p2 * (4 * (teta - 1) + 4 * teta))) + 
    0.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p2 * teta)) * 
        (2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - (1.5 * ((2 * 
            (((p1 + p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * 
            p2 * (4 * (teta - 1) + 4 * teta)) * (2 * ((p1 + p2) * 
            (teta - 1) + 1) - 4 * (p1 * teta)) * (teta - 1)/(((p1 + 
            p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 
            1)))) + p1 * (4 * (teta - 1) + 4 * teta))))) * (teta - 
    1)/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * 
    (teta - 1)))^1.5 + 8 * teta/sqrt(((p1 + p2) * (teta - 1) + 
    1)^2 - 4 * (p1 * p2 * teta * (teta - 1))))/(2 * (teta - 1))^2 + 
    ((4 * (teta * (2 * (((p1 + p2) * (teta - 1) + 1) * (p1 + 
        p2)) - p1 * p2 * (4 * (teta - 1) + 4 * teta))) - ((0.5 * 
        ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * 
            (4 * (p1 + p2) - 8 * p2)) + 0.5 * ((2 * ((p1 + p2) * 
        (teta - 1) + 1) - 4 * (p2 * teta)) * (4 * (p1 + p2) - 
        ((1.5 * ((2 * ((p1 + p2) * (teta - 1) + 1) - 4 * (p1 * 
            teta)) * (2 * (p1 + p2)^2 - (2.5 * ((2 * (((p1 + 
            p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * 
            (teta - 1) + 4 * teta))^2/(((p1 + p2) * (teta - 1) + 
            1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))) + 8 * 
            (p1 * p2))) * (teta - 1)) + 3 * ((2 * (((p1 + p2) * 
            (teta - 1) + 1) * (p1 + p2)) - p1 * p2 * (4 * (teta - 
            1) + 4 * teta)) * (2 * (1 + 2 * ((p1 + p2) * (teta - 
            1))) - p1 * (4 * (teta - 1) + 4 * teta))))/(((p1 + 
            p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 
            1))) + 8 * p1))) + 0.5 * ((2 * (p1 + p2)^2 - (1.5 * 
        ((2 * (((p1 + p2) * (teta - 1) + 1) * (p1 + p2)) - p1 * 
            p2 * (4 * (teta - 1) + 4 * teta))^2/(((p1 + p2) * 
            (teta - 1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 
            1)))) + 8 * (p1 * p2))) * (2 * (teta - 1) - 4 * teta))) * 
        (teta - 1) + (2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - 
        (1.5 * ((2 * (((p1 + p2) * (teta - 1) + 1) * (p1 + p2)) - 
            p1 * p2 * (4 * (teta - 1) + 4 * teta)) * (2 * ((p1 + 
            p2) * (teta - 1) + 1) - 4 * (p1 * teta)) * (teta - 
            1)/(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * 
            teta * (teta - 1)))) + p1 * (4 * (teta - 1) + 4 * 
            teta))) * (2 * (1 + 2 * ((p1 + p2) * (teta - 1))) - 
        p2 * (4 * (teta - 1) + 4 * teta))))/(((p1 + p2) * (teta - 
        1) + 1)^2 - 4 * (p1 * p2 * teta * (teta - 1)))^1.5 - 
        4/sqrt(((p1 + p2) * (teta - 1) + 1)^2 - 4 * (p1 * p2 * 
            teta * (teta - 1))))/(4 * (teta - 1))) 


}



     
     
     
    
     
     
     
if(BivD == "N"){



der2h.derteta.teta.st <- ((qnorm(p1) - teta * qnorm(p2)) * (qnorm(p2) - teta * 
    (qnorm(p1) - teta * qnorm(p2))/(1 - teta^2)) * (teta * (qnorm(p1) - 
    teta * qnorm(p2))/sqrt(1 - teta^2) - qnorm(p2) * sqrt(1 - 
    teta^2)) + (qnorm(p1) + teta * (2 * (teta * qnorm(p1)) - 
    3 * qnorm(p2)))/sqrt(1 - teta^2)) * dnorm((qnorm(p1) - teta * 
    qnorm(p2))/sqrt(1 - teta^2))/(1 - teta^2)^2
    
    

der2h.derp1teta <- exp(-(teta * (teta * (qnorm(p1)^2 + qnorm(p2)^2) - 
    2 * (qnorm(p1) * qnorm(p2)))/(2 * (1 - teta^2)))) * (qnorm(p1) * 
    qnorm(p2) + teta * (1 - (qnorm(p1)^2 + qnorm(p2)^2 + teta * 
    (teta * (qnorm(p1)^2 + qnorm(p2)^2) - 2 * (qnorm(p1) * qnorm(p2)))/(1 - 
    teta^2))))/((1 - teta^2) * sqrt(1 - teta^2))
                                                            

#c <- 1/sqrt(1 - teta^2)*exp(  - (teta^2*( qnorm(p1)^2 +  qnorm(p2)^2 ) - 2*teta*qnorm(p1)*qnorm(p2) ) / (2*(1 - teta^2)) ) 

			t6 = qnorm(p1) 
			t7 = qnorm(p2)
			t1 = t6*t7
			t2 = teta*teta
			t3 = 1-t2
			t4 = 4*t3*t3
			t5 = 1/t4
			t12 = t6*t6
			t13 = t7*t7
			t14 = 2*teta*t6*t7-t12-t13
			t21 = t14*t5
			t26 = 1/t3/2
			t29 = exp(t12/2+t13/2+t14*t26)
			t31 = sqrt(t3)
			t32 = 1/t31
			t38 = 2*t1*t26+4*t21*teta
			t39 = t38*t38
			t44 = 1/t31/t3
			t48 = t3*t3
der2c.derrho.derrho <- (16.0*t1*t5*teta+16.0*t14/t4/t3*t2+4.0*t21)*t29*t32+t39*t29*t32+2.0*t38*t29*t44*teta+3.0*t29/t31/t48*t2+t29*t44



}





if(BivD == "T"){

der2h.derteta.teta.st <- BiCopHfuncDeriv2(p1, p2, family = 2, teta, par2, deriv = "par")


#if(Cont == TRUE){
#
#eps <- 0.001 # because of num deriv below
#m.p <- 0.999
#
#p1 <- ifelse(p1 > m.p, m.p, p1) 
#p1 <- ifelse(p1 < eps, eps, p1) 
#
#p2 <- ifelse(p2 > m.p, m.p, p2) 
#p2 <- ifelse(p2 < eps, eps, p2) 
#
#}


der2h.derp1teta     <-  BiCopDeriv(p1, p2, family = 2, teta, par2, deriv = "par")  
der2c.derrho.derrho <- BiCopDeriv2(p1, p2, family = 2, teta, par2, deriv = "par")

}







if(BivD == "FGM"){

der2h.derteta.teta.st <- der3C.derp1tetateta <- 0


der2h.derp1teta <-  (1 - 2 * p1) * (1 - 2 * p2)

der2c.derrho.derrho <- 0


}











     


if(BivD == "AMH"){

 

der2h.derteta.teta.st <- -(p1 * (1 - p1)^2 * (1 - p2) * (2 * (p2 * (1 + teta * (1 - p1) * 
    (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) - 2 * (1 - p2 * 
    (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
        (1 - p2))))))/(1 - teta * (1 - p1) * (1 - p2))^3)



der2h.derp1teta <-  -((((1 - p2) * (2 * (teta * (p1 * (1 - p2 * (2 + 2 * (teta * 
    (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + 
    p2 * (1 - p1))/(1 - teta * (1 - p1) * (1 - p2))) - 1) + p2 * 
    (1 - 2 * (p1 * teta * (1 - p2) * (1 + teta * (1 - p1) * (1 - 
        p2)/(1 - teta * (1 - p1) * (1 - p2)))/(1 - teta * (1 - 
        p1) * (1 - p2))))) * (1 - p1) + p1 * (1 - p2 * (2 + 2 * 
    (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))))/(1 - 
    teta * (1 - p1) * (1 - p2))^2)



der2c.derrho.derrho <- -((((1 - p2) * (teta * (p1 * (1 - p2 * (2 + 2 * (teta * (1 - 
    p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2))))) + p2 * 
    (1 - p1))/(1 - teta * (1 - p1) * (1 - p2)) - 2 * (1 - teta * 
    (p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - 
    teta * (1 - p1) * (1 - p2)))) + 3 * (teta * (1 - p2) * (p1 * 
    (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - teta * 
        (1 - p1) * (1 - p2))))) + p2 * (1 - p1))/(1 - teta * 
    (1 - p1) * (1 - p2))) + 4 * p2 * (1 - 2 * (p1 * teta * (1 - 
    p2) * (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * 
    (1 - p2)))/(1 - teta * (1 - p1) * (1 - p2))))) * (1 - p1) + 
    p1 * (1 - p2 * (2 + 2 * (teta * (1 - p1) * (1 - p2)/(1 - 
        teta * (1 - p1) * (1 - p2))))) + p1 * (3 - p2 * (4 * 
    (1 + teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - 
        p2))) + 6 + teta * (1 - p1) * (1 - p2) * (2 * (2 + 2 * 
    (teta * (1 - p1) * (1 - p2)/(1 - teta * (1 - p1) * (1 - p2)))) + 
    6)/(1 - teta * (1 - p1) * (1 - p2))))) * (1 - p1) * (1 - 
    p2)/(1 - teta * (1 - p1) * (1 - p2))^3)



}




        


if(BivD %in% c("C0","C90","C180","C270")){


t1 = p1*p2;
t2 = -teta-1.0;
t3 = t1^(1.0*t2);
t4 = p1^(-1.0*teta);
t5 = p2^(-1.0*teta);
t6 = t4+t5-1.0;
t7 = -2.0-1/teta;
t8 = t6^(1.0*t7);
t9 = -t2*t3;
t10 = log(t1);
t11 = teta^2;
t12 = log(t6);
t13 = log(p1);
t14 = log(p2);

der2h.derp1teta <-  t3*t8-t9*t10*t8+t9*t8*(1/t11*t12+t7*(-t4*t13-t5*t14)/t6);



			t1 = p1*p2
			t2 = -teta-1
			t3 = t1^t2
			t4 = log(t1)
			t6 = p1^-teta
			t7 = p2^-teta
			t8 = t6+t7-1
			t10 = -2-1/teta
			t11 = t8^t10
			t15 = teta^2
			t16 = 1/t15
			t17 = log(t8)
			t19 = log(p1)
			t21 = log(p2)
			t24 = -t6*t19-t7*t21
			t26 = 1/t8
			t27 = t16*t17+t10*t24*t26
			t30 = -t2*t3
			t32 = t4*t4
			t14 = t27*t27
			t13 = t19*t19
			t12 = t21*t21
			t9 = t24*t24
			t5 = t8*t8
der2c.derrho.derrho <- -2*t3*t4*t11+2*t3*t11*t27+t30*t32*t11-2*t30*t4*t11*t27+t30*t11*t14+t30*t11*(-2/t15/teta*t17+2*t16*t24*t26+t10*(t6*t13+t7*t12)*t26-t10*t9/t5)


t2 = p2^(-1.0*teta-1.0);
t3 = log(p2);
t4 = t3*t3;
t6 = p1^(-1.0*teta);
t7 = p2^(-1.0*teta);
t8 = t6+t7-1.0;
t10 = -1.0-1/teta;
t11 = t8^(1.0*t10);
t14 = teta^2;
t15 = 1/t14;
t16 = log(t8);
t18 = log(p1);
t21 = -t6*t18-t7*t3;
t23 = 1/t8;
t25 = t15*t16+t10*t21*t23;
t29 = t2*t11;
t32 = t25*t25;
t12 = t18*t18;
t9 = t21*t21;
t5 = t8*t8;

der2h.derteta.teta.st <- t2*t4*t11-2.0*t2*t3*t11*t25+t29*t32+t29*(-2.0/t14/teta*t16+2.0*t15*t21*t23+t10*(t6*t12+t7*t4)*t23-t10*t9/t5)



}






if(BivD == "F"){



			t2 = exp(teta);
			t3 = t2-1.0;
			t4 = teta*p2;
			t5 = teta*p1;
			t7 = exp(t4+t5+teta);
			t10 = exp(t4+t5);
			t12 = exp(t4+teta);
			t14 = exp(t5+teta);
			t15 = t10-t12-t14+t2;
			t16 = t15*t15;
			t17 = 1/t16;
			t21 = teta*t3;
der2h.derp1teta <- t3*t7*t17+teta*t2*t7*t17+t21*(p2+p1+1.0)*t7*t17-2.0*t21*t7/t15/t16*((p2+p1)*t10-(p2+1.0)*t12-(p1+1.0)*t14+t2);
		



			t1 = exp(teta)
			t2 = teta*p2
			t3 = teta*p1
			t5 = exp(t2+t3+teta)
			t8 = exp(t2+t3)
			t10 = exp(t2+teta)
			t12 = exp(t3+teta)
			t13 = t8-t10-t12+t1
			t14 = t13*t13
			t15 = 1/t14
			t18 = t1-1.0
			t19 = p2+p1+1.0
			t21 = t5*t15
			t26 = 1/t14/t13
			t27 = p2 + p1
			t29 = p2+1.0
			t31 = p1+1.0
			t33 = t27*t8-t29*t10-t31*t12+t1
			t37 = teta*t1
			t43 = t5*t26
			t44 = t43*t33
			t47 = teta*t18
			t48 = t19*t19
			t11 = t14*t14
			t9 = t33*t33
			t7 = t27*t27
			t6 = t29*t29
			t4 = t31*t31
der2c.derrho.derrho <-	2.0*t1*t5*t15+2.0*t18*t19*t21-4.0*t18*t5*t26*t33+t37*t21+2.0*t37*t19*t5*t15-4.0*t37*t44+t47*t48*t5*t15-4.0*t47*t19*t44+6.0*t47*t5/t11*t9-2.0*t47*t43*(t7*t8-t6*t10-t4*t12+t1)




t1 = exp(teta);
t2 = teta*p1;
t3 = exp(t2);
t5 = t1*(t3-1.0);
t6 = teta*p2;
t8 = exp(t6+t2);
t10 = exp(t6+teta);
t12 = exp(t2+teta);
t13 = t8-t10-t12+t1;
t14 = 1/t13;
t16 = t1*p1;
t18 = t3*t14;
t20 = t13*t13;
t21 = 1/t20;
t23 = p2+p1;
t25 = p2+1.0;
t26 = p1+1.0;
t28 = t23*t8-t25*t10-t26*t12+t1;
t32 = p1*p1;
t42 = t28*t28;
t44 = t23*t23;
t47 = t25*t25;
t49 = t26*t26;

der2h.derteta.teta.st <- -t5*t14-2.0*t16*t18+2.0*t5*t21*t28-t1*t32*t18+2.0*t16*t3*t21*t28-2.0*t5/t20/t13*t42+t5*t21*(t44*t8-t47*t10-t49*t12+t1);





}






if(BivD %in% c("G0","G90","G180","G270")){

  
  
  			  t3 = log(p1);
  			  t4 = (-t3)^teta
  			  t5 = log(p2);
  			  t6 = (-t5)^teta
  			  t7 = t4+t6;
  			  t8 = 1/teta;
  			  t9 = t7^t8
  			  t10 = teta^2;
  			  t12 = log(t7);
  			  t13 = 1/t10*t12;
  			  t14 = log(-t3);
  			  t16 = log(-t5);
  			  t18 = t4*t14+t6*t16;
  			  t20 = 1/t7;
  			  t22 = -t13+t8*t18*t20;
  			  t24 = exp(-t9);
  			  t26 = t24/p1;
  			  t28 = 1/p2;
  			  t29 = -1.0+t8;
  			  t30 = t7^(2*t29)
  			  t32 = t3*t5;
  			  t33 = teta-1.0;
  			  t34 = t32^t33
  			  t35 = t7^(-1.0*t8)
  			  t36 = t33*t35;
  			  t17 = 1.0+t36;
  			  t15 = t34*t17;
  			  t11 = t26*t28;
  			  t2 = t30*t34;
  			  t1 = log(t32);
  			  
 der2h.derp1teta <- -t9*t22*t26*t28*t30*t15+t11*t30*(-2.0*t13+2.0*t29*t18*t20)*t15+t11*t2*t1*t17+t11*t2*(t35-t36*t22);
  




                        t3 = log(p1);
			t4 = (-t3)^teta
			t5 = log(p2)
			t6 = (-t5)^teta
			t7 = t4+t6
			t8 = 1/teta
			t9 = t7^t8
			t10 = teta*teta
			t11 = 1/t10
			t12 = log(t7)
			t13 = t11*t12
			t14 = log(-t3)
			t16 = log(-t5)
			t18 = t4*t14+t6*t16
			t20 = 1/t7
			t22 = -t13+t8*t18*t20
			t23 = t22*t22
			t25 = exp(-t9)
			t27 = t25/p1
			t29 = 1/p2
			t30 = -1.0+t8
			t31 = t7^(2*t30)
			t32 = t29*t31
			t33 = t3*t5
			t34 = teta-1.0
			t35 = t33^t34
			t36 = t7^-t8
			t37 = t34*t36
			t38 = 1.0+t37
			t39 = t35*t38
			t40 = t32*t39
			t44 = 1/t10/teta*t12
			t47 = t11*t18*t20
			t49 = t14*t14
			t51 = t16*t16
			t53 = t4*t49+t6*t51
			t56 = t18*t18
			t58 = t7*t7
			t59 = 1/t58
			t61 = 2.0*t44-2.0*t47+t8*t53*t20-t8*t56*t59
			t65 = t9*t9
			t70 = t9*t22*t27
			t74 = -2.0*t13+2.0*t30*t18*t20
			t75 = t74*t35
			t80 = log(t33)
			t87 = t36-t37*t22
			t88 = t35*t87
			t17 = t27*t29
			t15 = t74*t74
			t2 = t31*t35
			t1 = t80*t80
der2c.derrho.derrho <- -t9*t23*t27*t40-t9*t61*t27*t40+t65*t23*t27*t40-2.0*t70*t32*t75*t38-2.0*t70*t32*t35*t80*t38-2.0*t70*t32*t88+t17*t31*t15*t39+t17*t31*(4.0*t44-4.0*t47+2.0*t30*t53*t20-2.0*t30*t56*t59)*t39+2.0*t27*t32*t75*t80*t38+2.0*t17*t31*t74*t88+t17*t2*t1*t38+2.0*t17*t2*t80*t87+t17*t2*(-2.0*t36*t22+t37*t23-t37*t61)



  
			  t1 = log(p2);
			  t2 = (-t1)^teta
			  t3 = log(p1);
			  t4 = (-t3)^teta
			  t5 = t2+t4;
			  t6 = 1/teta;
			  t7 = t5^t6
			  t8 = teta^2
			  t9 = 1/t8;
			  t10 = log(t5);
			  t11 = t9*t10;
			  t12 = log(-t1);
			  t13 = t2*t12;
			  t14 = log(-t3);
			  t16 = t13+t4*t14;
			  t18 = 1/t5;
			  t20 = -t11+t6*t16*t18;
			  t21 = t20^2
			  t23 = exp(-t7);
			  t25 = t6-1.0;
			  t26 = t5^t25
			  t28 = 1/p2
			  t29 = 1/t1;
			  t30 = t28*t29;
			  t31 = t26*t2*t30;
			  t36 = 2.0/t8/teta*t10;
			  t39 = 2.0*t9*t16*t18;
			  t40 = t12^2;
			  t42 = t14^2;
			  t44 = t2*t40+t4*t42;
			  t47 = t16^2
			  t49 = t5^2;
			  t50 = 1/t49;
			  t56 = t7*t7;
			  t61 = t23*t26;
			  t62 = t7*t20*t61;
			  t65 = -t11+t25*t16*t18;
			  t70 = t13*t30;
			  t73 = t65*t65;
			  t15 = t2*t28*t29;
			  
der2h.derteta.teta.st <- t7*t21*t23*t31+t7*(t36-t39+t6*t44*t18-t6*t47*t50)*t23*t31-t56*t21*t23*t31+2.0*t62*t65*t2*t30+2.0*t62*t70-t61*t73*t15-t61*(t36-t39+t25*t44*t18-t25*t47*t50)*t15-2.0*t61*t65*t70-t61*t2*t40*t28*t29;  
  
  

}







if(BivD %in% c("J0","J90","J180","J270")){


t1 = 1.0-p1;
t2 = t1^(1.0*teta);
t3 = 1.0-p2;
t4 = t3^(1.0*teta);
t5 = t2*t4;
t6 = t2+t4-t5;
t8 = 1/teta-2.0;
t9 = t6^(1.0*t8);
t10 = teta^2;
t11 = log(t6);
t12 = log(t1);
t13 = t2*t12;
t14 = log(t3);
t15 = t4*t14;
t16 = t13*t4;
t19 = t5*t14;
t21 = teta-1.0;
t27 = t1^(1.0*t21);
t28 = t3^(1.0*t21);
t30 = teta-1.0+t2+t4-t5;
t33 = t9*t27;

der2h.derp1teta <- t9*(-1/t10*t11+t8*(t13+t15-t16-t19)/t6)*t27*t28*t30+t33*t12*t28*t30+t33*t28*t14*t30+t33*t28*(1.0+t13+t15-t16-t19);
            


                        t1 = 1.0-p1
			t2 = t1^teta
			t3 = 1.0-p2
			t4 = t3^teta
			t5 = t2*t4
			t6 = t2+t4-t5
			t8 = 1/teta-2.0
			t9 = t6^t8
			t10 = teta^2
			t11 = 1/t10
			t12 = log(t6)
			t14 = log(t1)
			t15 = t2*t14
			t16 = log(t3)
			t17 = t4*t16
			t18 = t15*t4
			t19 = t5*t16
			t20 = t15+t17-t18-t19
			t22 = 1/t6
			t24 = -t11*t12+t8*t20*t22
			t25 = t24*t24
			t27 = teta-1.0
			t28 = t1^t27
			t29 = t3^t27
			t30 = t28*t29
			t31 = teta-1.0+t2+t4-t5
			t32 = t30*t31
			t41 = t14*t14
			t42 = t2*t41
			t43 = t16*t16
			t49 = t42+t4*t43-t42*t4-2.0*t15*t17-t5*t43
			t51 = t20*t20
			t53 = t6*t6
			t60 = t9*t24
			t61 = t60*t28
			t62 = t14*t29
			t66 = t29*t16
			t67 = t66*t31
			t70 = 1.0+t15+t17-t18-t19
			t74 = t9*t28
der2c.derrho.derrho <- t9*t25*t32+t9*(2.0/t10/teta*t12-2.0*t11*t20*t22+t8*t49*t22-t8*t51/t53)*t32+2.0*t61*t62*t31+2.0*t61*t67+2.0*t60*t30*t70+t74*t41*t29*t31+2.0*t74*t14*t67+2.0*t74*t62*t70+t74*t29*t43*t31+2.0*t74*t66*t70+t74*t29*t49


 
 
 
 			  t1 = 1.0-p1;
 			  t2 = t1^teta
 			  t3 = 1.0-p2;
 			  t4 = t3^teta;
 			  t5 = t2*t4;
 			  t6 = t2+t4-t5;
 			  t8 = 1/teta-1.0;
 			  t9 = t6^t8;
 			  t10 = teta^2
 			  t11 = 1/t10;
 			  t12 = log(t6);
 			  t14 = log(t1);
 			  t15 = t2*t14;
 			  t16 = log(t3);
 			  t18 = t4*t16;
 			  t20 = t15+t18-t15*t4-t5*t16;
 			  t21 = 1/t6;
 			  t23 = -t11*t12+t8*t20*t21;
 			  t25 = t23^2;
 			  t28 = t3^(teta-1)
 			  t29 = 1 - t2;
 			  t30 = t28*t29;
 			  t39 = t14^2;
 			  t40 = t2*t39;
 			  t42 = t16^2;
 			  t50 = t20^2;
 			  t56 = t6^2;
 			  t13 = t9*t23;
 			  t7 = t9*t28;
 			  
 			  
 der2h.derteta.teta.st <- t9*t25*t30+t9*(2.0/t10/teta*t12-2.0*t11*t20*t21+t8*(t40+t4*t42-t40*t4-2.0*t15*t18-t5*t42)*t21-t8*t50/t56)*t30+2.0*t13*t28*t16*t29-2.0*t13*t28*t2*t14+t7*t42*t29-2.0*t7*t16*t2*t14-t7*t40;

}





if(BivD %in% c("C90","J90","G90","GAL90") ) {


der2h.derp1teta       <- -der2h.derp1teta
der2h.derteta.teta.st <- -der2h.derteta.teta.st 


}

if(BivD %in% c("C180","J180","G180","GAL180") ) {

derteta.derteta.st         = -derteta.derteta.st
der2teta.derteta.stteta.st = -der2teta.derteta.stteta.st 

der2h.derp1teta            = -der2h.derp1teta      
der2h.derteta.teta.st      = -der2h.derteta.teta.st



}




if(BivD %in% c("C270","J270","G270","GAL270") ) {

der2h.derp1teta <- -der2h.derp1teta



}




list(
derteta.derteta.st         = ifef(derteta.derteta.st), 
der2h.derteta.teta.st      = ifef(der2h.derteta.teta.st),  
der2teta.derteta.stteta.st = ifef(der2teta.derteta.stteta.st),  
der2h.derp1teta            = ifef(der2h.derp1teta),  
der2c.derrho.derrho        = ifef(der2c.derrho.derrho)
   
)     



}

