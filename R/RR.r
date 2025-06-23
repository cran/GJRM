RR <- function(x, trt, trt.val = NULL, int.var = NULL, joint = TRUE, n.sim = 100, prob.lev = 0.05, length.out = NULL){

if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous/discrete margins.")
if(x$Cont == "NO" && x$VC$ccss == "yes" ) stop("This function is not suitable for selection models.")

if(joint == TRUE)  type <- "joint"
if(joint == FALSE) type <- "univariate"

lbn <- paste(prob.lev/2*100, "%", sep = "")
ubn <- paste((1-(prob.lev/2))*100, "%", sep = "")

CIs <- est.AT <- NULL

etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- delta.AT <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- NULL

m2  <- x$VC$m2 
m3  <- x$VC$m3 
bin.link <- x$VC$bl  
end <- 0
est.ATb <- NA
indD <- list()

if(x$v1[1] %in% x$v2[-1]) {end <- 1; eq <- 2} 
if(x$v2[1] %in% x$v1[-1]) {end <- 2; eq <- 1}


if( !( type %in% c("univariate","joint") ) ) stop("Error in parameter type value. It should be univariate or joint.")
if( !(x$margins[2] %in% bin.link) && eq == 2) stop("Calculation of RR not allowed for.")
if( !(x$margins[1] %in% bin.link && x$margins[2] %in% c("probit", "logit", "cloglog", "N","tN")) ) stop("The margins have to be probit, logit, cloglog or (truncated) Gaussian for this measure to make sense.")
if(missing(trt)) stop("You must provide the name of the treatment variable.")

if(type == "naive" && !(x$margins[2] %in% bin.link)) stop("Please fit a bivariate model with intercept and endogenous variable only and then use RR with the univariate type option.")

if(x$Model=="BSS" || x$Model=="BPO" || x$Model=="BPO0" || end==0) stop("Calculation of this effect is valid for recursive models only.")
if(is.character(trt)==FALSE) stop("trt is not a character!")

 
 
 if( !is.null(int.var) && !is.list(int.var) ){
 
   if( length(int.var) != 2 )              stop("int.var must contain a name and a value for the interaction term.")
   if( is.character(int.var[1]) == FALSE ) stop("The first element of int.var must be the name of the interaction term.")
 
   int.var1 <- int.var[1]
   int.var2 <- as.numeric(int.var[2]) # as.numeric works for both numeric and factor vars
 
   if( !(int.var2 %in% c(0, 1)) ) stop("The interaction can only currently take value 0 or 1.")  
 
 }
 
 
 if( !is.null(int.var) && is.list(int.var) ){
 
 int.var1 <- int.var2 <- NA
 
  for(i in 1:length(int.var)){
 
   if( length(int.var[[i]]) != 2 )              stop("Each of the int.var objects in the list must contain a name and a value for the interaction term.")
   if( is.character(int.var[[i]][1]) == FALSE ) stop("The first element of each of the int.var objects must be the name of the interaction term.")
 
   int.var1[i] <- int.var[[i]][1]
   int.var2[i] <- as.numeric(int.var[[i]][2]) # as.numeric works for both numeric and factor vars
 
   if( !(int.var2[i] %in% c(0, 1)) ) stop("The interaction term can only currently take value 0 or 1.")  
  }
 
 
 }
 
 


##############################################################################################################

if(type == "naive" && x$margins[2] %in% bin.link){ # it looks ok from comparing the naive and univariate options

if(eq==2){
y1 <- x$y1
y2 <- x$y2
}

if(eq==1){
y1 <- x$y2 
y2 <- x$y1
}

tab2 <- table(y1, y2)                                  


pY1cT1 <- prop.table(tab2,1)[4] 
pY1cT0 <- prop.table(tab2,1)[3] 

est.AT <- (pY1cT1 / pY1cT0)

sv <- qnorm(prob.lev/2,lower.tail = FALSE) * sqrt( (1-pY1cT1)/(table(y1)[2]*pY1cT1) + (1-pY1cT0)/(table(y1)[1]*pY1cT0) )

CIs <- exp( c(log(est.AT) - sv, log(est.AT) + sv) )


est.ATb <- est.ATso <- NULL

}



##############################################################################################################



if(type != "naive" && x$margins[2] %in% bin.link) {


########################################################
# Set-up
########################################################



if(type == "joint"){

	indD[[1]] <- 1:x$X1.d2 
	indD[[2]] <- x$X1.d2+(1:x$X2.d2)
	
}


if(eq==1){ 
           ff <- reformulate(all.vars(x$gam1$terms)[-1]); tgam <- x$gam1
           if(type == "joint") ind.int <- indD[[1]]  
}

if(eq==2){ 
           ff <- reformulate(all.vars(x$gam2$terms)[-1]); tgam <- x$gam2
           if(type == "joint") ind.int <- indD[[2]]
}



d0 <- d1 <- model.frame(ff, data = eval(x$mcd), subset = eval(x$sbs))
attr(d0,"terms") <- attr(d1,"terms") <- NULL








if(type == "joint") coef.int <- x$coefficients[ind.int]
	   

if( is.logical(d0[, trt]) == TRUE) stop("The treatment variable must be a binary numeric or factor variable.")
             
d0[, trt] <- 0
d1[, trt] <- 1

d0 <- predict(tgam, d0, type = "lpmatrix")  
d1 <- predict(tgam, d1, type = "lpmatrix")  



if( !is.null(int.var) && !is.list(int.var)) {

   if( any(grepl(int.var1, dimnames(d1)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1)) == FALSE      )          stop("Check the name provided for the interaction term.")

if( int.var2 == 0) d1[, int.var1] <- 0
if( int.var2 == 1) d1[, int.var1] <- 1

}


if( !is.null(int.var) && is.list(int.var)) {

 for(i in 1:length(int.var)){

   if( any(grepl(int.var1[i], dimnames(d1)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1[i])) == FALSE      )          stop("Check the name provided for the interaction term.")

if( int.var2[i] == 0) d1[, int.var1[i]] <- 0
if( int.var2[i] == 1) d1[, int.var1[i]] <- 1

 }

}




if(type == "joint"){

	eti1 <- d1%*%coef.int 
	eti0 <- d0%*%coef.int 
  
}

if(type == "univariate"){

if(eq==1) ngam <- x$gam1
if(eq==2) ngam <- x$gam2

        eti1 <- d1%*%ngam$coefficients 
        eti0 <- d0%*%ngam$coefficients 

}




#############################################################################
# RR
#############################################################################

p.int1 <- probm(eti1, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 
p.int0 <- probm(eti0, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr  

est.AT <- mean(p.int1, na.rm = TRUE) / mean(p.int0, na.rm = TRUE) 
                        
#############################################################################
# CIs RR
#############################################################################


 if(type == "univariate") {bs <- rMVN(n.sim, mean = ngam$coefficients, sigma=ngam$Vp)
                           eti1s <- d1%*%t(bs)
                           eti0s <- d0%*%t(bs) 
                           }
 if(type == "joint")  {bs <- rMVN(n.sim, mean = x$coefficients, sigma=x$Vb)
                           eti1s <- d1%*%t(bs[,ind.int])
                           eti0s <- d0%*%t(bs[,ind.int]) 
                           } 


 peti1s <- probm(eti1s, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr   
 peti0s <- probm(eti0s, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 

 est.ATb <- colMeans(peti1s, na.rm = TRUE) / colMeans(peti0s, na.rm = TRUE) 


CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2,1-prob.lev/2), na.rm = TRUE))



                     
                     
 
}


















###############################################


if(type != "naive" && !(x$margins[2] %in% bin.link)) {

if(is.null(trt.val)) stop("You must provide a value for the endogenous variable.")


ff <- reformulate(all.vars(x$gam1$terms)[-1])
d0 <- d1 <- model.frame(ff, data = eval(x$mcd), subset = eval(x$sbs))
attr(d0,"terms") <- attr(d1,"terms") <- NULL

d0[, trt] <- trt.val
d1[, trt] <- trt.val + 1

d0 <- predict(x$gam1, d0, type = "lpmatrix")  
d1 <- predict(x$gam1, d1, type = "lpmatrix") 


if( !is.null(int.var) && !is.list(int.var)) {

   if( any(grepl(int.var1, dimnames(d1)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1)) == FALSE      )          stop("Check the name provided for the interaction term.")

if( int.var2 == 0) d0[, int.var1] <- d1[, int.var1] <- 0
if( int.var2 == 1) d0[, int.var1] <- d1[, int.var1] <- 1

}




if( !is.null(int.var) && is.list(int.var)) {

 for(i in 1:length(int.var)){

   if( any(grepl(int.var1[i], dimnames(d1)[[2]])) == FALSE ) stop("Check the name provided for the interaction term.")
   if( any(grepl(":", int.var1[i])) == FALSE      )          stop("Check the name provided for the interaction term.")

if( int.var2[i] == 0) d0[, int.var1[i]] <- d1[, int.var1[i]] <- 0
if( int.var2[i] == 1) d0[, int.var1[i]] <- d1[, int.var1[i]] <- 1

 }

}




ind.int <- 1:x$X1.d2
ngam <- x$gam1


if(type == "joint"){
	eti1 <- d1%*%x$coefficients[ind.int] 
	eti0 <- d0%*%x$coefficients[ind.int] 
                    }

if(type == "univariate"){
	eti1 <- d1%*%x$gam1$coefficients 
	eti0 <- d0%*%x$gam1$coefficients
                         }

#############################################################################
# RR
#############################################################################

eq <- 1

p.int1 <- probm(eti1, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 
p.int0 <- probm(eti0, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr  

est.AT <- mean(p.int1, na.rm = TRUE) / mean(p.int0, na.rm = TRUE) 
                        
#############################################################################
# CIs RR
#############################################################################


 if(type == "univariate") {bs <- rMVN(n.sim, mean = ngam$coefficients, sigma=ngam$Vp)
                           eti1s <- d1%*%t(bs)
                           eti0s <- d0%*%t(bs) 
                           }
 if(type == "joint")  {bs <- rMVN(n.sim, mean = x$coefficients, sigma=x$Vb)
                           eti1s <- d1%*%t(bs[,ind.int])
                           eti0s <- d0%*%t(bs[,ind.int]) 
                           } 


 peti1s <- probm(eti1s, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr   
 peti0s <- probm(eti0s, x$margins[eq], min.dn = x$VC$min.dn, min.pr = x$VC$min.pr, max.pr = x$VC$max.pr)$pr 

 est.ATb <- colMeans(peti1s, na.rm = TRUE) / colMeans(peti0s, na.rm = TRUE) 


CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2,1-prob.lev/2), na.rm = TRUE))



}



####################



res <- c(CIs[1], est.AT, CIs[2])
if(!(type != "naive" && !(x$margins[2] %in% bin.link)))   names(res) <- c(lbn, "RR", ubn)



rm(etap.noi, X.int, X.noi, eti1, eti0, etno, indS, bs, ind.excl, p.int1, p.int0, d.int1, d.int0,
   p.etn, d.etn, ass.p, ass.pst, C.11, C.10, sig2, peti1s, peti0s, sigma2.st, sigma2s, eti1s, eti0s, d0, d1,
   p.etns, etnos, etds, ass.ps)  
   
   
out <- list(res=res, prob.lev=prob.lev, sim.RR=est.ATb, mar2=x$margins[2], type = type,
            Ratios = Effects, treat = y2, eq = eq, bl = x$VC$bl) # Pr = Pr,
 
  
 
class(out) <- "RR"

out





}





















