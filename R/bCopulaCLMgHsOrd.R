bCopulaCLMgHsOrd <- function(params, respvec, VC, ps, AT = FALSE) {

# This function is based on bCopulaCLMgHsCon

p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA

n <- VC$n

# With the cut points we might need an extra parameter to be counted in the dimensions below: K1 and K2 will be 
# the dimension of the cut-points

c1.ti <- params[1:(VC$K1-1)]
c2.ti <- params[VC$K1 : (VC$K1 + VC$K2 - 2)]

eta1 <- VC$X1%*%params[(VC$K1 + VC$K2 - 1) : (VC$K1 + VC$K2 + VC$X1.d2 - 2)]
eta2 <- VC$X2%*%params[(VC$K1 + VC$K2 + VC$X1.d2 - 1) : (VC$K1 + VC$K2 + VC$X1.d2 + VC$X2.d2 - 2)]
etad <- l.ln <- NULL

epsilon <-  sqrt(.Machine$double.eps) # RECALL: results are sensible to the choice of epsilon (try also 0.00001)

if (   is.null(VC$X3)  ) { teta.st <- etad <- params[VC$K1 + VC$K2 + VC$X1.d2 + VC$X2.d2 - 1] }
if ( !(is.null(VC$X3)) ) { teta.st <- etad <- VC$X3%*%params[(VC$K1 + VC$K2 + VC$X1.d2 + VC$X2.d2 - 1) : (VC$K1 + VC$K2 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 - 2)] }


##############################

# Independence model

if (VC$ind.ord == "TRUE") teta.st <- etad <- log(epsilon) # This ensures that teta = 1 (i.e. independence for J0 copula)

##############################


# Cut points are transformed and the linear predictors created 

infty  <- 1e+25 ; e <- .Machine$double.eps * 10^6

c1  <- VC$c1 ; c1[1] <- c1.ti[1] ; for (i in 2 : (VC$K1 - 1)) c1[i] <- c1[i - 1] + c1.ti[i]^2
c1  <- sign(c1) * pmin(10000 * infty, abs(c1))

c2  <- VC$c2 ; c2[1] <- c2.ti[1] ; for (i in 2 : (VC$K2 - 1)) c2[i] <- c2[i - 1] + c2.ti[i]^2
c2  <- sign(c2) * pmin(10000 * infty, abs(c2))

c1.m   <- t(matrix(nrow = VC$K1 - 1, ncol = n        , c1  )) 
eta1.m <-   matrix(nrow = n        , ncol = VC$K1 - 1, eta1)

c2.m   <- t(matrix(nrow = VC$K2 - 1, ncol = n        , c2  )) 
eta2.m <-   matrix(nrow = n        , ncol = VC$K2 - 1, eta2)

lp1 <- c1.m - eta1.m ; lp1.p <- cbind(lp1, infty) ; lp1.m <- cbind(-infty, lp1)
lp2 <- c2.m - eta2.m ; lp2.p <- cbind(lp2, infty) ; lp2.m <- cbind(-infty, lp2)

sel1 <- VC$sel1 ; lp1.p.sel <- rowSums(lp1.p * sel1) ; lp1.m.sel <- rowSums(lp1.m * sel1)
sel2 <- VC$sel2 ; lp2.p.sel <- rowSums(lp2.p * sel2) ; lp2.m.sel <- rowSums(lp2.m * sel2)

lp1.p <- lp1.p.sel ; lp1.m <- lp1.m.sel
lp2.p <- lp2.p.sel ; lp2.m <- lp2.m.sel

##############################

####

#sstr1 <- esp.tr(sigma2.st, VC$margins[2]) #
#sigma2.st <- sstr1$vrb.st                 #
#sigma2    <- sstr1$vrb                    # In the ordinal-ordinal model there's not sd/variance for the margin

####

pd1.p <- probm(lp1.p, VC$margins[1], only.pr = FALSE, bc = TRUE, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) # I set only.pr = FALSE to obtain d.n 
pd1.m <- probm(lp1.m, VC$margins[1], only.pr = FALSE, bc = TRUE, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) #

pd2.p <- probm(lp2.p, VC$margins[2], only.pr = FALSE, bc = TRUE, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
pd2.m <- probm(lp2.m, VC$margins[2], only.pr = FALSE, bc = TRUE, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

p1.p <- pd1.p$pr ; p1.m <- pd1.m$pr #
p2.p <- pd2.p$pr ; p2.m <- pd2.m$pr # These are cdfs

####

# All the resT, teta.st1, teta1

resT <- teta.tr(VC, teta.st)

teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
teta1 <- teta2 <- teta <- resT$teta

#### 

Cop1 <- Cop2 <- VC$BivD 

teta.ind1 <- as.logical(c(1, 0, round(runif(VC$n-2)))) 
teta.ind2 <- teta.ind1 == FALSE 


if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

teta1 <- teta[teta.ind1]
teta2 <- teta[teta.ind2]

}


if(VC$BivD %in% VC$BivD2){

if(VC$BivD %in% VC$BivD2[c(1:4,13:16) ]) teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov)    , TRUE, FALSE)
if(VC$BivD %in% VC$BivD2[5 : 12]       ) teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov) + 1, TRUE, FALSE) 
teta.ind2 <- teta.ind1 == FALSE 

VC$my.env$signind <- ifelse(teta.ind1 == TRUE,  1, -1) 

teta1 <-   teta[teta.ind1]
teta2 <- - teta[teta.ind2]

teta.st1 <- teta.st[teta.ind1]
teta.st2 <- teta.st[teta.ind2]

if(length(teta) == 1) teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)  

Cop1Cop2R <- Cop1Cop2(VC$BivD)
Cop1 <- Cop1Cop2R$Cop1
Cop2 <- Cop1Cop2R$Cop2

} 


##### Log-likelihood function ##### 

h.pp <- h.pm <- h.mp <- h.mm <- NA

if (length(teta1) != 0) {

dH1.pp <- copgHs(p1.p[teta.ind1], p2.p[teta.ind1], eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
dH1.pm <- copgHs(p1.p[teta.ind1], p2.m[teta.ind1], eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
dH1.mp <- copgHs(p1.m[teta.ind1], p2.p[teta.ind1], eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
dH1.mm <- copgHs(p1.m[teta.ind1], p2.m[teta.ind1], eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

h.pp[teta.ind1] <- BiCDF(u1 = p1.p[teta.ind1], u2 = p2.p[teta.ind1], family = VC$nC, par1 = teta1, par2 = NULL, test = TRUE)
h.pm[teta.ind1] <- BiCDF(u1 = p1.p[teta.ind1], u2 = p2.m[teta.ind1], family = VC$nC, par1 = teta1, par2 = NULL, test = TRUE)
h.mp[teta.ind1] <- BiCDF(u1 = p1.m[teta.ind1], u2 = p2.p[teta.ind1], family = VC$nC, par1 = teta1, par2 = NULL, test = TRUE)
h.mm[teta.ind1] <- BiCDF(u1 = p1.m[teta.ind1], u2 = p2.m[teta.ind1], family = VC$nC, par1 = teta1, par2 = NULL, test = TRUE)

}

if (length(teta2) != 0) {

dH2.pp <- copgHs(p1.p[teta.ind2], p2.p[teta.ind2], eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
dH2.pm <- copgHs(p1.p[teta.ind2], p2.m[teta.ind2], eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
dH2.mp <- copgHs(p1.m[teta.ind2], p2.p[teta.ind2], eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
dH2.mm <- copgHs(p1.m[teta.ind2], p2.m[teta.ind2], eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, CLM = TRUE, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)

h.pp[teta.ind2] <- BiCDF(u1 = p1.p[teta.ind2], u2 = p2.p[teta.ind2], family = VC$nC, par1 = teta2, par2 = NULL, test = TRUE)
h.pm[teta.ind2] <- BiCDF(u1 = p1.p[teta.ind2], u2 = p2.m[teta.ind2], family = VC$nC, par1 = teta2, par2 = NULL, test = TRUE)
h.mp[teta.ind2] <- BiCDF(u1 = p1.m[teta.ind2], u2 = p2.p[teta.ind2], family = VC$nC, par1 = teta2, par2 = NULL, test = TRUE)
h.mm[teta.ind2] <- BiCDF(u1 = p1.m[teta.ind2], u2 = p2.m[teta.ind2], family = VC$nC, par1 = teta2, par2 = NULL, test = TRUE)

}

h.2p2m <- pmax(epsilon, h.pp - h.pm - h.mp + h.mm)

l.par <- VC$weights * log(h.2p2m)


##### Score vector ##### 

c.copula.be1.pp <- c.copula.be1.pm <- c.copula.be1.mp <- c.copula.be1.mm <- 
c.copula.be2.pp <- c.copula.be2.pm <- c.copula.be2.mp <- c.copula.be2.mm <- 
c.copula.th.pp  <- c.copula.th.pm  <- c.copula.th.mp  <- c.copula.th.mm  <- NA


if(length(teta1) != 0){

c.copula.be1.pp[teta.ind1] <- dH1.pp$c.copula.be1   ; c.copula.be1.pm[teta.ind1] <- dH1.pm$c.copula.be1
c.copula.be1.mp[teta.ind1] <- dH1.mp$c.copula.be1   ; c.copula.be1.mm[teta.ind1] <- dH1.mm$c.copula.be1
c.copula.be2.pp[teta.ind1] <- dH1.pp$c.copula.be2   ; c.copula.be2.pm[teta.ind1] <- dH1.pm$c.copula.be2
c.copula.be2.mp[teta.ind1] <- dH1.mp$c.copula.be2   ; c.copula.be2.mm[teta.ind1] <- dH1.mm$c.copula.be2
c.copula.th.pp[teta.ind1]  <- dH1.pp$c.copula.theta ; c.copula.th.pm[teta.ind1]  <- dH1.pm$c.copula.theta
c.copula.th.mp[teta.ind1]  <- dH1.mp$c.copula.theta ; c.copula.th.mm[teta.ind1]  <- dH1.mm$c.copula.theta

}

if(length(teta2) != 0){

c.copula.be1.pp[teta.ind2] <- dH2.pp$c.copula.be1   ; c.copula.be1.pm[teta.ind2] <- dH2.pm$c.copula.be1
c.copula.be1.mp[teta.ind2] <- dH2.mp$c.copula.be1   ; c.copula.be1.mm[teta.ind2] <- dH2.mm$c.copula.be1
c.copula.be2.pp[teta.ind2] <- dH2.pp$c.copula.be2   ; c.copula.be2.pm[teta.ind2] <- dH2.pm$c.copula.be2
c.copula.be2.mp[teta.ind2] <- dH2.mp$c.copula.be2   ; c.copula.be2.mm[teta.ind2] <- dH2.mm$c.copula.be2
c.copula.th.pp[teta.ind2]  <- dH2.pp$c.copula.theta ; c.copula.th.pm[teta.ind2]  <- dH2.pm$c.copula.theta
c.copula.th.mp[teta.ind2]  <- dH2.mp$c.copula.theta ; c.copula.th.mm[teta.ind2]  <- dH2.mm$c.copula.theta

}

c.copula.th.2p2m <- c.copula.th.pp - c.copula.th.pm - c.copula.th.mp + c.copula.th.mm

####

cc1.ti <- matrix(nrow = n, ncol = VC$K1 - 2, byrow = TRUE, c1.ti[2 : (VC$K1 - 1)])
cc2.ti <- matrix(nrow = n, ncol = VC$K2 - 2, byrow = TRUE, c2.ti[2 : (VC$K2 - 1)])

################# These quantities are extracted from CopulaCLM()

sel1.p <- VC$sel1.p ; sel1.m <- VC$sel1.m
sel2.p <- VC$sel2.p ; sel2.m <- VC$sel2.m

sel1.mm <- VC$sel1.mm ; sel1.pm <- VC$sel1.pm
sel2.mm <- VC$sel2.mm ; sel2.pm <- VC$sel2.pm

D11 <- VC$D11 ; D12 <- VC$D12
D21 <- VC$D21 ; D22 <- VC$D22

#################

D11.1 <- 2 * sel1.m * cc1.ti
D12.1 <- 2 * sel1.p * cc1.ti

D21.1 <- 2 * sel2.m * cc2.ti
D22.1 <- 2 * sel2.p * cc2.ti

#################

derp1.dereta1.p <- - pd1.p$derp1.dereta1 ; derp1.dereta1.m <- - pd1.m$derp1.dereta1 # RECALL: The function probm was originally constructed for 1 - pd1.p$pr. I use instead pd1.p$pr
derp2.dereta2.p <- - pd2.p$derp1.dereta1 ; derp2.dereta2.m <- - pd2.m$derp1.dereta1 # 

derh.dereta1.p  <- (c.copula.be1.pp - c.copula.be1.pm) * derp1.dereta1.p ; derh.dereta1.m <- (c.copula.be1.mp - c.copula.be1.mm) * derp1.dereta1.m
derh.dereta2.p  <- (c.copula.be2.pp - c.copula.be2.mp) * derp2.dereta2.p ; derh.dereta2.m <- (c.copula.be2.pm - c.copula.be2.mm) * derp2.dereta2.m

derh.dereta1.pm <- derh.dereta1.p - derh.dereta1.m
derh.dereta2.pm <- derh.dereta2.p - derh.dereta2.m


# Approximating quantities for the score 

c.copula.th.2p2m_app <- approx.CLM(c.copula.th.2p2m, h.2p2m, epsilon)

derh.dereta1.pm_app <- approx.CLM(derh.dereta1.pm, h.2p2m, epsilon)
derh.dereta2.pm_app <- approx.CLM(derh.dereta2.pm, h.2p2m, epsilon)


dl.dc1_1 <- (D11 * derh.dereta1.p - D12 * derh.dereta1.m) / h.2p2m
dl.dc1_2 <- (D11.1 * c(derh.dereta1.p) - D12.1 * c(derh.dereta1.m)) / h.2p2m
dl.dc2_1 <- (D21 * derh.dereta2.p - D22 * derh.dereta2.m) / h.2p2m
dl.dc2_2 <- (D21.1 * c(derh.dereta2.p) - D22.1 * c(derh.dereta2.m)) / h.2p2m

dl.dbe1 <- derh.dereta1.pm_app / h.2p2m
dl.dbe2 <- derh.dereta2.pm_app / h.2p2m

dl.dteta.st <- c.copula.th.2p2m_app / h.2p2m


# Adding the weigths

dl.dc1_1_w    <- VC$weights * dl.dc1_1
dl.dc1_2_w    <- VC$weights * dl.dc1_2
dl.dc2_1_w    <- VC$weights * dl.dc2_1
dl.dc2_2_w    <- VC$weights * dl.dc2_2
dl.dbe1_w     <- VC$weights * dl.dbe1
dl.dbe2_w     <- VC$weights * dl.dbe2
dl.dteta.st_w <- VC$weights * dl.dteta.st 


##### Hessian matrix #####

c.copula2.be1be1.pp <- c.copula2.be1be1.pm <- c.copula2.be1be1.mp <- c.copula2.be1be1.mm <-
c.copula2.be1be2.pp <- c.copula2.be1be2.pm <- c.copula2.be1be2.mp <- c.copula2.be1be2.mm <-
c.copula2.be1th.pp  <- c.copula2.be1th.pm  <- c.copula2.be1th.mp  <- c.copula2.be1th.mm  <-
c.copula2.be2be2.pp <- c.copula2.be2be2.pm <- c.copula2.be2be2.mp <- c.copula2.be2be2.mm <-
c.copula2.be2th.pp  <- c.copula2.be2th.pm  <- c.copula2.be2th.mp  <- c.copula2.be2th.mm  <-
c.copula2.th2.pp <- c.copula2.th2.pm <- c.copula2.th2.mp <- c.copula2.th2.mm <- NA

derteta.derteta.st      <- der2teta.derteta.stteta.st <- NA

der2p1.dereta1eta1.p <- - pd1.p$der2p1.dereta1eta1 # RECALL: The function probm was originally constructed for 1 - pd1.p$pr. I use instead pd1.p$pr
der2p1.dereta1eta1.m <- - pd1.m$der2p1.dereta1eta1 #
der2p2.dereta2eta2.p <- - pd2.p$der2p1.dereta1eta1 #
der2p2.dereta2eta2.m <- - pd2.m$der2p1.dereta1eta1 #


if(length(teta1) != 0){

c.copula2.be1be1.pp[teta.ind1] <- dH1.pp$c.copula2.be1    ; c.copula2.be1be1.pm[teta.ind1] <- dH1.pm$c.copula2.be1
c.copula2.be1be1.mp[teta.ind1] <- dH1.mp$c.copula2.be1    ; c.copula2.be1be1.mm[teta.ind1] <- dH1.mm$c.copula2.be1
c.copula2.be1be2.pp[teta.ind1] <- dH1.pp$c.copula2.be1be2 ; c.copula2.be1be2.pm[teta.ind1] <- dH1.pm$c.copula2.be1be2
c.copula2.be1be2.mp[teta.ind1] <- dH1.mp$c.copula2.be1be2 ; c.copula2.be1be2.mm[teta.ind1] <- dH1.mm$c.copula2.be1be2
c.copula2.be1th.pp[teta.ind1]  <- dH1.pp$c.copula2.be1th  ; c.copula2.be1th.pm[teta.ind1]  <- dH1.pm$c.copula2.be1th
c.copula2.be1th.mp[teta.ind1]  <- dH1.mp$c.copula2.be1th  ; c.copula2.be1th.mm[teta.ind1]  <- dH1.mm$c.copula2.be1th

c.copula2.be2be2.pp[teta.ind1] <- dH1.pp$c.copula2.be2   ; c.copula2.be2be2.pm[teta.ind1] <- dH1.pm$c.copula2.be2
c.copula2.be2be2.mp[teta.ind1] <- dH1.mp$c.copula2.be2   ; c.copula2.be2be2.mm[teta.ind1] <- dH1.mm$c.copula2.be2
c.copula2.be2th.pp[teta.ind1]  <- dH1.pp$c.copula2.be2th ; c.copula2.be2th.pm[teta.ind1]  <- dH1.pm$c.copula2.be2th
c.copula2.be2th.mp[teta.ind1]  <- dH1.mp$c.copula2.be2th ; c.copula2.be2th.mm[teta.ind1]  <- dH1.mm$c.copula2.be2th

c.copula2.th2.pp[teta.ind1] <- dH1.pp$bit1.th2ATE ; c.copula2.th2.pm[teta.ind1] <- dH1.pm$bit1.th2ATE
c.copula2.th2.mp[teta.ind1] <- dH1.mp$bit1.th2ATE ; c.copula2.th2.mm[teta.ind1] <- dH1.mm$bit1.th2ATE

derteta.derteta.st[teta.ind1]         <- dH1.pp$derteta.derteta.st # RECALL: this quantity does not depend on the margins, so the use of BITS.p and BITS.m is indifferent 
der2teta.derteta.stteta.st[teta.ind1] <- dH1.pp$der2teta.derteta.stteta.st

}

if(length(teta2) != 0){

c.copula2.be1be1.pp[teta.ind2] <- dH2.pp$c.copula2.be1    ; c.copula2.be1be1.pm[teta.ind2] <- dH2.pm$c.copula2.be1
c.copula2.be1be1.mp[teta.ind2] <- dH2.mp$c.copula2.be1    ; c.copula2.be1be1.mm[teta.ind2] <- dH2.mm$c.copula2.be1
c.copula2.be1be2.pp[teta.ind2] <- dH2.pp$c.copula2.be1be2 ; c.copula2.be1be2.pm[teta.ind2] <- dH2.pm$c.copula2.be1be2
c.copula2.be1be2.mp[teta.ind2] <- dH2.mp$c.copula2.be1be2 ; c.copula2.be1be2.mm[teta.ind2] <- dH2.mm$c.copula2.be1be2
c.copula2.be1th.pp[teta.ind2]  <- dH2.pp$c.copula2.be1th  ; c.copula2.be1th.pm[teta.ind2]  <- dH2.pm$c.copula2.be1th
c.copula2.be1th.mp[teta.ind2]  <- dH2.mp$c.copula2.be1th  ; c.copula2.be1th.mm[teta.ind2]  <- dH2.mm$c.copula2.be1th

c.copula2.be2be2.pp[teta.ind2] <- dH2.pp$c.copula2.be2   ; c.copula2.be2be2.pm[teta.ind2] <- dH2.pm$c.copula2.be2
c.copula2.be2be2.mp[teta.ind2] <- dH2.mp$c.copula2.be2   ; c.copula2.be2be2.mm[teta.ind2] <- dH2.mm$c.copula2.be2
c.copula2.be2th.pp[teta.ind2]  <- dH2.pp$c.copula2.be2th ; c.copula2.be2th.pm[teta.ind2]  <- dH2.pm$c.copula2.be2th
c.copula2.be2th.mp[teta.ind2]  <- dH2.mp$c.copula2.be2th ; c.copula2.be2th.mm[teta.ind2]  <- dH2.mm$c.copula2.be2th

c.copula2.th2.pp[teta.ind2] <- dH2.pp$bit1.th2ATE ; c.copula2.th2.pm[teta.ind2] <- dH2.pm$bit1.th2ATE
c.copula2.th2.mp[teta.ind2] <- dH2.mp$bit1.th2ATE ; c.copula2.th2.mm[teta.ind2] <- dH2.mm$bit1.th2ATE

derteta.derteta.st[teta.ind2]         <- dH2.pp$derteta.derteta.st # RECALL: this quantity does not depend on the margins, so the use of BITS.pp and BITS.mp etc. is indifferent 
der2teta.derteta.stteta.st[teta.ind2] <- dH2.pp$der2teta.derteta.stteta.st

}

c.copula2.be1be1.2pm <- c.copula2.be1be1.pp - c.copula2.be1be1.pm
c.copula2.be1be1.p2m <- c.copula2.be1be1.mp - c.copula2.be1be1.mm
c.copula2.be1th.2pm  <- c.copula2.be1th.pp  - c.copula2.be1th.pm
c.copula2.be1th.p2m  <- c.copula2.be1th.mp  - c.copula2.be1th.mm

c.copula2.be2be2.2pm <- c.copula2.be2be2.pp - c.copula2.be2be2.mp
c.copula2.be2be2.p2m <- c.copula2.be2be2.pm - c.copula2.be2be2.mm
c.copula2.be2th.2pm  <- c.copula2.be2th.pp  - c.copula2.be2th.mp
c.copula2.be2th.p2m  <- c.copula2.be2th.pm  - c.copula2.be2th.mm

c.copula2.th2.2p2m <- c.copula2.th2.pp - c.copula2.th2.pm - c.copula2.th2.mp + c.copula2.th2.mm

der2h.derc12.derc12.2pm <- c.copula2.be1be1.2pm * (derp1.dereta1.p)^2 + (c.copula.be1.pp - c.copula.be1.pm) * der2p1.dereta1eta1.p
der2h.derc12.derc12.p2m <- c.copula2.be1be1.p2m * (derp1.dereta1.m)^2 + (c.copula.be1.mp - c.copula.be1.mm) * der2p1.dereta1eta1.m
der2h.derc12.derc12.mag <- sel1.m * (c.copula.be1.pp - c.copula.be1.pm) * derp1.dereta1.p - sel1.p * (c.copula.be1.mp - c.copula.be1.mm) * derp1.dereta1.m
der2h.derc12.derc22.pp  <- c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p)
der2h.derc12.derc22.pm  <- c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m)
der2h.derc12.derc22.mp  <- c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p)
der2h.derc12.derc22.mm  <- c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m)

der2h.derc22.derc22.2pm <- c.copula2.be2be2.2pm * (derp2.dereta2.p)^2 + (c.copula.be2.pp - c.copula.be2.mp) * der2p2.dereta2eta2.p
der2h.derc22.derc22.p2m <- c.copula2.be2be2.p2m * (derp2.dereta2.m)^2 + (c.copula.be2.pm - c.copula.be2.mm) * der2p2.dereta2eta2.m
der2h.derc22.derc22.mag <- sel2.m * (c.copula.be2.pp - c.copula.be2.mp) * derp2.dereta2.p - sel2.p * (c.copula.be2.pm - c.copula.be2.mm) * derp2.dereta2.m

der2h.derc12.derc12.2pm_app <- approx.CLM(der2h.derc12.derc12.2pm, h.2p2m, epsilon)
der2h.derc12.derc12.p2m_app <- approx.CLM(der2h.derc12.derc12.p2m, h.2p2m, epsilon)
der2h.derc12.derc12.mag_app <- approx.CLM(der2h.derc12.derc12.mag, h.2p2m, epsilon)
der2h.derc12.derc22.pp_app  <- approx.CLM(der2h.derc12.derc22.pp , h.2p2m, epsilon)
der2h.derc12.derc22.pm_app  <- approx.CLM(der2h.derc12.derc22.pm , h.2p2m, epsilon)
der2h.derc12.derc22.mp_app  <- approx.CLM(der2h.derc12.derc22.mp , h.2p2m, epsilon)
der2h.derc12.derc22.mm_app  <- approx.CLM(der2h.derc12.derc22.mm , h.2p2m, epsilon)

der2h.derc22.derc22.2pm_app <- approx.CLM(der2h.derc22.derc22.2pm, h.2p2m, epsilon)
der2h.derc22.derc22.p2m_app <- approx.CLM(der2h.derc22.derc22.p2m, h.2p2m, epsilon)
der2h.derc22.derc22.mag_app <- approx.CLM(der2h.derc22.derc22.mag, h.2p2m, epsilon)


der2h.derc11.derc11.2p2m  <- D11         * (c.copula2.be1be1.2pm * (derp1.dereta1.p)^2 + (c.copula.be1.pp - c.copula.be1.pm) * der2p1.dereta1eta1.p) - 
                             D12         * (c.copula2.be1be1.p2m * (derp1.dereta1.m)^2 + (c.copula.be1.mp - c.copula.be1.mm) * der2p1.dereta1eta1.m)
der2h.derc11.derc12.2p2m  <- D11.1       * (c.copula2.be1be1.2pm * (derp1.dereta1.p)^2 + (c.copula.be1.pp - c.copula.be1.pm) * der2p1.dereta1eta1.p) - 
                             D12.1       * (c.copula2.be1be1.p2m * (derp1.dereta1.m)^2 + (c.copula.be1.mp - c.copula.be1.mm) * der2p1.dereta1eta1.m)
der2h.derc11.derc21.2p2m  <- D11 * D21   * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p)) - 
                             D11 * D22   * (c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m)) -
                             D12 * D21   * (c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p)) + 
                             D12 * D22   * (c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc11.derc22.2p2m  <- D11 * D21.1 * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p)) - 
                             D11 * D22.1 * (c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m)) -
                             D12 * D21.1 * (c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p)) + 
                             D12 * D22.1 * (c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc11.dereta2.2p2m <- D11         * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p) - c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m)) -
                             D12         * (c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p) - c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc11.derteta.2p2m <- D11         * c.copula2.be1th.2pm * derp1.dereta1.p - 
                             D12         * c.copula2.be1th.p2m * derp1.dereta1.m

der2h.derc12.derc12.2p2m  <- crossprod(D11.1 * (VC$weights * der2h.derc12.derc12.2pm_app / h.2p2m), D11.1) - 
                             crossprod(D12.1 * (VC$weights * der2h.derc12.derc12.p2m_app / h.2p2m), D12.1) +
                             diag(colSums(2 * (VC$weights * der2h.derc12.derc12.mag_app / h.2p2m)))
der2h.derc12.derc21.2p2m  <- D11.1 * D21 * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p)) - 
                             D11.1 * D22 * (c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m)) -
                             D12.1 * D21 * (c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p)) + 
                             D12.1 * D22 * (c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc12.derc22.2p2m  <- crossprod(D11.1 * (VC$weights * der2h.derc12.derc22.pp_app / h.2p2m), D21.1) - 
                             crossprod(D11.1 * (VC$weights * der2h.derc12.derc22.pm_app / h.2p2m), D22.1) -
                             crossprod(D12.1 * (VC$weights * der2h.derc12.derc22.mp_app / h.2p2m), D21.1) + 
                             crossprod(D12.1 * (VC$weights * der2h.derc12.derc22.mm_app / h.2p2m), D22.1)
der2h.derc12.dereta1.2p2m <- D11.1 * (c.copula2.be1be1.2pm * (derp1.dereta1.p)^2 + (c.copula.be1.pp - c.copula.be1.pm) * der2p1.dereta1eta1.p) - 
                             D12.1 * (c.copula2.be1be1.p2m * (derp1.dereta1.m)^2 + (c.copula.be1.mp - c.copula.be1.mm) * der2p1.dereta1eta1.m)
der2h.derc12.dereta2.2p2m <- D11.1 * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p) - c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m)) -
                             D12.1 * (c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p) - c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc12.derteta.2p2m <- D11.1 * c.copula2.be1th.2pm * derp1.dereta1.p - 
                             D12.1 * c.copula2.be1th.p2m * derp1.dereta1.m

der2h.derc21.derc21.2p2m  <- D21   * (c.copula2.be2be2.2pm * (derp2.dereta2.p)^2 + (c.copula.be2.pp - c.copula.be2.mp) * der2p2.dereta2eta2.p) - 
                             D22   * (c.copula2.be2be2.p2m * (derp2.dereta2.m)^2 + (c.copula.be2.pm - c.copula.be2.mm) * der2p2.dereta2eta2.m)
der2h.derc21.derc22.2p2m  <- D21.1 * (c.copula2.be2be2.2pm * (derp2.dereta2.p)^2 + (c.copula.be2.pp - c.copula.be2.mp) * der2p2.dereta2eta2.p) - 
                             D22.1 * (c.copula2.be2be2.p2m * (derp2.dereta2.m)^2 + (c.copula.be2.pm - c.copula.be2.mm) * der2p2.dereta2eta2.m)
der2h.derc21.dereta1.2p2m <- D21   * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p) - c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p)) -
                             D22   * (c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m) - c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc21.derteta.2p2m <- D21   * c.copula2.be2th.2pm * derp2.dereta2.p - 
                             D22   * c.copula2.be2th.p2m * derp2.dereta2.m

der2h.derc22.derc22.2p2m  <- crossprod(D21.1 * (VC$weights * der2h.derc22.derc22.2pm_app / h.2p2m), D21.1) -
                             crossprod(D22.1 * (VC$weights * der2h.derc22.derc22.p2m_app / h.2p2m), D22.1) +
                             diag(colSums(2 * (VC$weights * der2h.derc22.derc22.mag_app / h.2p2m)))
der2h.derc22.dereta1.2p2m <- D21.1 * (c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p) - c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p)) -
                             D22.1 * (c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m) - c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m))
der2h.derc22.dereta2.2p2m <- D21.1 * (c.copula2.be2be2.2pm * (derp2.dereta2.p)^2 + (c.copula.be2.pp - c.copula.be2.mp) * der2p2.dereta2eta2.p) - 
                             D22.1 * (c.copula2.be2be2.p2m * (derp2.dereta2.m)^2 + (c.copula.be2.pm - c.copula.be2.mm) * der2p2.dereta2eta2.m)
der2h.derc22.derteta.2p2m <- D21.1 * c.copula2.be2th.2pm * derp2.dereta2.p - 
                             D22.1 * c.copula2.be2th.p2m * derp2.dereta2.m

der2h.dereta1.dereta1.2p2m    <- c.copula2.be1be1.2pm * (derp1.dereta1.p)^2 + (c.copula.be1.pp - c.copula.be1.pm) * der2p1.dereta1eta1.p - 
                                 c.copula2.be1be1.p2m * (derp1.dereta1.m)^2 - (c.copula.be1.mp - c.copula.be1.mm) * der2p1.dereta1eta1.m
der2h.dereta1.dereta2.2p2m    <- c.copula2.be1be2.pp * (derp1.dereta1.p * derp2.dereta2.p) - c.copula2.be1be2.pm * (derp1.dereta1.p * derp2.dereta2.m) -
                                 c.copula2.be1be2.mp * (derp1.dereta1.m * derp2.dereta2.p) + c.copula2.be1be2.mm * (derp1.dereta1.m * derp2.dereta2.m)
der2h.dereta1.derteta.st.2p2m <- c.copula2.be1th.2pm * derp1.dereta1.p - c.copula2.be1th.p2m * derp1.dereta1.m

der2h.dereta2.dereta2.2p2m    <- c.copula2.be2be2.2pm * (derp2.dereta2.p)^2 + (c.copula.be2.pp - c.copula.be2.mp) * der2p2.dereta2eta2.p - 
                                 c.copula2.be2be2.p2m * (derp2.dereta2.m)^2 - (c.copula.be2.pm - c.copula.be2.mm) * der2p2.dereta2eta2.m
der2h.dereta2.derteta.st.2p2m <- c.copula2.be2th.2pm * derp2.dereta2.p - c.copula2.be2th.p2m * derp2.dereta2.m

der2h.derteta.st2.2p2m <- c.copula2.th2.2p2m * (derteta.derteta.st)^2 + (c.copula.th.2p2m / derteta.derteta.st) * der2teta.derteta.stteta.st   


# Approximating quantities for the Hessian 

der2h.derc11.derc11.2p2m_app  <- approx.CLM(der2h.derc11.derc11.2p2m , h.2p2m, epsilon)
der2h.derc11.derc12.2p2m_app  <- approx.CLM(der2h.derc11.derc12.2p2m , h.2p2m, epsilon)
der2h.derc11.derc21.2p2m_app  <- approx.CLM(der2h.derc11.derc21.2p2m , h.2p2m, epsilon)
der2h.derc11.derc22.2p2m_app  <- approx.CLM(der2h.derc11.derc22.2p2m , h.2p2m, epsilon)
der2h.derc11.dereta2.2p2m_app <- approx.CLM(der2h.derc11.dereta2.2p2m, h.2p2m, epsilon)
der2h.derc11.derteta.2p2m_app <- approx.CLM(der2h.derc11.derteta.2p2m, h.2p2m, epsilon)

der2h.derc12.derc12.2p2m_app  <- der2h.derc12.derc12.2p2m
der2h.derc12.derc21.2p2m_app  <- approx.CLM(der2h.derc12.derc21.2p2m , h.2p2m, epsilon)
der2h.derc12.dereta1.2p2m_app <- approx.CLM(der2h.derc12.dereta1.2p2m, h.2p2m, epsilon)
der2h.derc12.dereta2.2p2m_app <- approx.CLM(der2h.derc12.dereta2.2p2m, h.2p2m, epsilon)
der2h.derc12.derteta.2p2m_app <- approx.CLM(der2h.derc12.derteta.2p2m, h.2p2m, epsilon)

der2h.derc21.derc21.2p2m_app  <- approx.CLM(der2h.derc21.derc21.2p2m , h.2p2m, epsilon)
der2h.derc12.derc22.2p2m_app  <- der2h.derc12.derc22.2p2m
der2h.derc21.derc22.2p2m_app  <- approx.CLM(der2h.derc21.derc22.2p2m , h.2p2m, epsilon)
der2h.derc21.dereta1.2p2m_app <- approx.CLM(der2h.derc21.dereta1.2p2m, h.2p2m, epsilon)
der2h.derc21.derteta.2p2m_app <- approx.CLM(der2h.derc21.derteta.2p2m, h.2p2m, epsilon)

der2h.derc22.derc22.2p2m_app  <- der2h.derc22.derc22.2p2m
der2h.derc22.dereta1.2p2m_app <- approx.CLM(der2h.derc22.dereta1.2p2m, h.2p2m, epsilon)
der2h.derc22.dereta2.2p2m_app <- approx.CLM(der2h.derc22.dereta2.2p2m, h.2p2m, epsilon)
der2h.derc22.derteta.2p2m_app <- approx.CLM(der2h.derc22.derteta.2p2m, h.2p2m, epsilon)

der2h.dereta1.dereta1.2p2m_app    <- approx.CLM(der2h.dereta1.dereta1.2p2m   , h.2p2m, epsilon)
der2h.dereta1.dereta2.2p2m_app    <- approx.CLM(der2h.dereta1.dereta2.2p2m   , h.2p2m, epsilon)
der2h.dereta1.derteta.st.2p2m_app <- approx.CLM(der2h.dereta1.derteta.st.2p2m, h.2p2m, epsilon)

der2h.dereta2.dereta2.2p2m_app    <- approx.CLM(der2h.dereta2.dereta2.2p2m   , h.2p2m, epsilon)
der2h.dereta2.derteta.st.2p2m_app <- approx.CLM(der2h.dereta2.derteta.st.2p2m, h.2p2m, epsilon)

der2h.derteta.st2.2p2m_app <- approx.CLM(der2h.derteta.st2.2p2m, h.2p2m, epsilon)


# Components of the Hessian

d2l.c1_1.c1_1    <- der2h.derc11.derc11.2p2m_app  / h.2p2m - dl.dc1_1^2
d2l.c1_1.c1_2    <- der2h.derc11.derc12.2p2m_app  / h.2p2m - dl.dc1_1 * dl.dc1_2
d2l.c1_1.c2_1    <- der2h.derc11.derc21.2p2m_app  / h.2p2m - dl.dc1_1 * dl.dc2_1
d2l.c1_1.c2_2    <- der2h.derc11.derc22.2p2m_app  / h.2p2m - dl.dc1_1 * dl.dc2_2
d2l.c1_1.be1     <- der2h.derc11.derc11.2p2m_app  / h.2p2m - dl.dc1_1 * dl.dbe1
d2l.c1_1.be2     <- der2h.derc11.dereta2.2p2m_app / h.2p2m - dl.dc1_1 * dl.dbe2
d2l.c1_1.teta.st <- der2h.derc11.derteta.2p2m_app / h.2p2m - dl.dc1_1 * dl.dteta.st

d2l.c1_2.c1_2    <- der2h.derc12.derc12.2p2m_app - crossprod(dl.dc1_2, dl.dc1_2)
d2l.c1_2.c2_1    <- der2h.derc12.derc21.2p2m_app  / h.2p2m - dl.dc1_2 * dl.dc2_1
d2l.c1_2.c2_2    <- der2h.derc12.derc22.2p2m_app - crossprod(dl.dc1_2, dl.dc2_2)
d2l.c1_2.be1     <- der2h.derc12.dereta1.2p2m_app / h.2p2m - dl.dc1_2 * dl.dbe1
d2l.c1_2.be2     <- der2h.derc12.dereta2.2p2m_app / h.2p2m - dl.dc1_2 * dl.dbe2
d2l.c1_2.teta.st <- der2h.derc12.derteta.2p2m_app / h.2p2m - dl.dc1_2 * dl.dteta.st

d2l.c2_1.c2_1    <- der2h.derc21.derc21.2p2m_app  / h.2p2m - dl.dc2_1^2
d2l.c2_1.c2_2    <- der2h.derc21.derc22.2p2m_app  / h.2p2m - dl.dc2_1 * dl.dc2_2
d2l.c2_1.be1     <- der2h.derc21.dereta1.2p2m_app / h.2p2m - dl.dc2_1 * dl.dbe1
d2l.c2_1.be2     <- der2h.derc21.derc21.2p2m_app  / h.2p2m - dl.dc2_1 * dl.dbe2
d2l.c2_1.teta.st <- der2h.derc21.derteta.2p2m_app / h.2p2m - dl.dc2_1 * dl.dteta.st

d2l.c2_2.c2_2    <- der2h.derc22.derc22.2p2m_app - crossprod(dl.dc2_2, dl.dc2_2)
d2l.c2_2.be1     <- der2h.derc22.dereta1.2p2m_app / h.2p2m - dl.dc2_2 * dl.dbe1
d2l.c2_2.be2     <- der2h.derc22.dereta2.2p2m_app / h.2p2m - dl.dc2_2 * dl.dbe2
d2l.c2_2.teta.st <- der2h.derc22.derteta.2p2m_app / h.2p2m - dl.dc2_2 * dl.dteta.st

d2l.be1.be1     <- der2h.dereta1.dereta1.2p2m_app    / h.2p2m - dl.dbe1^2
d2l.be1.be2     <- der2h.dereta1.dereta2.2p2m_app    / h.2p2m - dl.dbe1 * dl.dbe2
d2l.be1.teta.st <- der2h.dereta1.derteta.st.2p2m_app / h.2p2m - dl.dbe1 * dl.dteta.st

d2l.be2.be2     <- der2h.dereta2.dereta2.2p2m_app    / h.2p2m - dl.dbe2^2
d2l.be2.teta.st <- der2h.dereta2.derteta.st.2p2m_app / h.2p2m - dl.dbe2 * dl.dteta.st

d2l.teta.st.teta.st <- der2h.derteta.st2.2p2m_app / h.2p2m - dl.dteta.st^2


# Adding the weigths

d2l.c1_1.c1_1_w    <- VC$weights * d2l.c1_1.c1_1
d2l.c1_1.c1_2_w    <- VC$weights * d2l.c1_1.c1_2
d2l.c1_1.c2_1_w    <- VC$weights * d2l.c1_1.c2_1
d2l.c1_1.c2_2_w    <- VC$weights * d2l.c1_1.c2_2
d2l.c1_1.be1_w     <- VC$weights * d2l.c1_1.be1
d2l.c1_1.be2_w     <- VC$weights * d2l.c1_1.be2
d2l.c1_1.teta.st_w <- VC$weights * d2l.c1_1.teta.st

d2l.c1_2.c1_2_w    <- d2l.c1_2.c1_2 # This quantity has already the weights from its construction 
d2l.c1_2.c2_1_w    <- VC$weights * d2l.c1_2.c2_1
d2l.c1_2.c2_2_w    <- d2l.c1_2.c2_2 # This quantity has already the weights from its construction 
d2l.c1_2.be1_w     <- VC$weights * d2l.c1_2.be1
d2l.c1_2.be2_w     <- VC$weights * d2l.c1_2.be2
d2l.c1_2.teta.st_w <- VC$weights * d2l.c1_2.teta.st

d2l.c2_1.c2_1_w    <- VC$weights * d2l.c2_1.c2_1
d2l.c2_1.c2_2_w    <- VC$weights * d2l.c2_1.c2_2
d2l.c2_1.be1_w     <- VC$weights * d2l.c2_1.be1
d2l.c2_1.be2_w     <- VC$weights * d2l.c2_1.be2
d2l.c2_1.teta.st_w <- VC$weights * d2l.c2_1.teta.st

d2l.c2_2.c2_2_w    <- d2l.c2_2.c2_2 # This quantity has already the weights from its construction 
d2l.c2_2.be1_w     <- VC$weights * d2l.c2_2.be1
d2l.c2_2.be2_w     <- VC$weights * d2l.c2_2.be2
d2l.c2_2.teta.st_w <- VC$weights * d2l.c2_2.teta.st

d2l.be1.be1_w     <- VC$weights * d2l.be1.be1
d2l.be1.be2_w     <- VC$weights * d2l.be1.be2
d2l.be1.teta.st_w <- VC$weights * d2l.be1.teta.st

d2l.be2.be2_w     <- VC$weights * d2l.be2.be2
d2l.be2.teta.st_w <- VC$weights * d2l.be2.teta.st

d2l.teta.st.teta.st_w <- VC$weights * d2l.teta.st.teta.st


##### Global score and Hessian matrix #####

rownames(d2l.c1_2.c1_2_w) <-
rownames(d2l.c1_2.c2_2_w) <-
rownames(d2l.c2_2.c2_2_w) <- NULL

colnames(dl.dc1_2_w        ) <- 
colnames(d2l.c1_1.c1_2_w   ) <- 
colnames(d2l.c1_1.c2_2_w   ) <- 
colnames(d2l.c1_2.c1_2_w   ) <-
colnames(d2l.c1_2.c2_1_w   ) <-
colnames(d2l.c1_2.c2_2_w   ) <-
colnames(d2l.c1_2.be1_w    ) <- 
colnames(d2l.c1_2.be2_w    ) <- 
colnames(d2l.c1_2.teta.st_w) <- 
colnames(d2l.c2_1.c2_2_w   ) <- 
colnames(dl.dc2_2_w        ) <- 
colnames(d2l.c2_2.c2_2_w   ) <- 
colnames(d2l.c2_2.be1_w    ) <-
colnames(d2l.c2_2.be2_w    ) <-
colnames(d2l.c2_2.teta.st_w) <- NULL


if( is.null(VC$X3) ) {

c11.c11  <- sum(d2l.c1_1.c1_1_w)
c11.c12  <- t(rowSums(t(d2l.c1_1.c1_2_w)))
c11.c21  <- sum(d2l.c1_1.c2_1_w)
c11.c22  <- t(rowSums(t(d2l.c1_1.c2_2_w)))
c11.be1  <- t(rowSums(t(-VC$X1 * c(d2l.c1_1.be1_w)))) 
c11.be2  <- t(rowSums(t(-VC$X2 * c(d2l.c1_1.be2_w))))
c11.teta <- sum(d2l.c1_1.teta.st_w)

c12.c12  <- d2l.c1_2.c1_2_w
c12.c21  <- t(t(rowSums(t(d2l.c1_2.c2_1_w))))
c12.c22  <- d2l.c1_2.c2_2_w
c12.be1  <- crossprod(d2l.c1_2.be1_w, -VC$X1)
c12.be2  <- crossprod(d2l.c1_2.be2_w, -VC$X2)
c12.teta <- t(t(rowSums(t(d2l.c1_2.teta.st_w))))

c21.c21  <- sum(d2l.c2_1.c2_1_w)
c21.c22  <- t(rowSums(t(d2l.c2_1.c2_2_w)))
c21.be1  <- t(rowSums(t(-VC$X1 * c(d2l.c2_1.be1_w)))) 
c21.be2  <- t(rowSums(t(-VC$X2 * c(d2l.c2_1.be2_w))))
c21.teta <- sum(d2l.c2_1.teta.st_w)

c22.c22  <- d2l.c2_2.c2_2_w
c22.be1  <- crossprod(d2l.c2_2.be1_w, -VC$X1)
c22.be2  <- crossprod(d2l.c2_2.be2_w, -VC$X2)
c22.teta <- t(t(rowSums(t(d2l.c2_2.teta.st_w))))

be1.be1  <- crossprod(VC$X1 * c(d2l.be1.be1_w), VC$X1)
be1.be2  <- crossprod(VC$X1 * c(d2l.be1.be2_w), VC$X2)
be1.teta <- t(t(rowSums(t(-VC$X1 * c(d2l.be1.teta.st_w)))))

be2.be2  <- crossprod(VC$X2 * c(d2l.be2.be2_w), VC$X2)
be2.teta <- t(t(rowSums(t(-VC$X2 * c(d2l.be2.teta.st_w)))))

teta.teta <- sum(d2l.teta.st.teta.st_w)

G <- - c(sum(      dl.dc1_1_w),
         colSums(  dl.dc1_2_w),
         sum(      dl.dc2_1_w),
         colSums(  dl.dc2_2_w),
         colSums(c(dl.dbe1_w) * -VC$X1),
         colSums(c(dl.dbe2_w) * -VC$X2),
         sum(      dl.dteta.st_w))

H <- - rbind(
cbind(  c11.c11  ,   c11.c12  ,   c11.c21  ,   c11.c22  ,   c11.be1  ,   c11.be2  , c11.teta),
cbind(t(c11.c12 ),   c12.c12  ,   c12.c21  ,   c12.c22  ,   c12.be1  ,   c12.be2  , c12.teta),
cbind(t(c11.c21 ), t(c12.c21 ),   c21.c21  ,   c21.c22  ,   c21.be1  ,   c21.be2  , c21.teta),
cbind(t(c11.c22 ), t(c12.c22 ), t(c21.c22 ),   c22.c22  ,   c22.be1  ,   c22.be2  , c22.teta),
cbind(t(c11.be1 ), t(c12.be1 ), t(c21.be1 ), t(c22.be1) ,   be1.be1  ,   be1.be2  , be1.teta),
cbind(t(c11.be2 ), t(c12.be2 ), t(c21.be2 ), t(c22.be2) , t(be1.be2 ),   be2.be2  , be2.teta),
cbind(t(c11.teta), t(c12.teta), t(c21.teta), t(c22.teta), t(be1.teta), t(be2.teta), teta.teta) )

}

if( !(is.null(VC$X3)) ){

c11.c11  <- sum(d2l.c1_1.c1_1_w)
c11.c12  <- t(rowSums(t(d2l.c1_1.c1_2_w)))
c11.c21  <- sum(d2l.c1_1.c2_1_w)
c11.c22  <- t(rowSums(t(d2l.c1_1.c2_2_w)))
c11.be1  <- t(rowSums(t(-VC$X1 * c(d2l.c1_1.be1_w)))) 
c11.be2  <- t(rowSums(t(-VC$X2 * c(d2l.c1_1.be2_w))))
c11.teta <- t(rowSums(t( VC$X3 * c(d2l.c1_1.teta.st_w))))

c12.c12  <- d2l.c1_2.c1_2_w
c12.c21  <- t(t(rowSums(t(d2l.c1_2.c2_1_w))))
c12.c22  <- d2l.c1_2.c2_2_w
c12.be1  <- crossprod(d2l.c1_2.be1_w,    -VC$X1)
c12.be2  <- crossprod(d2l.c1_2.be2_w,    -VC$X2)
c12.teta <- crossprod(d2l.c1_2.teta.st_w, VC$X3)

c21.c21  <- sum(d2l.c2_1.c2_1_w)
c21.c22  <- t(rowSums(t(d2l.c2_1.c2_2_w)))
c21.be1  <- t(rowSums(t(-VC$X1 * c(d2l.c2_1.be1_w)))) 
c21.be2  <- t(rowSums(t(-VC$X2 * c(d2l.c2_1.be2_w))))
c21.teta <- t(rowSums(t( VC$X3 * c(d2l.c2_1.teta.st_w))))

c22.c22  <- d2l.c2_2.c2_2_w
c22.be1  <- crossprod(d2l.c2_2.be1_w,     -VC$X1)
c22.be2  <- crossprod(d2l.c2_2.be2_w,     -VC$X2)
c22.teta <- crossprod(d2l.c2_2.teta.st_w,  VC$X3)

be1.be1  <- crossprod( VC$X1 * c(d2l.be1.be1_w),     VC$X1)
be1.be2  <- crossprod( VC$X1 * c(d2l.be1.be2_w),     VC$X2)
be1.teta <- crossprod(-VC$X1 * c(d2l.be1.teta.st_w), VC$X3)

be2.be2  <- crossprod( VC$X2 * c(d2l.be2.be2_w),     VC$X2)
be2.teta <- crossprod(-VC$X2 * c(d2l.be2.teta.st_w), VC$X3)

teta.teta <- crossprod(VC$X3 * c(d2l.teta.st.teta.st_w), VC$X3)

G <- - c(sum(      dl.dc1_1_w),
         colSums(  dl.dc1_2_w),
         sum(      dl.dc2_1_w),
         colSums(  dl.dc2_2_w),
         colSums(c(dl.dbe1_w)     * -VC$X1),
         colSums(c(dl.dbe2_w)     * -VC$X2),
         colSums(c(dl.dteta.st_w) *  VC$X3))

H <- - rbind(
cbind(  c11.c11  ,   c11.c12  ,   c11.c21  ,   c11.c22  ,   c11.be1  ,   c11.be2  , c11.teta),
cbind(t(c11.c12 ),   c12.c12  ,   c12.c21  ,   c12.c22  ,   c12.be1  ,   c12.be2  , c12.teta),
cbind(t(c11.c21 ), t(c12.c21 ),   c21.c21  ,   c21.c22  ,   c21.be1  ,   c21.be2  , c21.teta),
cbind(t(c11.c22 ), t(c12.c22 ), t(c21.c22 ),   c22.c22  ,   c22.be1  ,   c22.be2  , c22.teta),
cbind(t(c11.be1 ), t(c12.be1 ), t(c21.be1 ), t(c22.be1) ,   be1.be1  ,   be1.be2  , be1.teta),
cbind(t(c11.be2 ), t(c12.be2 ), t(c21.be2 ), t(c22.be2) , t(be1.be2 ),   be2.be2  , be2.teta),
cbind(t(c11.teta), t(c12.teta), t(c21.teta), t(c22.teta), t(be1.teta), t(be2.teta), teta.teta) )

}

##############################

res <- - sum(l.par)


########################################################################################## 

if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
S.h  <- ps$S.h  


##############################

# Independence model

if (VC$ind.ord == "TRUE") {
	G <- G[-VC$drop.ind]
	H <- H[-VC$drop.ind, -VC$drop.ind]

	if(length(S.h) != 1) S.h <- S.h[-VC$drop.ind, -VC$drop.ind]
}

##############################


if(length(S.h) != 1){
  
S.h1 <- 0.5 * crossprod(params, S.h) %*% params
S.h2 <- S.h %*% params
  
} else S.h <- S.h1 <- S.h2 <- 0   
 

S.res <- res
res   <- S.res + S.h1
G     <- G + S.h2
H     <- H + S.h  


if(VC$extra.regI == "sED") H <- regH(H, type = 2)


list(value = res, gradient = G, hessian = H, S.h = S.h, S.h1 = S.h1, S.h2 = S.h2, l = S.res, l.ln = l.ln, l.par = l.par, ps = ps,
     eta1 = eta1, eta2 = eta2, etad = etad, lp1 = lp1, lp2 = lp2,
     dl.dbe1 = dl.dbe1, dl.dbe2 = dl.dbe2, dl.dteta.st = dl.dteta.st,
     BivD = VC$BivD, p1.p = p1.p, p1.m = p1.m, p2 = p2, theta.star = teta.st,
     teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
     Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)

}