bCopulaCLMgHsCont <- function(params, respvec, VC, ps, AT = FALSE){

n <- VC$n
y1 <- respvec$y1

# With the cut points we might need an extra parameter to be counted in the dimensions below: K1 will be the 
# dimension of the cut-points

c1.ti <- params[1:(VC$K1-1)]

eta1 <- VC$X1%*%params[VC$K1:(VC$K1 + VC$X1.d2 - 1)]
eta2 <- VC$X2%*%params[(VC$K1 + VC$X1.d2):(VC$K1 + VC$X1.d2 + VC$X2.d2 - 1)]
etad <- etas <- l.ln <- NULL

epsilon <-  0.0000001 # RECALL: results are sensible to the choice of epsilon (try also 0.00001)
#max.p   <- 0.9999999

if ( is.null(VC$X3) ) {
sigma2.st <- etas <- params[VC$K1 + VC$X1.d2 + VC$X2.d2    ]
teta.st   <- etad <- params[VC$K1 + VC$X1.d2 + VC$X2.d2 + 1]
}

if ( !(is.null(VC$X3)) ) {
sigma2.st <- etas <- VC$X3%*%params[(VC$K1 + VC$X1.d2 + VC$X2.d2):(VC$K1 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 - 1)]
teta.st   <- etad <- VC$X4%*%params[(VC$K1 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2):(VC$K1 + VC$X1.d2 + VC$X2.d2 + VC$X3.d2 + VC$X4.d2 - 1)]
}


##############################

# Independence model

if (VC$ind.ord == "TRUE") teta.st <- etad <- log(epsilon) # This ensures that teta = 1 (i.e. independence for J0 copula)

##############################


# Cut points are transformed and the linear predictors created 

infty  <- 1e+25 ; e <- .Machine$double.eps * 10^6

c1  <- VC$c1 ; c1[1] <- c1.ti[1] ; for (i in 2 : (VC$K1 - 1)) c1[i] <- c1[i - 1] + c1.ti[i]^2
c1  <- sign(c1) * pmin(10000 * infty, abs(c1))

c1.m   <- t(matrix(nrow = VC$K1 - 1, ncol = n        , c1  )) 
eta1.m <-   matrix(nrow = n        , ncol = VC$K1 - 1, eta1)

lp1 <- c1.m - eta1.m ; lp.p <- cbind(lp1, infty) ; lp.m <- cbind(-infty, lp1)

sel <- VC$sel ; lp.p.sel <- rowSums(lp.p * sel) ; lp.m.sel <- rowSums(lp.m * sel)

lp.p <- lp.p.sel ; lp.m <- lp.m.sel

##############################

####

sstr1 <- esp.tr(sigma2.st, VC$margins[2])
sigma2.st <- sstr1$vrb.st
sigma2    <- sstr1$vrb

####

dHs <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2 = VC$margins[2], naive = FALSE) # WHY VC$margins[2]?

pdf2                         <- dHs$pdf2
p2                           <- dHs$p2
derpdf2.dereta2              <- dHs$derpdf2.dereta2
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st
derp2.dersigma.st            <- dHs$derp2.dersigma.st
derp2.dereta2                <- dHs$derp2.dereta2
der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st
der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2

####

# Definire qua tutte le quantita' che mi servono per lo score and hessian

####

pd1.p <- probm(lp.p, VC$margins[1], only.pr = FALSE, bc = TRUE, CLM = TRUE) # I set only.pr = FALSE to obtain d.n 
pd1.m <- probm(lp.m, VC$margins[1], only.pr = FALSE, bc = TRUE, CLM = TRUE) #

p1.p <- pd1.p$pr  # pd1.p$pr or 1 - pd1.p$pr??
p1.m <- pd1.m$pr

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

if(VC$BivD %in% VC$BivD2[1 : 4 ]) teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov)    , TRUE, FALSE)
if(VC$BivD %in% VC$BivD2[5 : 12]) teta.ind1 <- ifelse(VC$my.env$signind * teta > exp(VC$zerov) + 1, TRUE, FALSE) 
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

h.p <- h.m <- NA

if(length(teta1) != 0){

dH1.p <- copgHs(p1.p[teta.ind1], p2[teta.ind1], eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, CLM = TRUE)
dH1.m <- copgHs(p1.m[teta.ind1], p2[teta.ind1], eta1 = NULL, eta2 = NULL, teta1, teta.st1, Cop1, VC$dof, CLM = TRUE)

h.p[teta.ind1]  <- dH1.p$c.copula.be2 ; h.m[teta.ind1] <- dH1.m$c.copula.be2

}

if(length(teta2) != 0){

dH2.p <- copgHs(p1.p[teta.ind2], p2[teta.ind2], eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, CLM = TRUE)
dH2.m <- copgHs(p1.m[teta.ind2], p2[teta.ind2], eta1 = NULL, eta2 = NULL, teta2, teta.st2, Cop2, VC$dof, CLM = TRUE)

h.p[teta.ind2]  <- dH2.p$c.copula.be2 ; h.m[teta.ind2] <- dH2.m$c.copula.be2

}

h.pm <- pmax(epsilon, h.p - h.m)


l.par <- VC$weights * (log(h.pm) + log(pdf2))


##### Score vector #####

c.copula2.be2.p <- c.copula2.be1be2.p <- c.copula2.be2th.p <- 
c.copula2.be2.m <- c.copula2.be1be2.m <- c.copula2.be2th.m <- NA

if(length(teta1) != 0){

c.copula2.be2.p[teta.ind1]    <- dH1.p$c.copula2.be2    ; c.copula2.be2.m[teta.ind1]    <- dH1.m$c.copula2.be2
c.copula2.be1be2.p[teta.ind1] <- dH1.p$c.copula2.be1be2 ; c.copula2.be1be2.m[teta.ind1] <- dH1.m$c.copula2.be1be2
c.copula2.be2th.p[teta.ind1]  <- dH1.p$c.copula2.be2th  ; c.copula2.be2th.m[teta.ind1]  <- dH1.m$c.copula2.be2th

}

if(length(teta2) != 0){

c.copula2.be2.p[teta.ind2]    <- dH2.p$c.copula2.be2    ; c.copula2.be2.m[teta.ind2]    <- dH2.m$c.copula2.be2
c.copula2.be1be2.p[teta.ind2] <- dH2.p$c.copula2.be1be2 ; c.copula2.be1be2.m[teta.ind2] <- dH2.m$c.copula2.be1be2
c.copula2.be2th.p[teta.ind2]  <- dH2.p$c.copula2.be2th  ; c.copula2.be2th.m[teta.ind2]  <- dH2.m$c.copula2.be2th

}


# DO WE NEED TO ADJUST THE APPROX OF THE QUANTITIES BELOW?

c.copula2.be2.pm   <- c.copula2.be2.p   - c.copula2.be2.m   # Do we need any correction of the approx here???
c.copula2.be2th.pm <- c.copula2.be2th.p - c.copula2.be2th.m 


cc1.ti <- matrix(nrow = n, ncol = VC$K1 - 2, byrow = TRUE, c1.ti[2 : (VC$K1 - 1)])


#################
#################

sel.p <- VC$sel.p 
sel.m <- VC$sel.m

sel.mm <- VC$sel.mm
sel.pm <- VC$sel.pm

#################
#################

#if ( VC$K1 > 3 ) {
#for (l in 1 : (VC$K1 - 3)) {sel.pm [, l] <- rowSums(sel.p[, l : (VC$K1 - 2)])} ; sel.p[, 1 : (VC$K1 - 3)] <- sel.pm
#for (l in 1 : (VC$K1 - 3)) {sel.mm [, l] <- rowSums(sel.m[, l : (VC$K1 - 2)])} ; sel.m[, 1 : (VC$K1 - 3)] <- sel.mm
#}


#d21 <- d22 <- matrix(nrow = dim(sel.m)[1], ncol = dim(sel.m)[2], 0)
#for(j in 1 : dim(d21)[2]) {d21[, j] <- rowSums(sel.m) - j + 1 ; d21[, j][d21[, j] < 0] <- rep(0)}
#for(j in 1 : dim(d22)[2]) {d22[, j] <- rowSums(sel.p) - j + 1 ; d22[, j][d22[, j] < 0] <- rep(0)}

#D11 <-  rowSums(sel[, 1 : (VC$K1 - 1)]) ; D12 <- rowSums(sel[, 2 : VC$K1])  # there was no () for 1:K1-1 

D11 <- VC$D11
D12 <- VC$D12


D21 <-  2 * sel.m * cc1.ti
D22 <- 2 * sel.p * cc1.ti


derp1.dereta1.p <- - pd1.p$derp1.dereta1                  ; derp1.dereta1.m <- - pd1.m$derp1.dereta1 # RECALL: The function probm was originally constructed for 1 - pd1.p$pr. I use instead pd1.p$pr
derh.dereta1.p  <-   c.copula2.be1be2.p * derp1.dereta1.p ; derh.dereta1.m  <-   c.copula2.be1be2.m * derp1.dereta1.m

derh.dereta2.pm      <- c.copula2.be2.pm * derp2.dereta2
derh.dersigma2.st.pm <- c.copula2.be2.pm * derp2.dersigma.st


# Approximating quantities for the score 

c.copula2.be2th.pm   <- approx.CLM(c.copula2.be2th.pm  , h.pm, epsilon)
derh.dereta2.pm      <- approx.CLM(derh.dereta2.pm     , h.pm, epsilon)
derh.dersigma2.st.pm <- approx.CLM(derh.dersigma2.st.pm, h.pm, epsilon)

derpdf2.dereta2      <- approx.CLM(derpdf2.dereta2     , pdf2, epsilon)
derpdf2.dersigma2.st <- approx.CLM(derpdf2.dersigma2.st, pdf2, epsilon)


dl.dc1_1 <-   1/(h.pm)  * (D11 *   derh.dereta1.p  - D12 *   derh.dereta1.m)
dl.dc1_2 <- c(1/(h.pm)) * (D21 * c(derh.dereta1.p) - D22 * c(derh.dereta1.m))

dl.dbe1 <- 1/(h.pm) * (derh.dereta1.p - derh.dereta1.m)

dl.dbe2      <- 1/(h.pm) * derh.dereta2.pm      + derpdf2.dereta2      / pdf2
dl.dsigma.st <- 1/(h.pm) * derh.dersigma2.st.pm + derpdf2.dersigma2.st / pdf2

dl.dteta.st <- 1/(h.pm) * c.copula2.be2th.pm


# Removing the approximation of the score

dl.dbe2_1      <- 1/(h.pm) * derh.dereta2.pm
dl.dsigma.st_1 <- 1/(h.pm) * derh.dersigma2.st.pm


c.copula2.be2th.pm  <- c.copula2.be2th.p - c.copula2.be2th.m
derh.dereta2.pm     <- c.copula2.be2.pm * derp2.dereta2
derh.dersigma.st.pm <- derp2.dersigma.st * c.copula2.be2.pm

#derpdf2.dereta2      <- dHs$derpdf2.dereta2
#derpdf2.dersigma2.st <- dHs$derpdf2.dersigma2.st


# Adding the weigths

dl.dc1_1     <- VC$weights * dl.dc1_1
dl.dc1_2     <- VC$weights * dl.dc1_2
dl.dbe1      <- VC$weights * dl.dbe1
dl.dbe2      <- VC$weights * dl.dbe2
dl.dsigma.st <- VC$weights * dl.dsigma.st


##### Hessian matrix #####

der2h.derp1p1.p         <- der2h.derp1p1.m            <-
der2h.derp1p2.p         <- der2h.derp1p2.m            <-
der2h.derp1teta.p       <- der2h.derp1teta.m          <-
der2h.derp2p2.p         <- der2h.derp2p2.m            <-
der2h.derp2teta.p       <- der2h.derp2teta.m          <-
der2h.derteta.teta.st.p <- der2h.derteta.teta.st.m    <-
derteta.derteta.st      <- der2teta.derteta.stteta.st <- NA

der2p1.dereta1eta1.p <- - pd1.p$der2p1.dereta1eta1 # RECALL: The function probm was originally constructed for 1 - pd1.p$pr. I use instead pd1.p$pr
der2p1.dereta1eta1.m <- - pd1.m$der2p1.dereta1eta1


if(length(teta1) != 0){

BITS1.p <- copgHsCont(p1.p[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, par2 = VC$dof, nu.st = log(VC$dof - 2))
BITS1.m <- copgHsCont(p1.m[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, par2 = VC$dof, nu.st = log(VC$dof - 2))

der2h.derp1p1.p[teta.ind1]         <- BITS1.p$der2h.derp1p1         ; der2h.derp1p1.m[teta.ind1]         <- BITS1.m$der2h.derp1p1
der2h.derp1p2.p[teta.ind1]         <- BITS1.p$der2h.derp1p2         ; der2h.derp1p2.m[teta.ind1]         <- BITS1.m$der2h.derp1p2
der2h.derp1teta.p[teta.ind1]       <- BITS1.p$der2h.derp1teta       ; der2h.derp1teta.m[teta.ind1]       <- BITS1.m$der2h.derp1teta
der2h.derp2p2.p[teta.ind1]         <- BITS1.p$der2h.derp2p2         ; der2h.derp2p2.m[teta.ind1]         <- BITS1.m$der2h.derp2p2
der2h.derp2teta.p[teta.ind1]       <- BITS1.p$der2h.derp2teta       ; der2h.derp2teta.m[teta.ind1]       <- BITS1.m$der2h.derp2teta    
der2h.derteta.teta.st.p[teta.ind1] <- BITS1.p$der2h.derteta.teta.st ; der2h.derteta.teta.st.m[teta.ind1] <- BITS1.m$der2h.derteta.teta.st

}

if(length(teta2) != 0){

BITS2.p <- copgHsCont(p1.p[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, par2 = VC$dof, nu.st = log(VC$dof - 2)) #
BITS2.m <- copgHsCont(p1.m[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, par2 = VC$dof, nu.st = log(VC$dof - 2)) # teta2 & teta.st2 or teta1 & teta.st1?

der2h.derp1p1.p[teta.ind2]         <- BITS2.p$der2h.derp1p1         ; der2h.derp1p1.m[teta.ind2]         <- BITS2.m$der2h.derp1p1
der2h.derp1p2.p[teta.ind2]         <- BITS2.p$der2h.derp1p2         ; der2h.derp1p2.m[teta.ind2]         <- BITS2.m$der2h.derp1p2
der2h.derp1teta.p[teta.ind2]       <- BITS2.p$der2h.derp1teta       ; der2h.derp1teta.m[teta.ind2]       <- BITS2.m$der2h.derp1teta
der2h.derp2p2.p[teta.ind2]         <- BITS2.p$der2h.derp2p2         ; der2h.derp2p2.m[teta.ind2]         <- BITS2.m$der2h.derp2p2
der2h.derp2teta.p[teta.ind2]       <- BITS2.p$der2h.derp2teta       ; der2h.derp2teta.m[teta.ind2]       <- BITS2.m$der2h.derp2teta    
der2h.derteta.teta.st.p[teta.ind2] <- BITS2.p$der2h.derteta.teta.st ; der2h.derteta.teta.st.m[teta.ind2] <- BITS2.m$der2h.derteta.teta.st

}


der2h.derp2p2.pm           <- der2h.derp2p2.p - der2h.derp2p2.m
der2h.derp2dersigma2.st.pm <- der2h.derp2p2.pm * derp2.dersigma.st
der2h.derp2teta.pm         <- der2h.derp2teta.p - der2h.derp2teta.m
der2h.derteta.teta.st.pm   <- der2h.derteta.teta.st.p - der2h.derteta.teta.st.m


if(length(teta1) != 0){

derteta.derteta.st[teta.ind1]         <- BITS1.p$derteta.derteta.st # RECALL: this quantity does not depend on the margins, so the use of BITS.p and BITS.m is indifferent 
der2teta.derteta.stteta.st[teta.ind1] <- BITS1.p$der2teta.derteta.stteta.st

}

if(length(teta2) != 0){

derteta.derteta.st[teta.ind2]         <- BITS2.p$derteta.derteta.st # RECALL: this quantity does not depend on the margins, so the use of BITS.p and BITS.m is indifferent 
der2teta.derteta.stteta.st[teta.ind2] <- BITS2.p$der2teta.derteta.stteta.st

}


der2h.dereta1.dereta1.p <- der2h.derp1p1.p * (derp1.dereta1.p^2) + c.copula2.be1be2.p * der2p1.dereta1eta1.p
der2h.dereta1.dereta1.m <- der2h.derp1p1.m * (derp1.dereta1.m^2) + c.copula2.be1be2.m * der2p1.dereta1eta1.m

der2h.dereta1.dereta2.pm      <- (der2h.derp1p2.p   * derp1.dereta1.p - der2h.derp1p2.m   * derp1.dereta1.m) * derp2.dereta2
der2h.dereta1.dersigma2.st.pm <- (der2h.derp1p2.p   * derp1.dereta1.p - der2h.derp1p2.m   * derp1.dereta1.m) * derp2.dersigma.st
der2h.dereta1.derteta.st.pm   <- (der2h.derp1teta.p * derp1.dereta1.p - der2h.derp1teta.m * derp1.dereta1.m) * derteta.derteta.st

der2h.dereta1.derp2.c1.pm   <- D11 * der2h.derp1p2.p   * derp1.dereta1.p - D12 * der2h.derp1p2.m   * derp1.dereta1.m
der2h.dereta1.derteta.c1.pm <- D11 * der2h.derp1teta.p * derp1.dereta1.p - D12 * der2h.derp1teta.m * derp1.dereta1.m

der2h.dereta2.dereta2.pm      <- (der2h.derp2p2.pm  * (derp2.dereta2)^2 + c.copula2.be2.pm * der2p2.dereta2eta2)
der2h.dereta2.dersigma2.st.pm <-  der2h.derp2dersigma2.st.pm * derp2.dereta2 + c.copula2.be2.pm * der2p2.dereta2dersigma2.st
der2h.dereta2.derteta.st.pm   <-  der2h.derp2teta.pm * derp2.dereta2 * derteta.derteta.st

der2h.dersigma2.st2.pm           <- der2h.derp2dersigma2.st.pm * derp2.dersigma.st + c.copula2.be2.pm * der2p2.dersigma2.st2
der2h.derteta.st.dersigma2.st.pm <- der2h.derp2teta.pm * derteta.derteta.st * derp2.dersigma.st

der2h.derteta.st2.pm <- der2h.derteta.teta.st.pm * (derteta.derteta.st)^2 + (c.copula2.be2th.pm/derteta.derteta.st) * der2teta.derteta.stteta.st


# Approximating quantities for the Hessian 

der2h.dereta1.dereta2.pm      <- approx.CLM(der2h.dereta1.dereta2.pm     , h.pm, epsilon)
der2h.dereta1.dersigma2.st.pm <- approx.CLM(der2h.dereta1.dersigma2.st.pm, h.pm, epsilon)
der2h.dereta1.derteta.st.pm   <- approx.CLM(der2h.dereta1.derteta.st.pm  , h.pm, epsilon)

der2h.dereta2.dereta2.pm      <- approx.CLM(der2h.dereta2.dereta2.pm     , h.pm, epsilon)
der2h.dereta2.dersigma2.st.pm <- approx.CLM(der2h.dereta2.dersigma2.st.pm, h.pm, epsilon)
der2h.dereta2.derteta.st.pm   <- approx.CLM(der2h.dereta2.derteta.st.pm  , h.pm, epsilon)

der2h.dersigma2.st2.pm           <- approx.CLM(der2h.dersigma2.st2.pm          , h.pm, epsilon)
der2h.derteta.st.dersigma2.st.pm <- approx.CLM(der2h.derteta.st.dersigma2.st.pm, h.pm, epsilon)

der2h.derteta.st2.pm <- approx.CLM(der2h.derteta.st2.pm, h.pm, epsilon)

der2pdf2.dereta2             <- approx.CLM(der2pdf2.dereta2            , pdf2, epsilon)
der2pdf2.dereta2dersigma2.st <- approx.CLM(der2pdf2.dereta2dersigma2.st, pdf2, epsilon)

# Approximating (composite) quantities for the Hessian 

d2l.c1_1.c1_1.NUM_1     <- approx.CLM(D11 *   der2h.dereta1.dereta1.p  - D12 *   der2h.dereta1.dereta1.m , h.pm, epsilon)
d2l.c1_1.c1_2.NUM_1     <- approx.CLM(D21 * c(der2h.dereta1.dereta1.p) - D22 * c(der2h.dereta1.dereta1.m), h.pm, epsilon) 
d2l.c1_1.be1.NUM_1      <- approx.CLM(D11 *   der2h.dereta1.dereta1.p  - D12 *   der2h.dereta1.dereta1.m , h.pm, epsilon)
d2l.c1_1.be2.NUM_1      <- approx.CLM(        der2h.dereta1.derp2.c1.pm      *   derp2.dereta2           , h.pm, epsilon)
d2l.c1_1.sigma.st.NUM_1 <- approx.CLM(        der2h.dereta1.derp2.c1.pm      *   derp2.dersigma.st       , h.pm, epsilon)
d2l.c1_1.teta.st.NUM_1  <- approx.CLM(        der2h.dereta1.derteta.c1.pm    *   derteta.derteta.st      , h.pm, epsilon)

d2l.c1_2.be1.NUM_1      <- approx.CLM( D21 * c(der2h.dereta1.dereta1.p            ) - D22 * c(der2h.dereta1.dereta1.m            )                         , h.pm, epsilon)
d2l.c1_2.be2.NUM_1      <- approx.CLM((D21 * c(der2h.derp1p2.p   * derp1.dereta1.p) - D22 * c(der2h.derp1p2.m * derp1.dereta1.m  )) * c(derp2.dereta2     ), h.pm, epsilon)
d2l.c1_2.sigma.st.NUM_1 <- approx.CLM((D21 * c(der2h.derp1p2.p   * derp1.dereta1.p) - D22 * c(der2h.derp1p2.m * derp1.dereta1.m  )) * c(derp2.dersigma.st ), h.pm, epsilon)
d2l.c1_2.teta.st.NUM_1  <- approx.CLM((D21 * c(der2h.derp1teta.p * derp1.dereta1.p) - D22 * c(der2h.derp1teta.m * derp1.dereta1.m)) * c(derteta.derteta.st), h.pm, epsilon)

d2l.be1.be1.NUM_1 <- approx.CLM(der2h.dereta1.dereta1.p - der2h.dereta1.dereta1.m, h.pm, epsilon)


d2l.c1_1.c1_1     <-   1/(h.pm)  * d2l.c1_1.c1_1.NUM_1     -   1/(h.pm)^2 * (D11 * derh.dereta1.p - D12 * derh.dereta1.m)^2
d2l.c1_1.c1_2     <- c(1/(h.pm)) * d2l.c1_1.c1_2.NUM_1     - c(1/(h.pm)^2 * (D11 * derh.dereta1.p - D12 * derh.dereta1.m))  * (D21 * c(derh.dereta1.p) - D22 * c(derh.dereta1.m))
d2l.c1_1.be1      <-   1/(h.pm)  * d2l.c1_1.be1.NUM_1      -   1/(h.pm)^2 * (D11 * derh.dereta1.p - D12 * derh.dereta1.m)   * (        derh.dereta1.p  -         derh.dereta1.m )
d2l.c1_1.be2      <-   1/(h.pm)  * d2l.c1_1.be2.NUM_1      -   1/(h.pm)^2 * (D11 * derh.dereta1.p - D12 * derh.dereta1.m)   *  derh.dereta2.pm
d2l.c1_1.sigma.st <-   1/(h.pm)  * d2l.c1_1.sigma.st.NUM_1 -   1/(h.pm)^2 * (D11 * derh.dereta1.p - D12 * derh.dereta1.m)   *  derh.dersigma2.st.pm
d2l.c1_1.teta.st  <-   1/(h.pm)  * d2l.c1_1.teta.st.NUM_1  -   1/(h.pm)^2 * (D11 * derh.dereta1.p - D12 * derh.dereta1.m)   *  c.copula2.be2th.pm

d2l.c1_2.be1      <- c(1/(h.pm)) * d2l.c1_2.be1.NUM_1      - c(1/(h.pm))^2 * (D21 * c(derh.dereta1.p) - D22 * c(derh.dereta1.m)) * c(derh.dereta1.p - derh.dereta1.m     )
d2l.c1_2.be2      <- c(1/(h.pm)) * d2l.c1_2.be2.NUM_1      - c(1/(h.pm))^2 * (D21 * c(derh.dereta1.p) - D22 * c(derh.dereta1.m)) * c(derh.dereta2.pm                     )
d2l.c1_2.sigma.st <- c(1/(h.pm)) * d2l.c1_2.sigma.st.NUM_1 - c(1/(h.pm))^2 * (D21 * c(derh.dereta1.p) - D22 * c(derh.dereta1.m)) * c(c.copula2.be2.pm * derp2.dersigma.st)
d2l.c1_2.teta.st  <- c(1/(h.pm)) * d2l.c1_2.teta.st.NUM_1  - c(1/(h.pm))^2 * (D21 * c(derh.dereta1.p) - D22 * c(derh.dereta1.m)) * c(c.copula2.be2th.pm                  )

d2l.be1.be1      <- 1/(h.pm) * d2l.be1.be1.NUM_1             - dl.dbe1^2
d2l.be1.be2      <- 1/(h.pm) * der2h.dereta1.dereta2.pm      - dl.dbe1   * dl.dbe2_1
d2l.be1.sigma.st <- 1/(h.pm) * der2h.dereta1.dersigma2.st.pm - dl.dbe1   * dl.dsigma.st_1
d2l.be1.teta.st  <- 1/(h.pm) * der2h.dereta1.derteta.st.pm   - dl.dbe1   * dl.dteta.st

d2l.be2.be2      <- 1/(h.pm) * der2h.dereta2.dereta2.pm      - dl.dbe2_1^2 +
                    (der2pdf2.dereta2 / pdf2 - (derpdf2.dereta2 / pdf2)^2)
d2l.be2.sigma.st <- 1/(h.pm) * der2h.dereta2.dersigma2.st.pm - dl.dbe2_1 * dl.dsigma.st_1 +
                    (der2pdf2.dereta2dersigma2.st / pdf2 - (derpdf2.dereta2 / pdf2) * (derpdf2.dersigma2.st / pdf2))
d2l.be2.teta.st  <- 1/(h.pm) * der2h.dereta2.derteta.st.pm   - dl.dbe2_1 * dl.dteta.st

d2l.sigma.st.sigma.st <- 1/(h.pm) * der2h.dersigma2.st2.pm           - dl.dsigma.st_1^2 + 
                         (der2pdf2.dersigma2.st2 / pdf2 - (derpdf2.dersigma2.st / pdf2)^2)
d2l.sigma.st.teta.st  <- 1/(h.pm) * der2h.derteta.st.dersigma2.st.pm - dl.dsigma.st_1 * dl.dteta.st

d2l.teta.st.teta.st <- 1/(h.pm) * der2h.derteta.st2.pm - dl.dteta.st^2


### d2l.c1_2.c1_2

D21.c2c2 <- D22.c2c2 <- list() ; D.c2c2 <- rep(0, dim(D21)[2])
der2c2.derc2c2.p <- der2c2.derc2c2.m <- der2c2.derc2c2.pm <- list()

der2h.derc2.derc2.p <- crossprod(c(1/(h.pm)) * D21 * c(der2h.dereta1.dereta1.p), D21)
der2h.derc2.derc2.m <- crossprod(c(1/(h.pm)) * D22 * c(der2h.dereta1.dereta1.m), D22)

der2h.derc2.derc2.pm <- der2h.derc2.derc2.p - der2h.derc2.derc2.m


y1.mod1 <- y1.mod2 <- y1 ; y1.mod1[which(y1.mod1 == max(y1)    )] <- rep(min(y1))
                           y1.mod2[which(y1.mod2 == min(y1) + 1)] <- rep(min(y1))

for(i in 1 : dim(D21)[1]) {D1.c2c2 <- D2.c2c2 <- D.c2c2 ; D1.w2 <- max(y1.mod1[i] - 1, 0)
                                                          D2.w2 <- max(y1.mod2[i] - 2, 0)

                           if(D1.w2 != 0) D1.c2c2[1 : D1.w2] <- rep(2)
                           if(D2.w2 != 0) D2.c2c2[1 : D2.w2] <- rep(2)

                           D21.c2c2[[i]] <- D1.c2c2
                           D22.c2c2[[i]] <- D2.c2c2}

for(i in 1 : dim(D21)[1]) {der2c2.derc2c2.p [[i]] <- derh.dereta1.p[i] * D21.c2c2[[i]]
                           der2c2.derc2c2.m [[i]] <- derh.dereta1.m[i] * D22.c2c2[[i]]
                           der2c2.derc2c2.pm[[i]] <- 1/(h.pm)[i] * (der2c2.derc2c2.p[[i]] - der2c2.derc2c2.m[[i]])}


#der2c2.derc2c2.pm <- diag(Reduce('+', der2c2.derc2c2.pm))


der2c2.derc2c2.pm <- Reduce('+', der2c2.derc2c2.pm)
if(length(der2c2.derc2c2.pm) > 1) der2c2.derc2c2.pm <- diag(der2c2.derc2c2.pm)




d2l.c1_2.c1_2 <- der2h.derc2.derc2.pm + der2c2.derc2c2.pm -
                 crossprod(dl.dc1_2, dl.dc1_2)


# Adding the weigths

d2l.c1_1.c1_1         <- VC$weights * d2l.c1_1.c1_1
d2l.c1_1.c1_2         <- VC$weights * d2l.c1_1.c1_2
d2l.c1_1.be1          <- VC$weights * d2l.c1_1.be1
d2l.c1_1.be2          <- VC$weights * d2l.c1_1.be2
d2l.c1_1.sigma.st     <- VC$weights * d2l.c1_1.sigma.st
d2l.c1_1.teta.st      <- VC$weights * d2l.c1_1.teta.st

d2l.c1_2.be1          <- VC$weights * d2l.c1_2.be1
d2l.c1_2.be2          <- VC$weights * d2l.c1_2.be2
d2l.c1_2.sigma.st     <- VC$weights * d2l.c1_2.sigma.st
d2l.c1_2.teta.st      <- VC$weights * d2l.c1_2.teta.st

d2l.be1.be1           <- VC$weights * d2l.be1.be1
d2l.be1.be2           <- VC$weights * d2l.be1.be2
d2l.be1.sigma.st      <- VC$weights * d2l.be1.sigma.st
d2l.be1.teta.st       <- VC$weights * d2l.be1.teta.st

d2l.be2.be2           <- VC$weights * d2l.be2.be2
d2l.be2.sigma.st      <- VC$weights * d2l.be2.sigma.st
d2l.be2.teta.st       <- VC$weights * d2l.be2.teta.st

d2l.sigma.st.sigma.st <- VC$weights * d2l.sigma.st.sigma.st
d2l.sigma.st.teta.st  <- VC$weights * d2l.sigma.st.teta.st

d2l.teta.st.teta.st   <- VC$weights * d2l.teta.st.teta.st


##### Global score and Hessian matrix #####

colnames(dl.dc1_2        ) <- colnames(d2l.c1_1.c1_2) <- colnames(d2l.c1_2.c1_2    ) <-
colnames(d2l.c1_2.be1    ) <- colnames(d2l.c1_2.be2 ) <- colnames(d2l.c1_2.sigma.st) <-
colnames(d2l.c1_2.teta.st) <- NULL

rownames(d2l.c1_2.c1_2   ) <- NULL


if( is.null(VC$X3)){

c11.c11   <-      sum (             d2l.c1_1.c1_1    )
c11.c12   <- t(rowSums(t(           d2l.c1_1.c1_2  )))
c11.be1   <- t(rowSums(t(-VC$X1 * c(d2l.c1_1.be1  )))) 
c11.be2   <- t(rowSums(t( VC$X2 * c(d2l.c1_1.be2  ))))
c11.sigma <-      sum (             d2l.c1_1.sigma.st)
c11.teta  <-      sum (             d2l.c1_1.teta.st )

c12.c12   <-               d2l.c1_2.c1_2
c12.be1   <- crossprod(    d2l.c1_2.be1     , -VC$X1)
c12.be2   <- crossprod(    d2l.c1_2.be2     ,  VC$X2)
c12.sigma <- t(t(rowSums(t(d2l.c1_2.sigma.st     ))))
c12.teta  <- t(t(rowSums(t(d2l.c1_2.teta.st      ))))

be1.be1   <- crossprod(     VC$X1 * c(d2l.be1.be1     ), VC$X1)
be1.be2   <- crossprod(    -VC$X1 * c(d2l.be1.be2     ), VC$X2)
be1.sigma <- t(t(rowSums(t(-VC$X1 * c(d2l.be1.sigma.st)    ))))
be1.teta  <- t(t(rowSums(t(-VC$X1 * c(d2l.be1.teta.st )    ))))

be2.be2   <- crossprod(    VC$X2 * c(d2l.be2.be2     ), VC$X2)
be2.sigma <- t(t(rowSums(t(VC$X2 * c(d2l.be2.sigma.st)    ))))
be2.teta  <- t(t(rowSums(t(VC$X2 * c(d2l.be2.teta.st )    ))))

sigma.sigma <- sum(d2l.sigma.st.sigma.st)
sigma.teta  <- sum(d2l.sigma.st.teta.st )

teta.teta <- sum(d2l.teta.st.teta.st)


G <- - c(   sum (  dl.dc1_1              ),
         colSums(  dl.dc1_2              ),
         colSums(c(dl.dbe1     ) * -VC$X1),
         colSums(c(dl.dbe2     ) *  VC$X2),
            sum (  dl.dsigma.st          ), 
            sum (  dl.dteta.st           ) )

H <- - rbind(
cbind(  c11.c11   ,   c11.c12   ,   c11.be1   ,   c11.be2   ,   c11.sigma    , c11.teta  ),
cbind(t(c11.c12  ),   c12.c12   ,   c12.be1   ,   c12.be2   ,   c12.sigma    , c12.teta  ),
cbind(t(c11.be1  ), t(c12.be1  ),   be1.be1   ,   be1.be2   ,   be1.sigma    , be1.teta  ),
cbind(t(c11.be2  ), t(c12.be2  ), t(be1.be2  ),   be2.be2   ,   be2.sigma    , be2.teta  ),
cbind(t(c11.sigma), t(c12.sigma), t(be1.sigma), t(be2.sigma),   sigma.sigma  , sigma.teta),
cbind(t(c11.teta ), t(c12.teta ), t(be1.teta ), t(be2.teta ), t(sigma.teta  ), teta.teta ) )

}


if( !(is.null(VC$X3)) ){

c11.c11   <-      sum (             d2l.c1_1.c1_1       )
c11.c12   <- t(rowSums(t(           d2l.c1_1.c1_2     )))
c11.be1   <- t(rowSums(t(-VC$X1 * c(d2l.c1_1.be1     )))) 
c11.be2   <- t(rowSums(t( VC$X2 * c(d2l.c1_1.be2     ))))
c11.sigma <- t(rowSums(t( VC$X3 * c(d2l.c1_1.sigma.st))))
c11.teta  <- t(rowSums(t( VC$X4 * c(d2l.c1_1.teta.st ))))

c12.c12   <-           d2l.c1_2.c1_2
c12.be1   <- crossprod(d2l.c1_2.be1     , -VC$X1)
c12.be2   <- crossprod(d2l.c1_2.be2     ,  VC$X2)
c12.sigma <- crossprod(d2l.c1_2.sigma.st,  VC$X3)
c12.teta  <- crossprod(d2l.c1_2.teta.st ,  VC$X4)

be1.be1   <- crossprod( VC$X1 * c(d2l.be1.be1     ), VC$X1)
be1.be2   <- crossprod(-VC$X1 * c(d2l.be1.be2     ), VC$X2)
be1.sigma <- crossprod(-VC$X1 * c(d2l.be1.sigma.st), VC$X3)
be1.teta  <- crossprod(-VC$X1 * c(d2l.be1.teta.st ), VC$X4)

be2.be2   <- crossprod(VC$X2 * c(d2l.be2.be2     ), VC$X2)
be2.sigma <- crossprod(VC$X2 * c(d2l.be2.sigma.st), VC$X3)
be2.teta  <- crossprod(VC$X2 * c(d2l.be2.teta.st ), VC$X4)

sigma.sigma <- crossprod(VC$X3 * c(d2l.sigma.st.sigma.st), VC$X3)
sigma.teta  <- crossprod(VC$X3 * c(d2l.sigma.st.teta.st ), VC$X4)

teta.teta <- crossprod(VC$X4 * c(d2l.teta.st.teta.st), VC$X4)


G <- - c(   sum (dl.dc1_1                ),
         colSums(dl.dc1_2                ),
         colSums(c(dl.dbe1)      * -VC$X1),
         colSums(c(dl.dbe2)      *  VC$X2),
         colSums(c(dl.dsigma.st) *  VC$X3), 
         colSums(c(dl.dteta.st)  *  VC$X4) )

H <- - rbind(
cbind(  c11.c11   ,   c11.c12   ,   c11.be1   ,   c11.be2   ,   c11.sigma    , c11.teta  ),
cbind(t(c11.c12  ),   c12.c12   ,   c12.be1   ,   c12.be2   ,   c12.sigma    , c12.teta  ),
cbind(t(c11.be1  ), t(c12.be1  ),   be1.be1   ,   be1.be2   ,   be1.sigma    , be1.teta  ),
cbind(t(c11.be2  ), t(c12.be2  ), t(be1.be2  ),   be2.be2   ,   be2.sigma    , be2.teta  ),
cbind(t(c11.sigma), t(c12.sigma), t(be1.sigma), t(be2.sigma),   sigma.sigma  , sigma.teta),
cbind(t(c11.teta ), t(c12.teta ), t(be1.teta ), t(be2.teta ), t(sigma.teta  ), teta.teta ) )

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

	S.h <- S.h[-VC$drop.ind, -VC$drop.ind]
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


if(VC$margins[2] == "LN"){

dHs  <- distrHsAT(exp(respvec$y2), eta2, sigma2, 1, margin2 = VC$margins[2])
pdf2 <- dHs$pdf2
p2   <- dHs$p2 


h.p <- h.m <- NA

if(length(teta1) != 0){

dH1.p <- copgHsAT(p1.p[teta.ind1], p2[teta.ind1], teta1, Cop1) 
dH1.m <- copgHsAT(p1.m[teta.ind1], p2[teta.ind1], teta1, Cop1) 

h.p[teta.ind1] <- dH1.p$c.copula.be2 ; h.m[teta.ind1] <- dH1.m$c.copula.be2

}

if(length(teta2) != 0){

dH2.p <- copgHsAT(p1.p[teta.ind2], p2[teta.ind2], teta2, Cop2) 
dH2.m <- copgHsAT(p1.m[teta.ind2], p2[teta.ind2], teta2, Cop2)

h.p[teta.ind2] <- dH2.p$c.copula.be2 ; h.m[teta.ind2] <- dH2.m$c.copula.be2

}

h.pm  <- pmax(epsilon, h.p - h.m)


l.ln <- - sum( VC$weights*(log(h.pm) + log(pdf2)) )

}


list(value = res, gradient = G, hessian = H, S.h = S.h, S.h1 = S.h1, S.h2 = S.h2, l = S.res, l.ln = l.ln, l.par = l.par, ps = ps, etas = etas,
     eta1 = eta1, eta2 = eta2, etad = etad, lp1 = lp1,
     dl.dbe1 = dl.dbe1, dl.dbe2 = dl.dbe2, dl.dsigma.st = dl.dsigma.st, dl.dteta.st = dl.dteta.st,
     BivD = VC$BivD, p1.p = p1.p, p1.m = p1.m, p2 = p2, theta.star = teta.st, # why 1-p1? Then p1 has been substituted with p1.p and p1.m
     teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
     Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2)

}