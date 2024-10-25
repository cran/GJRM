bcontSurvG_extended <- function(params, respvec, VC, ps, AT = FALSE){
        
        p1 <- p2 <- pdf1 <- pdf2 <- c.copula.be2 <- c.copula.be1 <- c.copula2.be1be2 <- NA
        
        monP <- monP1 <- k1 <- k2 <- 0; Veq1 <- Veq2 <- list()
        monP2 <- matrix(0, length(params),length(params))
        
        rotConst <- 1
        
        params1 <- params[1:VC$X1.d2]
        params2 <- params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
        
        params1[VC$mono.sm.pos1] <- exp( params1[VC$mono.sm.pos1] )
        params2[VC$mono.sm.pos2] <- exp( params2[VC$mono.sm.pos2] )
        
        ##########################
        
        # eta is the design multiplied by the beta vector
        eta1 <- VC$X1%*%params1
        eta2 <- VC$X2%*%params2
        
        # NEW quantities [for the parts interval cencored]
        
        eta1.2 <- VC$X1.2%*%params1
        eta2.2 <- VC$X2.2%*%params2
        
        #### NEW ######
        
        etad <- etas1 <- etas2 <- l.ln <- NULL
        
        # derivatives of eta 1 and eta 2
        Xd1P <- VC$Xd1%*%params1
        Xd2P <- VC$Xd2%*%params2
        
        
        etad <- etas1 <- etas2 <- l.ln <- NULL
        
        if( is.null(VC$X3) ){
                X3 <- matrix(1, VC$n, 1)
                teta.st <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
        }
        
        if( !is.null(VC$X3) ){
                X3 <- VC$X3
                teta.st <- etad <- X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
        }
        
        ##################
        
        indNeq1 <- as.numeric(Xd1P < 0)
        indNeq2 <- as.numeric(Xd2P < 0)
        
        Xd1P <- ifelse(Xd1P < VC$min.dn, VC$min.dn, Xd1P ) # check that the derivatives are lower that a certain value
        Xd2P <- ifelse(Xd2P < VC$min.dn, VC$min.dn, Xd2P )
        
        ##################
        
        
        
        ##################
        ## Transformations
        ##################
        
        resT  <- teta.tr(VC, teta.st)
        
        teta.st1 <- teta.st2 <- teta.st <- resT$teta.st
        teta1 <- teta2 <- teta <- resT$teta
        
        ##################
        
        Cop1 <- Cop2 <- VC$BivD
        nC1 <- nC2 <- VC$nC
        
        teta.ind1 <- as.logical(c(1,0,round(runif(VC$n-2))) )
        teta.ind2 <- teta.ind1 == FALSE
        
        
        if(!(VC$BivD %in% VC$BivD2) && length(teta.st) > 1){
                
                teta.st1 <- teta.st[teta.ind1]
                teta.st2 <- teta.st[teta.ind2]
                
                teta1 <- teta[teta.ind1]
                teta2 <- teta[teta.ind2]
                
        }
        
        ###
        
        if(VC$BivD %in% VC$BivD2){
                
                if(VC$BivD %in% VC$BivD2[c(1:4,13:16)])  teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov), TRUE, FALSE)
                if(VC$BivD %in% VC$BivD2[5:12]) teta.ind1 <- ifelse(VC$my.env$signind*teta > exp(VC$zerov) + 1, TRUE, FALSE)
                teta.ind2 <- teta.ind1 == FALSE
                
                VC$my.env$signind <- ifelse(teta.ind1 == TRUE,  1, -1)
                
                teta1 <-  teta[teta.ind1]
                teta2 <- -teta[teta.ind2]
                
                teta.st1 <- teta.st[teta.ind1]
                teta.st2 <- teta.st[teta.ind2]
                
                if(length(teta) == 1) teta.ind2 <- teta.ind1 <- rep(TRUE, VC$n)
                
                Cop1Cop2R <- Cop1Cop2(VC$BivD)
                Cop1 <- Cop1Cop2R$Cop1
                Cop2 <- Cop1Cop2R$Cop2
                
                nC1 <- VC$ct[which(VC$ct[,1] == Cop1),2]
                nC2 <- VC$ct[which(VC$ct[,1] == Cop2),2]
                
        }
        
        
        
        ##################
        
        pd1 <- probmS(eta1, VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        pd2 <- probmS(eta2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        
        ###################
        ### NEW ###########
        ###################
        ## new quantities derived from eta1.2 and eta2.2##
        ## for the interval cencored parts              ##
        
        # the interval quantities of probms are for control only###
        
        pd1.2 <- probmS(eta1.2, VC$margins[1], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        pd2.2 <- probmS(eta2.2, VC$margins[2], min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr) #modificato errore su margine
        
        
        ###
        
        
        p1 <- pd1$pr
        p2 <- pd2$pr
        
        ### NEW ###
        # quantities derived from eta1.2 and eta2.2
        
        #p1.2 survival primo e secondo argomento calcolate per individui interval nel tempo right.
        
        
        p1.2 <- pd1.2$pr
        p2.2 <- pd2.2$pr
        
        ####### derivatives respect to time eta1 and eta2
        
        dS1eta1 <- pd1$dS
        dS2eta2 <- pd2$dS
        
        ##NEW#####
        # 1st derivative eta1.2 and eta2.2
        
        dS1eta1.2 <- pd1.2$dS
        dS2eta2.2 <- pd2.2$dS
        #
        
        d2S1eta1 <- pd1$d2S
        d2S2eta2 <- pd2$d2S
        
        ##NEW###
        #2nd derivative###
        
        d2S1eta1.2 <- pd1.2$d2S
        d2S2eta2.2 <- pd2.2$d2S
        
        #
        
        d3S1eta1 <- pd1$d3S
        d3S2eta2 <- pd2$d3S
        
        ### NEW ######
        # 3rd derivative ##
        
        d3S1eta1.2 <- pd1.2$d3S
        d3S2eta2.2 <- pd2.2$d3S
        
        ##################
        
        if( length(teta1) != 0) dH1 <- copgHs(p1[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        if( length(teta2) != 0) dH2 <- copgHs(p1[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        
        ## NEW####
        ## Quantities addded
        # dH1.2 {p1.2[theta1], p2.2[theta1]}
        # dH1.mix1 {p1[theta1], p2.2[theta1]}
        # dH1.mix2 {p1.2[theta1], p2[theta1]}
        
        ######
        #dH2.2 {p1.2[theta2], p2.2[theta2]}
        #dH2.mix1 {p1[theta2], p2.2[theta2]}
        #dH2.mix2 {p.1[theta2], p2[theta2]}
        
        
        if( length(teta1) != 0) dH1.2 <- copgHs(p1.2[teta.ind1], p2.2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        if( length(teta2) != 0) dH2.2 <- copgHs(p1.2[teta.ind2], p2.2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        
        ##NEW##########
        ####### queste sono le quantita miste
        if( length(teta1) != 0) dH1.mix1 <- copgHs(p1[teta.ind1], p2.2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        if( length(teta2) != 0) dH2.mix1 <- copgHs(p1[teta.ind2], p2.2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        
        if( length(teta1) != 0) dH1.mix2 <- copgHs(p1.2[teta.ind1], p2[teta.ind1], eta1=NULL, eta2=NULL, teta1, teta.st1, Cop1, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        if( length(teta2) != 0) dH2.mix2 <- copgHs(p1.2[teta.ind2], p2[teta.ind2], eta1=NULL, eta2=NULL, teta2, teta.st2, Cop2, VC$dof, min.dn = VC$min.dn, min.pr = VC$min.pr, max.pr = VC$max.pr)
        
        
        c.copula2.be1be2 <- c.copula.be1 <- c.copula.be2 <- p00 <- c.copula.theta <- c.copula.thet <- bit1.th2ATE <- NA
        c.copula2.be1be2.2 <- c.copula.be1.2 <- c.copula.be2.2 <- p00.2 <- c.copula.theta.2 <- c.copula.thet.2 <- bit1.th2ATE.2 <- NA
        c.copula2.be1be2.mix1 <- c.copula.be1.mix1 <- c.copula.be2.mix1 <- p00.mix1 <- c.copula.theta.mix1 <- c.copula.thet.mix1 <- bit1.th2ATE.mix1 <- NA
        c.copula2.be1be2.mix2 <- c.copula.be1.mix2 <- c.copula.be2.mix2 <- p00.mix2 <- c.copula.theta.mix2 <- c.copula.thet.mix2 <- bit1.th2ATE.mix2 <- NA
        #  dipendenze non simmetriche copula [ rotazione 180 gradi]
        
        
        if( length(teta1) != 0){
                c.copula2.be1be2[teta.ind1] <- dH1$c.copula2.be1be2
                c.copula.be1[teta.ind1]     <- dH1$c.copula.be1
                c.copula.be2[teta.ind1]     <- dH1$c.copula.be2
                c.copula.theta[teta.ind1]   <- dH1$c.copula.theta
                c.copula.thet[teta.ind1]    <- dH1$c.copula.thet
                bit1.th2ATE[teta.ind1]      <- dH1$bit1.th2ATE
                p00[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
        }
        
        ###NEW####
        #Extract from: dH1.mix1
        # prefix :     SOMETHING.mix1
        
        if( length(teta1) != 0){
                c.copula2.be1be2.mix1[teta.ind1] <- dH1.mix1$c.copula2.be1be2
                c.copula.be1.mix1[teta.ind1]     <- dH1.mix1$c.copula.be1
                c.copula.be2.mix1[teta.ind1]     <- dH1.mix1$c.copula.be2
                c.copula.theta.mix1[teta.ind1]   <- dH1.mix1$c.copula.theta
                c.copula.thet.mix1[teta.ind1]    <- dH1.mix1$c.copula.thet
                bit1.th2ATE.mix1[teta.ind1]      <- dH1.mix1$bit1.th2ATE
                p00.mix1[teta.ind1] <- mm(BiCDF(p1[teta.ind1], p2.2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
                
        }
        
        ###NEW#####
        #Extract from: dH1.mix2
        # prefix: SOMETHING.mix2
        if( length(teta1) != 0){
                c.copula2.be1be2.mix2[teta.ind1] <- dH1.mix2$c.copula2.be1be2
                c.copula.be1.mix2[teta.ind1]     <- dH1.mix2$c.copula.be1
                c.copula.be2.mix2[teta.ind1]     <- dH1.mix2$c.copula.be2
                c.copula.theta.mix2[teta.ind1]   <- dH1.mix2$c.copula.theta
                c.copula.thet.mix2[teta.ind1]    <- dH1.mix2$c.copula.thet
                bit1.th2ATE.mix2[teta.ind1]      <- dH1.mix2$bit1.th2ATE
                p00.mix2[teta.ind1] <- mm(BiCDF(p1.2[teta.ind1], p2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
                
        }
        
        ##NEW####
        #Extract from: dH1.2
        # prefix: SOMETHING.2
        if( length(teta1) != 0){
                c.copula2.be1be2.2[teta.ind1] <- dH1.2$c.copula2.be1be2
                c.copula.be1.2[teta.ind1]     <- dH1.2$c.copula.be1
                c.copula.be2.2[teta.ind1]     <- dH1.2$c.copula.be2
                c.copula.theta.2[teta.ind1]   <- dH1.2$c.copula.theta
                c.copula.thet.2[teta.ind1]    <- dH1.2$c.copula.thet
                bit1.th2ATE.2[teta.ind1]      <- dH1.2$bit1.th2ATE
                p00.2[teta.ind1] <- mm(BiCDF(p1.2[teta.ind1], p2.2[teta.ind1], nC1, teta1, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
                
        }
        
        
        
        if( length(teta2) != 0){
                c.copula2.be1be2[teta.ind2] <- dH2$c.copula2.be1be2
                c.copula.be1[teta.ind2]     <- dH2$c.copula.be1
                c.copula.be2[teta.ind2]     <- dH2$c.copula.be2
                c.copula.theta[teta.ind2]   <- dH2$c.copula.theta
                c.copula.thet[teta.ind2]    <- dH2$c.copula.thet
                bit1.th2ATE[teta.ind2]      <- dH2$bit1.th2ATE
                p00[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
        }
        
        ###NEW###
        #Extract from: dH2.mix1
        # Prefix: SOMETHING.mix1
        if( length(teta2) != 0){
                c.copula2.be1be2.mix1[teta.ind2] <- dH2.mix1$c.copula2.be1be2
                c.copula.be1.mix1[teta.ind2]     <- dH2.mix1$c.copula.be1
                c.copula.be2.mix1[teta.ind2]     <- dH2.mix1$c.copula.be2
                c.copula.theta.mix1[teta.ind2]   <- dH2.mix1$c.copula.theta
                c.copula.thet.mix1[teta.ind2]    <- dH2.mix1$c.copula.thet
                bit1.th2ATE.mix1[teta.ind2]      <- dH2.mix1$bit1.th2ATE
                p00.mix1[teta.ind2] <- mm(BiCDF(p1[teta.ind2], p2.2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
        }
        
        ###NEW####
        #Extract from: dH2.mix2
        # Prefix: SOMETHING.mix2
        if( length(teta2) != 0){
                c.copula2.be1be2.mix2[teta.ind2] <- dH2.mix2$c.copula2.be1be2
                c.copula.be1.mix2[teta.ind2]     <- dH2.mix2$c.copula.be1
                c.copula.be2.mix2[teta.ind2]     <- dH2.mix2$c.copula.be2
                c.copula.theta.mix2[teta.ind2]   <- dH2.mix2$c.copula.theta
                c.copula.thet.mix2[teta.ind2]    <- dH2.mix2$c.copula.thet
                bit1.th2ATE.mix2[teta.ind2]      <- dH2.mix2$bit1.th2ATE
                p00.mix2[teta.ind2] <- mm(BiCDF(p1.2[teta.ind2], p2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
        }
        ##
        #Extract from: dH2.2
        # Prefix: SOMETHING.2
        if( length(teta2) != 0){
                c.copula2.be1be2.2[teta.ind2] <- dH2.2$c.copula2.be1be2
                c.copula.be1.2[teta.ind2]     <- dH2.2$c.copula.be1
                c.copula.be2.2[teta.ind2]     <- dH2.2$c.copula.be2
                c.copula.theta.2[teta.ind2]   <- dH2.2$c.copula.theta
                c.copula.thet.2[teta.ind2]    <- dH2.2$c.copula.thet
                bit1.th2ATE.2[teta.ind2]      <- dH2.2$bit1.th2ATE
                p00.2[teta.ind2] <- mm(BiCDF(p1.2[teta.ind2], p2.2[teta.ind2], nC2, teta2, VC$dof), min.pr = VC$min.pr, max.pr = VC$max.pr  )
        }
        ################################
        ################################
        ################################
        
        #Pieces for the Gradient 
        ##################
        
        der.par1 <- der2.par1 <- params1; der.par2 <- der2.par2 <- params2
        
        der.par1[-c( VC$mono.sm.pos1 )] <- 1
        der.par2[-c( VC$mono.sm.pos2 )] <- 1
        
        der2.par1[-c( VC$mono.sm.pos1 )] <- 0
        der2.par2[-c( VC$mono.sm.pos2 )] <- 0
        
        
        der2eta1dery1b1 <- t(t(VC$Xd1)*der.par1)
        der2eta2dery2b2 <- t(t(VC$Xd2)*der.par2)
        
        dereta1derb1    <- t(t(VC$X1)*der.par1)
        dereta2derb2    <- t(t(VC$X2)*der.par2)
        
        #New##
        
        dereta1derb1.2    <- t(t(VC$X1.2)*der.par1)
        dereta2derb2.2    <- t(t(VC$X2.2)*der.par2)
        
        ##################
        ##################
        
        # STANDARD
        if( length(teta1) != 0) BITS1 <- copgHsCont(p1[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        if( length(teta2) != 0) BITS2 <- copgHsCont(p1[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        
        ###NEW####
        # p1 p2.2
        # prefix :     SOMETHING.mix1
        if( length(teta1) != 0) BITS1.mix1 <- copgHsCont(p1[teta.ind1], p2.2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        if( length(teta2) != 0) BITS2.mix1 <- copgHsCont(p1[teta.ind2], p2.2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        
        
        ###NEW####
        # p1.2 p2
        # prefix :     SOMETHING.mix2
        if( length(teta1) != 0) BITS1.mix2 <- copgHsCont(p1.2[teta.ind1], p2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        if( length(teta2) != 0) BITS2.mix2 <- copgHsCont(p1.2[teta.ind2], p2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        
        
        ###NEW####
        # p1.2 p2.2
        # prefix :     SOMETHING.2
        
        if( length(teta1) != 0) BITS1.2 <- copgHsCont(p1.2[teta.ind1], p2.2[teta.ind1], teta1, teta.st1, Cop1, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        if( length(teta2) != 0) BITS2.2 <- copgHsCont(p1.2[teta.ind2], p2.2[teta.ind2], teta2, teta.st2, Cop2, Cont = TRUE, par2 = VC$dof, nu.st = log(VC$dof - 2))
        
        ######## Modifica der2h.derp1p1 ##########
        der2h.derp1p1 <- NA
        if( length(teta1) != 0) der2h.derp1p1[teta.ind1]   <- BITS1$der2h.derp1p1
        if( length(teta2) != 0) der2h.derp1p1[teta.ind2]   <- BITS2$der2h.derp1p1
        
        # new ###
        ###NEW####
        # p1 p2.2
        # prefix :     SOMETHING.mix1
        
        der2h.derp1p1.mix1 <- NA
        if( length(teta1) != 0) der2h.derp1p1.mix1[teta.ind1]   <- BITS1.mix1$der2h.derp1p1
        if( length(teta2) != 0) der2h.derp1p1.mix1[teta.ind2]   <- BITS2.mix1$der2h.derp1p1
        
        #NEW#
        # p1.2 p2
        #prefix : something.mix2
        der2h.derp1p1.mix2 <- NA
        if( length(teta1) != 0) der2h.derp1p1.mix2[teta.ind1]   <- BITS1.mix2$der2h.derp1p1
        if( length(teta2) != 0) der2h.derp1p1.mix2[teta.ind2]   <- BITS2.mix2$der2h.derp1p1
        
        #NEW#
        #p1.2 p2.2
        # prefix: something:2
        der2h.derp1p1.2 <- NA
        if( length(teta1) != 0) der2h.derp1p1.2[teta.ind1]   <- BITS1.2$der2h.derp1p1
        if( length(teta2) != 0) der2h.derp1p1.2[teta.ind2]   <- BITS2.2$der2h.derp1p1
        
        
        
        ######## Modifica der2h.derp1p2 ##########
        
        der2h.derp1p2 <- NA
        if( length(teta1) != 0) der2h.derp1p2[teta.ind1] <- BITS1$der2h.derp1p2
        if( length(teta2) != 0) der2h.derp1p2[teta.ind2] <- BITS2$der2h.derp1p2
        
        ## NEW###
        # prefix: mix1
        der2h.derp1p2.mix1 <- NA
        if( length(teta1) != 0) der2h.derp1p2.mix1[teta.ind1] <- BITS1.mix1$der2h.derp1p2
        if( length(teta2) != 0) der2h.derp1p2.mix1[teta.ind2] <- BITS2.mix1$der2h.derp1p2
        
        # NEW ####
        # prefix: mix2
        
        der2h.derp1p2.mix2 <- NA
        if( length(teta1) != 0) der2h.derp1p2.mix2[teta.ind1] <- BITS1.mix2$der2h.derp1p2
        if( length(teta2) != 0) der2h.derp1p2.mix2[teta.ind2] <- BITS2.mix2$der2h.derp1p2
        
        # NEW ##
        # Prefix: 2
        der2h.derp1p2.2 <- NA
        if( length(teta1) != 0) der2h.derp1p2.2[teta.ind1] <- BITS1.2$der2h.derp1p2
        if( length(teta2) != 0) der2h.derp1p2.2[teta.ind2] <- BITS2.2$der2h.derp1p2
        
        ######## Modifica der2h.derp1teta & derteta.derteta.st   ##########
        
        der2h.derp1teta            <- NA
        derteta.derteta.st         <- NA
        if( length(teta1) != 0) der2h.derp1teta[teta.ind1]    <- BITS1$der2h.derp1teta
        if( length(teta2) != 0) der2h.derp1teta[teta.ind2]    <- BITS2$der2h.derp1teta
        if( length(teta1) != 0) derteta.derteta.st[teta.ind1] <- BITS1$derteta.derteta.st
        if( length(teta2) != 0) derteta.derteta.st[teta.ind2] <- BITS2$derteta.derteta.st
        
        der2h.derp1teta.st  <- der2h.derp1teta * derteta.derteta.st # new bit
        
        # NEW ##
        # prefix: mix 1
        
        der2h.derp1teta.mix1           <- NA
        derteta.derteta.st.mix1         <- NA
        if( length(teta1) != 0) der2h.derp1teta.mix1[teta.ind1]    <- BITS1.mix1$der2h.derp1teta
        if( length(teta2) != 0) der2h.derp1teta.mix1[teta.ind2]    <- BITS2.mix1$der2h.derp1teta
        if( length(teta1) != 0) derteta.derteta.st.mix1[teta.ind1] <- BITS1.mix1$derteta.derteta.st
        if( length(teta2) != 0) derteta.derteta.st.mix1[teta.ind2] <- BITS2.mix1$derteta.derteta.st
        
        der2h.derp1teta.st.mix1  <- der2h.derp1teta.mix1 * derteta.derteta.st.mix1 # new bit
        
        # NEW ##
        # prefix: mix2
        
        der2h.derp1teta.mix2            <- NA
        derteta.derteta.st.mix2         <- NA
        if( length(teta1) != 0) der2h.derp1teta.mix2[teta.ind1]    <- BITS1.mix2$der2h.derp1teta
        if( length(teta2) != 0) der2h.derp1teta.mix2[teta.ind2]    <- BITS2.mix2$der2h.derp1teta
        if( length(teta1) != 0) derteta.derteta.st.mix2[teta.ind1] <- BITS1.mix2$derteta.derteta.st
        if( length(teta2) != 0) derteta.derteta.st.mix2[teta.ind2] <- BITS2.mix2$derteta.derteta.st
        
        der2h.derp1teta.st.mix2  <- der2h.derp1teta.mix2 * derteta.derteta.st.mix2 # new bit
        
        #New##
        # prefix : .2
        
        der2h.derp1teta.2            <- NA
        derteta.derteta.st.2         <- NA
        if( length(teta1) != 0) der2h.derp1teta.2[teta.ind1]    <- BITS1.2$der2h.derp1teta
        if( length(teta2) != 0) der2h.derp1teta.2[teta.ind2]    <- BITS2.2$der2h.derp1teta
        if( length(teta1) != 0) derteta.derteta.st.2[teta.ind1] <- BITS1.2$derteta.derteta.st
        if( length(teta2) != 0) derteta.derteta.st.2[teta.ind2] <- BITS2.2$derteta.derteta.st
        
        der2h.derp1teta.st.2  <- der2h.derp1teta.2 * derteta.derteta.st.2 # new bit
        
        #################
        #################
        # creare le quantit mix1 mix2 e 2
        
        c.copula2.be1 <- c.copula2.be2 <- c.copula2.be1th <- c.copula2.be2th <- bit1.th2 <- c.copula2.be1t <- c.copula2.be2t <- NA
        
        if( length(teta1) != 0){
                
                c.copula2.be1[teta.ind1]    <- dH1$c.copula2.be1
                c.copula2.be2[teta.ind1]    <- dH1$c.copula2.be2
                c.copula2.be1th[teta.ind1]  <- dH1$c.copula2.be1th
                c.copula2.be2th[teta.ind1]  <- dH1$c.copula2.be2th
                c.copula2.be1t[teta.ind1]   <- dH1$c.copula2.be1t
                c.copula2.be2t[teta.ind1]   <- dH1$c.copula2.be2t
                bit1.th2[teta.ind1]         <- dH1$bit1.th2
                
        }
        
        #NEW#
        #Prefix: mix1
        
        c.copula2.be1.mix1 <- c.copula2.be2.mix1 <- c.copula2.be1th.mix1 <- c.copula2.be2th.mix1 <- bit1.th2.mix1 <- c.copula2.be1t.mix1 <- c.copula2.be2t.mix1 <- NA
        if( length(teta1) != 0){
                
                c.copula2.be1.mix1[teta.ind1]    <- dH1.mix1$c.copula2.be1
                c.copula2.be2.mix1[teta.ind1]    <- dH1.mix1$c.copula2.be2
                c.copula2.be1th.mix1[teta.ind1]  <- dH1.mix1$c.copula2.be1th
                c.copula2.be2th.mix1[teta.ind1]  <- dH1.mix1$c.copula2.be2th
                c.copula2.be1t.mix1[teta.ind1]   <- dH1.mix1$c.copula2.be1t
                c.copula2.be2t.mix1[teta.ind1]   <- dH1.mix1$c.copula2.be2t
                bit1.th2.mix1[teta.ind1]         <- dH1.mix1$bit1.th2
                
        }
        
        #NEW#
        #Prefix : mix2
        c.copula2.be1.mix2 <- c.copula2.be2.mix2 <- c.copula2.be1th.mix2 <- c.copula2.be2th.mix2 <- bit1.th2.mix2 <- c.copula2.be1t.mix2 <- c.copula2.be2t.mix2 <- NA
        
        if( length(teta1) != 0){
                
                c.copula2.be1.mix2[teta.ind1]    <- dH1.mix2$c.copula2.be1
                c.copula2.be2.mix2[teta.ind1]    <- dH1.mix2$c.copula2.be2
                c.copula2.be1th.mix2[teta.ind1]  <- dH1.mix2$c.copula2.be1th
                c.copula2.be2th.mix2[teta.ind1]  <- dH1.mix2$c.copula2.be2th
                c.copula2.be1t.mix2[teta.ind1]   <- dH1.mix2$c.copula2.be1t
                c.copula2.be2t.mix2[teta.ind1]   <- dH1.mix2$c.copula2.be2t
                bit1.th2.mix2[teta.ind1]         <- dH1.mix2$bit1.th2
                
        }
        
        #NEW#
        #prefix: 2
        c.copula2.be1.2 <- c.copula2.be2.2 <- c.copula2.be1th.2 <- c.copula2.be2th.2 <- bit1.th2.2 <- c.copula2.be1t.2 <- c.copula2.be2t.2 <- NA
        
        if( length(teta1) != 0){
                
                c.copula2.be1.2[teta.ind1]    <- dH1.2$c.copula2.be1
                c.copula2.be2.2[teta.ind1]    <- dH1.2$c.copula2.be2
                c.copula2.be1th.2[teta.ind1]  <- dH1.2$c.copula2.be1th
                c.copula2.be2th.2[teta.ind1]  <- dH1.2$c.copula2.be2th
                c.copula2.be1t.2[teta.ind1]   <- dH1.2$c.copula2.be1t
                c.copula2.be2t.2[teta.ind1]   <- dH1.2$c.copula2.be2t
                bit1.th2.2[teta.ind1]         <- dH1.2$bit1.th2
                
        }
        
        
        if( length(teta2) != 0){
                
                c.copula2.be1[teta.ind2]    <- dH2$c.copula2.be1
                c.copula2.be2[teta.ind2]    <- dH2$c.copula2.be2
                c.copula2.be1th[teta.ind2]  <- dH2$c.copula2.be1th
                c.copula2.be2th[teta.ind2]  <- dH2$c.copula2.be2th
                c.copula2.be1t[teta.ind2]   <- dH2$c.copula2.be1t
                c.copula2.be2t[teta.ind2]   <- dH2$c.copula2.be2t
                bit1.th2[teta.ind2]         <- dH2$bit1.th2
                
        }
        
        #NEW#
        #prefix: mix1
        
        if( length(teta2) != 0){
                
                c.copula2.be1.mix1[teta.ind2]    <- dH2.mix1$c.copula2.be1
                c.copula2.be2.mix1[teta.ind2]    <- dH2.mix1$c.copula2.be2
                c.copula2.be1th.mix1[teta.ind2]  <- dH2.mix1$c.copula2.be1th
                c.copula2.be2th.mix1[teta.ind2]  <- dH2.mix1$c.copula2.be2th
                c.copula2.be1t.mix1[teta.ind2]   <- dH2.mix1$c.copula2.be1t
                c.copula2.be2t.mix1[teta.ind2]   <- dH2.mix1$c.copula2.be2t
                bit1.th2.mix1[teta.ind2]         <- dH2.mix1$bit1.th2
                
        }
        
        #NEW#
        #prefix: mix2
        if( length(teta2) != 0){
                
                c.copula2.be1.mix2[teta.ind2]    <- dH2.mix2$c.copula2.be1
                c.copula2.be2.mix2[teta.ind2]    <- dH2.mix2$c.copula2.be2
                c.copula2.be1th.mix2[teta.ind2]  <- dH2.mix2$c.copula2.be1th
                c.copula2.be2th.mix2[teta.ind2]  <- dH2.mix2$c.copula2.be2th
                c.copula2.be1t.mix2[teta.ind2]   <- dH2.mix2$c.copula2.be1t
                c.copula2.be2t.mix2[teta.ind2]   <- dH2.mix2$c.copula2.be2t
                bit1.th2.mix2[teta.ind2]         <- dH2.mix2$bit1.th2
                
        }
        
        #NEW#
        #prefix: 2
        
        if( length(teta2) != 0){
                
                c.copula2.be1.2[teta.ind2]    <- dH2.2$c.copula2.be1
                c.copula2.be2.2[teta.ind2]    <- dH2.2$c.copula2.be2
                c.copula2.be1th.2[teta.ind2]  <- dH2.2$c.copula2.be1th
                c.copula2.be2th.2[teta.ind2]  <- dH2.2$c.copula2.be2th
                c.copula2.be1t.2[teta.ind2]   <- dH2.2$c.copula2.be1t
                c.copula2.be2t.2[teta.ind2]   <- dH2.2$c.copula2.be2t
                bit1.th2.2[teta.ind2]         <- dH2.2$bit1.th2
                
        }
        #################
        #################
        ### Pieces for the Hessian 
        
        
        
        
        
        der2c.derrho.derrho    <- NA
        der2c.derp1.derp1      <- NA
        der2c.derp2.derp2      <- NA
        der2c.derp1.derp2      <- NA 
        der2c.derp1.derrho     <- NA
        der2c.derp2.derrho     <- NA
        der2teta.derteta.stteta.st <- NA 
        
        if( length(teta1) != 0){ 
                der2c.derrho.derrho[teta.ind1]    <- BITS1$der2c.derrho.derrho
                der2c.derp1.derp1[teta.ind1]      <- BITS1$der2c.derp1.derp1  
                der2c.derp2.derp2[teta.ind1]      <- BITS1$der2c.derp2.derp2  
                der2c.derp1.derp2[teta.ind1]      <- BITS1$der2c.derp1.derp2  
                der2c.derp1.derrho[teta.ind1]     <- BITS1$der2c.derp1.derrho 
                der2c.derp2.derrho[teta.ind1]     <- BITS1$der2c.derp2.derrho
        }
        
        # New -theta1#
        # Prefix .mix1
        der2c.derrho.derrho.mix1    <- NA
        der2c.derp1.derp1.mix1      <- NA
        der2c.derp2.derp2.mix1      <- NA
        der2c.derp1.derp2.mix1      <- NA 
        der2c.derp1.derrho.mix1     <- NA
        der2c.derp2.derrho.mix1     <- NA
        der2teta.derteta.stteta.st.mix1 <- NA 
        
        if( length(teta1) != 0){ 
                der2c.derrho.derrho.mix1[teta.ind1]    <- BITS1.mix1$der2c.derrho.derrho
                der2c.derp1.derp1.mix1[teta.ind1]      <- BITS1.mix1$der2c.derp1.derp1  
                der2c.derp2.derp2.mix1[teta.ind1]      <- BITS1.mix1$der2c.derp2.derp2  
                der2c.derp1.derp2.mix1[teta.ind1]      <- BITS1.mix1$der2c.derp1.derp2  
                der2c.derp1.derrho.mix1[teta.ind1]     <- BITS1.mix1$der2c.derp1.derrho 
                der2c.derp2.derrho.mix1[teta.ind1]     <- BITS1.mix1$der2c.derp2.derrho 
        }
        
        
        #New-theta1#
        #prefix .mix2
        
        der2c.derrho.derrho.mix2    <- NA
        der2c.derp1.derp1.mix2      <- NA
        der2c.derp2.derp2.mix2      <- NA
        der2c.derp1.derp2.mix2      <- NA 
        der2c.derp1.derrho.mix2     <- NA
        der2c.derp2.derrho.mix2     <- NA
        der2teta.derteta.stteta.st.mix2 <- NA 
        
        if( length(teta1) != 0){ 
                der2c.derrho.derrho.mix2[teta.ind1]    <- BITS1.mix2$der2c.derrho.derrho
                der2c.derp1.derp1.mix2[teta.ind1]      <- BITS1.mix2$der2c.derp1.derp1  
                der2c.derp2.derp2.mix2[teta.ind1]      <- BITS1.mix2$der2c.derp2.derp2  
                der2c.derp1.derp2.mix2[teta.ind1]      <- BITS1.mix2$der2c.derp1.derp2  
                der2c.derp1.derrho.mix2[teta.ind1]     <- BITS1.mix2$der2c.derp1.derrho 
                der2c.derp2.derrho.mix2[teta.ind1]     <- BITS1.mix2$der2c.derp2.derrho 
        }
        
        #New-theta1#
        #prefix .2
        
        der2c.derrho.derrho.2    <- NA
        der2c.derp1.derp1.2      <- NA
        der2c.derp2.derp2.2      <- NA
        der2c.derp1.derp2.2      <- NA 
        der2c.derp1.derrho.2     <- NA
        der2c.derp2.derrho.2     <- NA
        der2teta.derteta.stteta.st.2 <- NA 
        
        if( length(teta1) != 0){ 
                der2c.derrho.derrho.2[teta.ind1]    <- BITS1.2$der2c.derrho.derrho
                der2c.derp1.derp1.2[teta.ind1]      <- BITS1.2$der2c.derp1.derp1  
                der2c.derp2.derp2.2[teta.ind1]      <- BITS1.2$der2c.derp2.derp2  
                der2c.derp1.derp2.2[teta.ind1]      <- BITS1.2$der2c.derp1.derp2  
                der2c.derp1.derrho.2[teta.ind1]     <- BITS1.2$der2c.derp1.derrho 
                der2c.derp2.derrho.2[teta.ind1]     <- BITS1.2$der2c.derp2.derrho 
        }
        
        
        
        if( length(teta2) != 0){ 
                der2c.derrho.derrho[teta.ind2] <- BITS2$der2c.derrho.derrho
                der2c.derp1.derp1[teta.ind2]      <- BITS2$der2c.derp1.derp1  
                der2c.derp2.derp2[teta.ind2]      <- BITS2$der2c.derp2.derp2  
                der2c.derp1.derp2[teta.ind2]      <- BITS2$der2c.derp1.derp2  
                der2c.derp1.derrho[teta.ind2]     <- BITS2$der2c.derp1.derrho 
                der2c.derp2.derrho[teta.ind2]     <- BITS2$der2c.derp2.derrho
        }
        
        
        #New- theta 2#
        # prefix .mix1
        
        if( length(teta2) != 0){ 
                der2c.derrho.derrho.mix1[teta.ind2]    <- BITS2.mix1$der2c.derrho.derrho
                der2c.derp1.derp1.mix1[teta.ind2]      <- BITS2.mix1$der2c.derp1.derp1  
                der2c.derp2.derp2.mix1[teta.ind2]      <- BITS2.mix1$der2c.derp2.derp2  
                der2c.derp1.derp2.mix1[teta.ind2]      <- BITS2.mix1$der2c.derp1.derp2  
                der2c.derp1.derrho.mix1[teta.ind2]     <- BITS2.mix1$der2c.derp1.derrho 
                der2c.derp2.derrho.mix1[teta.ind2]     <- BITS2.mix1$der2c.derp2.derrho 
        }
        
        #New- theta2
        # prefix mix2
        
        if( length(teta2) != 0){ 
                der2c.derrho.derrho.mix2[teta.ind2]    <- BITS2.mix2$der2c.derrho.derrho
                der2c.derp1.derp1.mix2[teta.ind2]      <- BITS2.mix2$der2c.derp1.derp1  
                der2c.derp2.derp2.mix2[teta.ind2]      <- BITS2.mix2$der2c.derp2.derp2  
                der2c.derp1.derp2.mix2[teta.ind2]      <- BITS2.mix2$der2c.derp1.derp2  
                der2c.derp1.derrho.mix2[teta.ind2]     <- BITS2.mix2$der2c.derp1.derrho 
                der2c.derp2.derrho.mix2[teta.ind2]     <- BITS2.mix2$der2c.derp2.derrho
        }
        
        #New-theta2
        # prefix .2
        
        if( length(teta2) != 0){ 
                der2c.derrho.derrho.2[teta.ind2] <- BITS2.2$der2c.derrho.derrho
                der2c.derp1.derp1.2[teta.ind2]      <- BITS2.2$der2c.derp1.derp1  
                der2c.derp2.derp2.2[teta.ind2]      <- BITS2.2$der2c.derp2.derp2  
                der2c.derp1.derp2.2[teta.ind2]      <- BITS2.2$der2c.derp1.derp2  
                der2c.derp1.derrho.2[teta.ind2]     <- BITS2.2$der2c.derp1.derrho 
                der2c.derp2.derrho.2[teta.ind2]     <- BITS2.2$der2c.derp2.derrho
        }
        
        if( length(teta1) != 0) der2teta.derteta.stteta.st[teta.ind1] <- BITS1$der2teta.derteta.stteta.st
        #New-theta1
        # prefix mix1
        if( length(teta1) != 0) der2teta.derteta.stteta.st.mix1[teta.ind1] <- BITS1.mix1$der2teta.derteta.stteta.st
        #New-theta1
        # prefix mix2
        if( length(teta1) != 0) der2teta.derteta.stteta.st.mix2[teta.ind1] <- BITS1.mix2$der2teta.derteta.stteta.st
        #New-theta1
        # prefix .2
        if( length(teta1) != 0) der2teta.derteta.stteta.st.2[teta.ind1] <- BITS1.2$der2teta.derteta.stteta.st
        
        if( length(teta2) != 0) der2teta.derteta.stteta.st[teta.ind2] <- BITS2$der2teta.derteta.stteta.st 
        #New-theta2
        #prefix.mix1
        if( length(teta2) != 0) der2teta.derteta.stteta.st.mix1[teta.ind2] <- BITS2.mix1$der2teta.derteta.stteta.st 
        #New-theta2
        #prefix .mix2
        if( length(teta2) != 0) der2teta.derteta.stteta.st.mix2[teta.ind2] <- BITS2.mix2$der2teta.derteta.stteta.st 
        #New-theta2
        #prefix .2
        if( length(teta2) != 0) der2teta.derteta.stteta.st.2[teta.ind2] <- BITS2.2$der2teta.derteta.stteta.st 
        ########################
        
        der3C.derp1p1p1 <- der3C.derp1tetateta <- der2h.derteta.teta.st <- der3C.p1p1teta <- der2h.derp2teta <- der2h.derp2p2 <- NA
        
        der3C.derp1p1p1.mix1 <- der3C.derp1tetateta.mix1 <- der2h.derteta.teta.st.mix1 <- der3C.p1p1teta.mix1 <- der2h.derp2teta.mix1 <- der2h.derp2p2.mix1 <- NA
        
        der3C.derp1p1p1.mix2 <- der3C.derp1tetateta.mix2 <- der2h.derteta.teta.st.mix2 <- der3C.p1p1teta.mix2 <- der2h.derp2teta.mix2 <- der2h.derp2p2.mix2 <- NA
        
        der3C.derp1p1p1.2 <- der3C.derp1tetateta.2 <- der2h.derteta.teta.st.2 <- der3C.p1p1teta.2 <- der2h.derp2teta.2 <- der2h.derp2p2.2 <- NA
        
        if( length(teta1) != 0){der3C.derp1p1p1[teta.ind1]       <- BITS1$der3C.derp1p1p1 
        der2h.derteta.teta.st[teta.ind1] <- BITS1$der2h.derteta.teta.st
        der3C.derp1tetateta[teta.ind1]   <- BITS1$der3C.derp1tetateta
        der3C.p1p1teta[teta.ind1]        <- BITS1$der3C.p1p1teta  
        der2h.derp2teta[teta.ind1]       <- BITS1$der2h.derp2teta
        der2h.derp2p2[teta.ind1]         <- BITS1$der2h.derp2p2
        der2h.derp1teta[teta.ind1]       <- BITS1$der2h.derp1teta
        }
        
        #New- theta 1
        #prefix mix1
        
        if( length(teta1) != 0){
                der3C.derp1p1p1.mix1[teta.ind1]       <- BITS1.mix1$der3C.derp1p1p1 
                der2h.derteta.teta.st.mix1[teta.ind1] <- BITS1.mix1$der2h.derteta.teta.st
                der3C.derp1tetateta.mix1[teta.ind1]   <- BITS1.mix1$der3C.derp1tetateta
                der3C.p1p1teta.mix1[teta.ind1]        <- BITS1.mix1$der3C.p1p1teta  
                der2h.derp2teta.mix1[teta.ind1]       <- BITS1.mix1$der2h.derp2teta
                der2h.derp2p2.mix1[teta.ind1]         <- BITS1.mix1$der2h.derp2p2
                der2h.derp1teta.mix1[teta.ind1]       <- BITS1.mix1$der2h.derp1teta
        }
        
        #New-theta1
        #prefix mix2
        
        if( length(teta1) != 0){
                der3C.derp1p1p1.mix2[teta.ind1]       <- BITS1.mix2$der3C.derp1p1p1 
                der2h.derteta.teta.st.mix2[teta.ind1] <- BITS1.mix2$der2h.derteta.teta.st
                der3C.derp1tetateta.mix2[teta.ind1]   <- BITS1.mix2$der3C.derp1tetateta
                der3C.p1p1teta.mix2[teta.ind1]        <- BITS1.mix2$der3C.p1p1teta  
                der2h.derp2teta.mix2[teta.ind1]       <- BITS1.mix2$der2h.derp2teta
                der2h.derp2p2.mix2[teta.ind1]         <- BITS1.mix2$der2h.derp2p2
                der2h.derp1teta.mix2[teta.ind1]       <- BITS1.mix2$der2h.derp1teta
        }
        
        #New-theta1
        #prefix .2
        
        if( length(teta1) != 0){
                der3C.derp1p1p1.2[teta.ind1]       <- BITS1.2$der3C.derp1p1p1 
                der2h.derteta.teta.st.2[teta.ind1] <- BITS1.2$der2h.derteta.teta.st
                der3C.derp1tetateta.2[teta.ind1]   <- BITS1.2$der3C.derp1tetateta
                der3C.p1p1teta.2[teta.ind1]        <- BITS1.2$der3C.p1p1teta  
                der2h.derp2teta.2[teta.ind1]       <- BITS1.2$der2h.derp2teta
                der2h.derp2p2.2[teta.ind1]         <- BITS1.2$der2h.derp2p2
                der2h.derp1teta.2[teta.ind1]       <- BITS1.2$der2h.derp1teta
        }
        
        if( length(teta2) != 0){
                der3C.derp1p1p1[teta.ind2]       <- BITS2$der3C.derp1p1p1
                der2h.derteta.teta.st[teta.ind2] <- BITS2$der2h.derteta.teta.st
                der3C.derp1tetateta[teta.ind2]   <- BITS2$der3C.derp1tetateta
                der3C.p1p1teta[teta.ind2]        <- BITS2$der3C.p1p1teta
                der2h.derp2teta[teta.ind2]       <- BITS2$der2h.derp2teta
                der2h.derp2p2[teta.ind2]         <- BITS2$der2h.derp2p2
                der2h.derp1teta[teta.ind2]       <- BITS2$der2h.derp1teta
        }
        
        #New-theta2
        #prefix .mix1
        
        if( length(teta2) != 0){
                der3C.derp1p1p1.mix1[teta.ind2]       <- BITS2.mix1$der3C.derp1p1p1
                der2h.derteta.teta.st.mix1[teta.ind2] <- BITS2.mix1$der2h.derteta.teta.st
                der3C.derp1tetateta.mix1[teta.ind2]   <- BITS2.mix1$der3C.derp1tetateta
                der3C.p1p1teta.mix1[teta.ind2]        <- BITS2.mix1$der3C.p1p1teta
                der2h.derp2teta.mix1[teta.ind2]       <- BITS2.mix1$der2h.derp2teta
                der2h.derp2p2.mix1[teta.ind2]         <- BITS2.mix1$der2h.derp2p2
                der2h.derp1teta.mix1[teta.ind2]       <- BITS2.mix1$der2h.derp1teta
        }
        
        #New-theta2
        #prefix .mix2
        
        if( length(teta2) != 0){
                der3C.derp1p1p1.mix2[teta.ind2]       <- BITS2.mix2$der3C.derp1p1p1
                der2h.derteta.teta.st.mix2[teta.ind2] <- BITS2.mix2$der2h.derteta.teta.st
                der3C.derp1tetateta.mix2[teta.ind2]   <- BITS2.mix2$der3C.derp1tetateta
                der3C.p1p1teta.mix2[teta.ind2]        <- BITS2.mix2$der3C.p1p1teta
                der2h.derp2teta.mix2[teta.ind2]       <- BITS2.mix2$der2h.derp2teta
                der2h.derp2p2.mix2[teta.ind2]         <- BITS2.mix2$der2h.derp2p2
                der2h.derp1teta.mix2[teta.ind2]       <- BITS2.mix2$der2h.derp1teta
        }
        
        #New-theta2
        #prefix .2
        
        if( length(teta2) != 0){
                der3C.derp1p1p1.2[teta.ind2]       <- BITS2.2$der3C.derp1p1p1
                der2h.derteta.teta.st.2[teta.ind2] <- BITS2.2$der2h.derteta.teta.st
                der3C.derp1tetateta.2[teta.ind2]   <- BITS2.2$der3C.derp1tetateta
                der3C.p1p1teta.2[teta.ind2]        <- BITS2.2$der3C.p1p1teta
                der2h.derp2teta.2[teta.ind2]       <- BITS2.2$der2h.derp2teta
                der2h.derp2p2.2[teta.ind2]         <- BITS2.2$der2h.derp2p2
                der2h.derp1teta.2[teta.ind2]       <- BITS2.2$der2h.derp1teta
        }
        
        ########################
        ########################                   
        
        #Likelihood, Gradient and Hessian 



#initialization 
likelihood <- 0
G <- 0
H <- 0
#UU
if(sum(VC$indUU)>1){
        #Likelihood
        l.par <- VC$weights*( VC$indUU*( log(c.copula2.be1be2) + log(-dS1eta1) + log(-dS2eta2) + log(Xd1P) + log(Xd2P) ))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indUU*(c(c.copula2.be1be2^(-1)*der2h.derp1p1*dS1eta1) * dereta1derb1+
                                                  c(dS1eta1^(-1)*d2S1eta1)*dereta1derb1 
                                                  +c(Xd1P)^(-1)* der2eta1dery1b1
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indUU*(c(c.copula2.be1be2^(-1)*der2h.derp1p2*dS2eta2)*dereta2derb2+
                                                  c((dS2eta2)^(-1)*d2S2eta2)*dereta2derb2
                                                  +c(Xd2P)^(-1)*der2eta2dery2b2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indUU*(c.copula2.be1be2^(-1)*der2h.derp1teta.st
        ))*X3
        
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        #Hessian
        
        be1.be1 <-  -( 
                
                crossprod(VC$weights*VC$indUU*c(-c.copula2.be1be2^-2*der2h.derp1p1^2*dS1eta1^2 
                                                + c.copula2.be1be2^-1*der2c.derp1.derp1*dS1eta1^2 
                                                + c.copula2.be1be2^-1*der2h.derp1p1*d2S1eta1 -dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUU*c(c.copula2.be1be2^-1*der2h.derp1p1*dS1eta1 
                                                                  + dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) +
                        
                        crossprod(VC$weights*VC$indUU*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUU*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )      )
        
        be2.be2 <-  -( 
                
                crossprod(VC$weights*VC$indUU*c(-c.copula2.be1be2^-2*der2h.derp1p2^2*dS2eta2^2 
                                                + c.copula2.be1be2^-1*der2c.derp2.derp2*dS2eta2^2 
                                                + c.copula2.be1be2^-1*der2h.derp1p2*d2S2eta2 -dS2eta2^-2*d2S2eta2^2 
                                                + dS2eta2^-1*d3S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUU*c(c.copula2.be1be2^-1*der2h.derp1p2*dS2eta2 + dS2eta2^-1*d2S2eta2)*VC$X2)*der2.par2 ) ) ) +
                        
                        crossprod(VC$weights*VC$indUU*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUU*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) ) 
        )
                
        be1.be2 <- -( 
                
                crossprod(VC$weights*VC$indUU*c((-c.copula2.be1be2^-2*der2h.derp1p2*der2h.derp1p1 
                                                 + c.copula2.be1be2^-1*der2c.derp1.derp2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
                
                )
        

        d2l.rho.rho <- -( 
                
                VC$weights*VC$indUU*( -c.copula2.be1be2^-2*der2h.derp1teta^2*derteta.derteta.st^2 
                                      + c.copula2.be1be2^-1*der2c.derrho.derrho*derteta.derteta.st^2 
                                      + c.copula2.be1be2^-1*der2h.derp1teta*der2teta.derteta.stteta.st)
        )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        
        be1.rho <- -( 
                
                crossprod(VC$weights*VC$indUU*c((-c.copula2.be1be2^-2*der2h.derp1p1*der2h.derp1teta 
                                                 + c.copula2.be1be2^-1*der2c.derp1.derrho)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3) 
        )
        be2.rho <- -( 
                
                crossprod(VC$weights*VC$indUU*c((-c.copula2.be1be2^-2*der2h.derp1p2*der2h.derp1teta 
                                                 + c.copula2.be1be2^-1*der2c.derp2.derrho)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3) 
        )
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                    cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                    cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
}


if(sum(VC$indRR)>1){
        #Likelihood
        l.par <- VC$weights*( VC$indRR*log(mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)) )
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indRR*(c(mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula.be1*dS1eta1) *dereta1derb1))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*( VC$indRR*(c(mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(c.copula.be2*dS2eta2))*dereta2derb2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indRR*(mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula.theta
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indRR*c(-mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.be1^2*dS1eta1^2 
                                                + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1*dS1eta1^2  
                                                + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.be1*d2S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRR*c( mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.be1*dS1eta1  )*VC$X1)*der2.par1 ) ) )
        )
        be2.be2 <-  -( 
                
                
                crossprod(VC$weights*VC$indRR*c(-mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.be2^2*dS2eta2^2 
                                                + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be2*dS2eta2^2 
                                                + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.be2*d2S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRR*c( mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.be2*dS2eta2  )*VC$X2)*der2.par2 ) ) )
        )
        
        be1.be2 <- -(
               
                crossprod(VC$weights*VC$indRR*c((-mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.be2*c.copula.be1 
                                                 + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        
        d2l.rho.rho <- -( 
                    
                VC$weights*VC$indRR*( -mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.thet^2*derteta.derteta.st^2 
                                      + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*bit1.th2ATE*derteta.derteta.st^2 
                                      + rotConst*mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.thet*der2teta.derteta.stteta.st )
                
        )
        
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                
                crossprod(VC$weights*VC$indRR*c(rotConst*(-mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.be1*c.copula.thet 
                                                          + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1t)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)
        )
        be2.rho <- -( 
                
                crossprod(VC$weights*VC$indRR*c(rotConst*(-mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.be2*c.copula.thet 
                                                          + mm(p00,min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be2t)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
        )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
        
}

#LL
#mm() has been added according to what has been discussed in the document
if(sum(VC$indLL)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indLL*log(mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indLL*(c(mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*c(((-dS1eta1)+c.copula.be1*dS1eta1))*dereta1derb1
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indLL*(c(mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*((-dS2eta2)+c.copula.be2*dS2eta2))*dereta2derb2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indLL*(c(mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula.theta)
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indLL*c(-mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-dS1eta1+c.copula.be1*dS1eta1)^2 + 
                                                        mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1*dS1eta1^2 + c.copula.be1*d2S1eta1-d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLL*c( mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be1*dS1eta1 -dS1eta1) )*VC$X1)*der2.par1 ) ) )
        )
        be2.be2 <-  -( 
        
                crossprod(VC$weights*VC$indLL*c(-mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-dS2eta2+c.copula.be2*dS2eta2)^2 
                                                + mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2*dS2eta2^2 + c.copula.be2*d2S2eta2-d2S2eta2) )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLL*c( mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2*dS2eta2 - dS2eta2)  )*VC$X2)*der2.par2 ) ) )
        )
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indLL*c((-mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*((c.copula.be2-1)*(c.copula.be1-1)) 
                                                 + mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2)*dS2eta2*dS1eta1)*dereta1derb1, dereta2derb2)
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -(
                
                VC$weights*VC$indLL*( -mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula.thet^2*derteta.derteta.st^2 
                                      + mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*bit1.th2ATE*derteta.derteta.st^2 
                                      + rotConst*mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.thet*der2teta.derteta.stteta.st )
        )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                
                crossprod(VC$weights*VC$indLL*c(rotConst*(-mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1-1)*c.copula.thet 
                                                          + mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1t)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)
        )
        
        be2.rho <- -( 
                
                crossprod(VC$weights*VC$indLL*c(rotConst*(-mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2-1)*c.copula.thet 
                                                          + mm(1-p1-p2+p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be2t)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
                
                )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
        
        
        
}


#UR
#mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)
#mm() added only to c.copula.be1
if(sum(VC$indUR)>1){
        
        #Likelihood
        l.par <- VC$weights*( VC$indUR*(log(mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr))+log(-dS1eta1)+log(Xd1P)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        
        dl.dbe1 <- -VC$weights*( VC$indUR*(c(mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula2.be1*dS1eta1)*dereta1derb1+
                                                   c((dS1eta1)^(-1)*(d2S1eta1)) *dereta1derb1
                                           +c(Xd1P)^(-1)*der2eta1dery1b1
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indUR*(c(mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula2.be1be2*dS2eta2)*dereta2derb2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indUR*(mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula2.be1th
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                    
                crossprod(VC$weights*VC$indUR*c(-mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1^2*dS1eta1^2  
                                                + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der3C.derp1p1p1*dS1eta1^2 
                                                + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1*d2S1eta1 -dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUR*c(mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1*dS1eta1  + dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) +
                        
                        crossprod(VC$weights*VC$indUR*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUR*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )     
                
        )
        be2.be2 <-  -(
                
                crossprod(VC$weights*VC$indUR*c(-mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2^2*dS2eta2^2 
                                                + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp1p2*dS2eta2^2 
                                                + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*d2S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUR*c( mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*dS2eta2  )*VC$X2)*der2.par2 ) ) ) 
        )
        be1.be2 <- -(
                  
                crossprod(VC$weights*VC$indUR*c((-mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2*c.copula2.be1 + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp1p1)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -(
                      
                VC$weights*VC$indUR*( -mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1t^2*derteta.derteta.st^2  
                                      + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der3C.derp1tetateta*derteta.derteta.st^2 
                                      + rotConst*mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1t*der2teta.derteta.stteta.st ) 
                
        )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                       
                crossprod(VC$weights*VC$indUR*c((rotConst*-mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1*c.copula2.be1t 
                                                 + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der3C.p1p1teta)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)  
        )
        be2.rho <- -( 
                    
                crossprod(VC$weights*VC$indUR*c((rotConst*-mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2*c.copula2.be1t 
                                                 + mm(c.copula.be1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp1teta)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
        )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
        
}


#RU
# mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)
#mm() added only to c.copula.be2
if(sum(VC$indRU)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indRU*(log(mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr))+ log(-dS2eta2)+ log(Xd2P)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indRU*(c(mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula2.be1be2*dS1eta1)*dereta1derb1))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indRU*(c(mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula2.be2*dS2eta2)*dereta2derb2+
                                                  c((dS2eta2)^(-1)*(d2S2eta2))*dereta2derb2
                                          +c(Xd2P)^(-1)*der2eta2dery2b2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indRU*(mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*c.copula2.be2th
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indRU*c(-mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2^2*dS1eta1^2 
                                                + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp1p1*dS1eta1^2 
                                                + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*d2S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRU*c( mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*dS1eta1  )*VC$X1)*der2.par1 ) ) )
        )
        
        be2.be2 <-  -( 
               
                
                crossprod(VC$weights*VC$indRU*c(-mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be2^2*dS2eta2^2 
                                                + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp2p2*dS2eta2^2 
                                                + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be2*d2S2eta2 
                                                -dS2eta2^-2*d2S2eta2^2 + dS2eta2^-1*d3S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRU*c(mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be2*dS2eta2 
                                                                  + dS2eta2^-1*d2S2eta2)*VC$X2)*der2.par2 ) ) ) +
                        
                        crossprod(VC$weights*VC$indRU*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRU*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) ) 
        )
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indRU*c((-mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2*c.copula2.be2  
                                                 + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp1p2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -( 
                     
                VC$weights*VC$indRU*( -mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be2t^2*derteta.derteta.st^2  
                                      + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derteta.teta.st*derteta.derteta.st^2 
                                      + rotConst*mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be2t*der2teta.derteta.stteta.st )
        
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        
        be1.rho <- -(
                   
                crossprod(VC$weights*VC$indRU*c((rotConst*-mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2*c.copula2.be2t  
                                                 + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp1teta)*dS1eta1*derteta.derteta.st)*dereta1derb1, X3) 
        )
        
        be2.rho <- -( 
                   
                crossprod(VC$weights*VC$indRU*c((rotConst*-mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be2*c.copula2.be2t  
                                                 + mm(c.copula.be2 , min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*der2h.derp2teta)*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
        )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
}



#UL
if(sum(VC$indUL)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indUL*(log( (c.copula.be1-1) * (dS1eta1) * Xd1P)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #grad
        dl.dbe1 <- -VC$weights*(VC$indUL*(c((c.copula.be1-1)^(-1)*c.copula2.be1*dS1eta1)*dereta1derb1+
                                                  c((dS1eta1^(-1))*d2S1eta1)*dereta1derb1
                                          +c(Xd1P)^(-1)*der2eta1dery1b1
        ))
        
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*( VC$indUL*(c((c.copula.be1-1)^(-1)*c.copula2.be1be2*dS2eta2)*dereta2derb2
        ))
        
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indUL*((c.copula.be1-1)^(-1)*c.copula2.be1th
        ))*X3
        
        dl.dteta.st <- colSums( dl.dteta.st)
        
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indUL*c(-(c.copula.be1-1)^-2*c.copula2.be1^2*dS1eta1^2 + (c.copula.be1-1)^-1*der3C.derp1p1p1*dS1eta1^2 
                                                + (c.copula.be1-1)^-1*c.copula2.be1*d2S1eta1 
                                                -dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUL*c((c.copula.be1-1)^-1*c.copula2.be1*dS1eta1 + dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) +
                        
                        crossprod(VC$weights*VC$indUL*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUL*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )
                
        )
        be2.be2 <-  -( 
                
                crossprod(VC$weights*VC$indUL*c(-(c.copula.be1-1)^-2*c.copula2.be1be2^2*dS2eta2^2 
                                                + (c.copula.be1-1)^-1*der2h.derp1p2*dS2eta2^2 
                                                + (c.copula.be1-1)^-1*c.copula2.be1be2*d2S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUL*c( (c.copula.be1-1)^-1*c.copula2.be1be2*dS2eta2  )*VC$X2)*der2.par2 ) ) )
                
                )
        
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indUL*c((-(c.copula.be1-1)^-2*c.copula2.be1*c.copula2.be1be2 + (c.copula.be1-1)^-1*der2h.derp1p1)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
        )
        
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        
        d2l.rho.rho <- -( 
               
                VC$weights*VC$indUL*( -(c.copula.be1-1)^-2*c.copula2.be1t^2*derteta.derteta.st^2 
                                      + (c.copula.be1-1)^-1*der3C.derp1tetateta*derteta.derteta.st^2 
                                      + rotConst*(c.copula.be1-1)^-1*c.copula2.be1t*der2teta.derteta.stteta.st )
                
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        
        
    be1.rho <- -(crossprod(VC$weights * VC$indUL * c( (rotConst *-(c.copula.be1 - 1)^-2 * c.copula2.be1 * c.copula2.be1t + 
                                                          (c.copula.be1 - 1)^-1 * der3C.p1p1teta) * dS1eta1 * 
                                                       derteta.derteta.st) * dereta1derb1, X3))
    be2.rho <- -(crossprod(VC$weights * VC$indUL * c( (rotConst *-(c.copula.be1 - 1)^-2 * (c.copula2.be1be2) * c.copula2.be1t + 
                                                          (c.copula.be1 - 1)^-1 * der2h.derp1teta) * dS2eta2 * 
                                                       derteta.derteta.st) * dereta2derb2, X3))        
        
        
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
        
}
#LU  
if(sum(VC$indLU)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indLU*(log( (c.copula.be2-1) * (dS2eta2) * Xd2P)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indLU*(c((c.copula.be2-1)^(-1)*c.copula2.be1be2*dS1eta1)*dereta1derb1))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indLU*(c((c.copula.be2-1)^(-1)*c.copula2.be2*dS2eta2)*dereta2derb2+
                                                  c((dS2eta2)^(-1)*d2S2eta2)*dereta2derb2
                                          +c(Xd2P)^(-1)*der2eta2dery2b2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*( VC$indLU*((c.copula.be2-1)^(-1)*c.copula2.be2th
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                
                crossprod(VC$weights*VC$indLU*c(-(c.copula.be2-1)^-2*c.copula2.be1be2^2*dS1eta1^2 
                                                + (c.copula.be2-1)^-1*der2h.derp1p1*dS1eta1^2 
                                                + (c.copula.be2-1)^-1*c.copula2.be1be2*d2S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLU*c( (c.copula.be2-1)^-1*c.copula2.be1be2*dS1eta1  )*VC$X1)*der2.par1 ) ) )
        )
        be2.be2 <-  -(
                
                
                crossprod(VC$weights*VC$indLU*c(-(c.copula.be2-1)^-2*c.copula2.be2^2*dS2eta2^2 
                                                + (c.copula.be2-1)^-1*der2h.derp2p2*dS2eta2^2 
                                                + (c.copula.be2-1)^-1*c.copula2.be2*d2S2eta2 
                                                -dS2eta2^-2*d2S2eta2^2 + dS2eta2^-1*d3S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLU*c((c.copula.be2-1)^-1*c.copula2.be2*dS2eta2 
                                                                  + dS2eta2^-1*d2S2eta2)*VC$X2)*der2.par2 ) ) ) +
                        
                        crossprod(VC$weights*VC$indLU*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLU*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) )
        )
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indLU*c((-(c.copula.be2-1)^-2*c.copula2.be1be2*c.copula2.be2 
                                                 + (c.copula.be2-1)^-1*der2h.derp1p2)*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
        )
        
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        
        d2l.rho.rho <- -( 
                
                VC$weights*VC$indLU*( -(c.copula.be2-1)^-2*c.copula2.be2t^2*derteta.derteta.st^2 
                                      + (c.copula.be2-1)^-1*der2h.derteta.teta.st*derteta.derteta.st^2 
                                      + rotConst*(c.copula.be2-1)^-1*c.copula2.be2t*der2teta.derteta.stteta.st )
                )
        
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        
    be1.rho <- -(crossprod(VC$weights * VC$indLU * c( 
                                                       (rotConst *-(c.copula.be2 - 1)^-2 * c.copula2.be1be2 * c.copula2.be2t + 
                                                          (c.copula.be2 - 1)^-1 * der2h.derp1teta) * dS1eta1 * 
                                                       derteta.derteta.st) * dereta1derb1, X3))
    be2.rho <- -(crossprod(VC$weights * VC$indLU * c(
                                                       (rotConst * -(c.copula.be2 - 1)^-2 * (c.copula2.be2) * c.copula2.be2t + 
                                                          (c.copula.be2 - 1)^-1 * der2h.derp2teta) * dS2eta2 * 
                                                       derteta.derteta.st) * dereta2derb2, X3))
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
       
}
#RL

# #mm() has been added according to what has been discussed in the document
if(sum(VC$indRL)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indRL*log(mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        dl.dbe1 <- -VC$weights*( VC$indRL*(c( (mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)) *(
                dS1eta1
                -c.copula.be1*dS1eta1))*dereta1derb1
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*( VC$indRL*(c(mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula.be2*dS2eta2))*dereta2derb2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*( VC$indRL*(c(mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula.theta))
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indRL*c(-mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(dS1eta1-c.copula.be1*dS1eta1)^2 
                                                + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(d2S1eta1-c.copula2.be1*dS1eta1^2-c.copula.be1*d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRL*c( mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be1*dS1eta1+dS1eta1)  )*VC$X1)*der2.par1 ) ) )
        )
        be2.be2 <-  -( 
                
                crossprod(VC$weights*VC$indRL*c(-mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2)^2*dS2eta2^2 
                                                + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2)*dS2eta2^2 
                                                + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be2)*d2S2eta2)*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRL*c( mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be2)*dS2eta2  )*VC$X2)*der2.par2 ) ) )
        )
        be1.be2 <- -(
                
                
                crossprod(VC$weights*VC$indRL*c((-mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be1)*(-c.copula.be2) 
                                                 + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1be2))*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
        
                )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -(
                
                VC$weights*VC$indRL*( -mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.thet)^2*derteta.derteta.st^2 
                                      + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-bit1.th2ATE)*derteta.derteta.st^2 
                                      + rotConst*mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.thet)*der2teta.derteta.stteta.st )
        
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                
                crossprod(VC$weights*VC$indRL*c(rotConst*(-mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be1)*(-c.copula.thet) 
                                                          + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1t))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)
        )
        be2.rho <- -( 
                
                crossprod(VC$weights*VC$indRL*c(rotConst*(-mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2)*(-c.copula.thet) 
                                                          + mm(p1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2t))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
        
                )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
}


#LR

#mm() has been added according to what has been discussed in the document
if(sum(VC$indLR)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indLR*log(mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        
        dl.dbe1 <- -VC$weights*( VC$indLR*(c((mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(-c.copula.be1*dS1eta1)) * dereta1derb1
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indLR*(c(mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(dS2eta2)*dereta2derb2
                -c(c.copula.be2*dS2eta2)*dereta2derb2)
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*( VC$indLR*(c(mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula.theta))
        ))*X3
        
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indLR*c(-mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1*dS1eta1)^2 
                                                + mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1*dS1eta1^2-c.copula.be1*d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLR*c( mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be1*dS1eta1)  )*VC$X1)*der2.par1 ) ) )
                
        )
        be2.be2 <-  -(
        
                crossprod(VC$weights*VC$indLR*c(-mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(dS2eta2-c.copula.be2*dS2eta2)^2 
                                                + mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2*dS2eta2^2 -c.copula.be2*d2S2eta2+ d2S2eta2) )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLR*c( mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be2*dS2eta2 + dS2eta2)  )*VC$X2)*der2.par2 ) ) )
        )
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indLR*c((-mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be2)*(-c.copula.be1) 
                                                 + mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1be2))*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -( 
        
                VC$weights*VC$indLR*( -mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.thet)^2*derteta.derteta.st^2 
                                      + mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-bit1.th2ATE)*derteta.derteta.st^2 
                                      + rotConst*mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.thet)*der2teta.derteta.stteta.st )
        )
        
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                
                crossprod(VC$weights*VC$indLR*c(rotConst*(-mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1)*(-c.copula.thet) 
                                                          + mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1t))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)
        )
        
        be2.rho <- -( 
                
                crossprod(VC$weights*VC$indLR*c(rotConst*(-mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be2)*(-c.copula.thet) 
                                                          + mm(p2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2t))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
        
                )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
}

#RI
if(sum(VC$indRI)>1){
        
        #Likelihood
        l.par <- VC$weights*(VC$indRI*log(mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        
        dl.dbe1 <- -VC$weights*(VC$indRI*(mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(
                c(c.copula.be1*dS1eta1)*dereta1derb1
                -c(c.copula.be1.mix1*dS1eta1)*dereta1derb1)
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indRI*(c(mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula.be2*dS2eta2)*dereta2derb2
                -c(c.copula.be2.mix1*dS2eta2.2)*dereta2derb2.2)
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indRI*(c(mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c.copula.theta
                -c.copula.theta.mix1)
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indRI*c(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1*dS1eta1-c.copula.be1.mix1*dS1eta1)^2 
                                                + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1*dS1eta1^2+c.copula.be1*d2S1eta1
                                                                                                               -c.copula2.be1.mix1*dS1eta1^2-c.copula.be1.mix1*d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRI*c( mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be1*dS1eta1-c.copula.be1.mix1*dS1eta1)  )*VC$X1)*der2.par1 ) ) )
                
        )
        be2.be2 <-  -(
                
                
                crossprod(VC$weights*VC$indRI*c(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2*dS2eta2)^2 
                                                + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2*dS2eta2^2 + 
                                                                                                                       c.copula.be2*d2S2eta2) )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRI*c( mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula.be2*dS2eta2 )*VC$X2)*der2.par2 ) ) )+
                       
                        crossprod(VC$weights*VC$indRI*c(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2.mix1*dS2eta2.2)^2 
                                                        + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2.mix1*dS2eta2.2^2  
                                                                                                                       -c.copula.be2.mix1*d2S2eta2.2) )*dereta2derb2.2, dereta2derb2.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indRI*c( mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be2.mix1)*dS2eta2.2
                        )*VC$X2.2)*der2.par2 ) ) )  
                +
                        
                       
                        crossprod(VC$weights*VC$indRI*c(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2.mix1*c.copula.be2)*dS2eta2*dS2eta2.2 )*dereta2derb2, dereta2derb2.2)+
                        crossprod(VC$weights*VC$indRI*c(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2.mix1*c.copula.be2)*dS2eta2*dS2eta2.2 )*dereta2derb2.2, dereta2derb2)
                
        )
        
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indRI*c((-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2)*(c.copula.be1-c.copula.be1.mix1)+
                                                         mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(c.copula2.be1be2)
                )*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indRI*c((-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2.mix1)*(c.copula.be1-c.copula.be1.mix1) +
                                                        mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula2.be1be2.mix1))*dS1eta1*dS2eta2.2)*dereta1derb1, dereta2derb2.2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -(
                
                
                VC$weights*VC$indRI*( -mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.thet*derteta.derteta.st-c.copula.thet.mix1*derteta.derteta.st)^2 
                                      + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(bit1.th2ATE-bit1.th2ATE.mix1)*derteta.derteta.st^2 
                                      + rotConst*mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.thet-c.copula.thet.mix1)*der2teta.derteta.stteta.st )
        )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                crossprod(VC$weights*VC$indRI*c(rotConst*(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1-c.copula.be1.mix1)*(c.copula.thet-c.copula.thet.mix1) 
                                                          + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t-c.copula2.be1t.mix1))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)
                )
        be2.rho <- -(
                
                crossprod(VC$weights*VC$indRI*c(rotConst*(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2)*(c.copula.thet-c.copula.thet.mix1) 
                                                          + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)+
                        crossprod(VC$weights*VC$indRI*c(rotConst*(-mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be2.mix1)*(c.copula.thet-c.copula.thet.mix1) 
                                                                  + mm(p00-p00.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2t.mix1))*dS2eta2.2*derteta.derteta.st)*dereta2derb2.2, X3)
        )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
}
#IR
if(sum(VC$indIR)>1){
        #Likelihood
        l.par <- VC$weights*( VC$indIR*log(mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)) )
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indIR*(mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(
                c(c.copula.be1*dS1eta1)*dereta1derb1
                -c(c.copula.be1.mix2*dS1eta1.2)*dereta1derb1.2)
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*( VC$indIR*(c(mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula.be2*dS2eta2)*dereta2derb2
                -c(c.copula.be2.mix2*dS2eta2)*dereta2derb2)
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indIR*(c(mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c.copula.theta
                -c.copula.theta.mix2)
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                
                crossprod(VC$weights*VC$indIR*c(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1*dS1eta1)^2 
                                                + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1*dS1eta1^2+c.copula.be1*d2S1eta1
                                                ))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIR*c( mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be1*dS1eta1)  )*VC$X1)*der2.par1 ) ) )+
                        
                        
                crossprod(VC$weights*VC$indIR*c(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1.mix2*dS1eta1.2)^2 
                                                        + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1.mix2*dS1eta1.2^2-c.copula.be1.mix2*d2S1eta1.2
                                                        ))*dereta1derb1.2, dereta1derb1.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIR*c( mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula.be1.mix2*dS1eta1.2)  )*VC$X1.2)*der2.par1 ) ) )+
                        
                crossprod(VC$weights*VC$indIR*c(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1.mix2*dS1eta1.2)*(c.copula.be1*dS1eta1) 
                        )*dereta1derb1, dereta1derb1.2) +
                crossprod(VC$weights*VC$indIR*c(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1.mix2*dS1eta1.2)*(c.copula.be1*dS1eta1) 
                        )*dereta1derb1.2, dereta1derb1)
        )
        
        be2.be2 <-  -( 
                
                crossprod(VC$weights*VC$indIR*c(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2*dS2eta2-c.copula.be2.mix2*dS2eta2)^2 
                                                + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2-c.copula2.be2.mix2)*dS2eta2^2 
                                                + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2-c.copula.be2.mix2)*d2S2eta2 )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIR*c( mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2-c.copula.be2.mix2)*dS2eta2 
                        )*VC$X2)*der2.par2 ) ) )
                )
        
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indIR*c((-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1)*(c.copula.be2-c.copula.be2.mix2)
                                                 +mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(c.copula2.be1be2) )*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indIR*c((-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1.mix2)*(c.copula.be2-c.copula.be2.mix2)
                                                         +mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula2.be1be2.mix2) )*dS1eta1.2*dS2eta2)*dereta1derb1.2, dereta2derb2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        
        d2l.rho.rho <- -( 
                
                VC$weights*VC$indIR*( -mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.thet*derteta.derteta.st-c.copula.thet.mix2*derteta.derteta.st)^2 
                                      + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(bit1.th2ATE-bit1.th2ATE.mix2)*derteta.derteta.st^2 
                                      + rotConst*mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.thet-c.copula.thet.mix2)*der2teta.derteta.stteta.st )
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -(
                
                crossprod(VC$weights*VC$indIR*c(rotConst*(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1)*(c.copula.thet-c.copula.thet.mix2) 
                                                          + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)+
                        
                        crossprod(VC$weights*VC$indIR*c(rotConst*(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1.mix2)*(c.copula.thet-c.copula.thet.mix2)
                                                                  + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1t.mix2))*dS1eta1.2*derteta.derteta.st)*dereta1derb1.2, X3)
        )
        be2.rho <- -(
                
                
                crossprod(VC$weights*VC$indIR*c(rotConst*(-mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2-c.copula.be2.mix2)*(c.copula.thet-c.copula.thet.mix2) 
                                                          + mm(p00-p00.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t-c.copula2.be2t.mix2))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
                
        )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
        
}

#LI
if(sum(VC$indLI)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indLI*log(mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        #Gradient
        dl.dbe1 <- -VC$weights*( VC$indLI*(c(mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula.be1.mix1*dS1eta1
                  -c.copula.be1*dS1eta1)*dereta1derb1)
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indLI*(c(mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(dS2eta2)*dereta2derb2
                -c(dS2eta2.2)*dereta2derb2.2
                +c(c.copula.be2.mix1*dS2eta2.2)*dereta2derb2.2
                -c(c.copula.be2*dS2eta2)*dereta2derb2)
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indLI*(c(mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(
                c.copula.theta.mix1
                -c.copula.theta))
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indLI*c(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.mix1*dS1eta1-c.copula.be1*dS1eta1)^2 
                                                + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1.mix1*dS1eta1^2+c.copula.be1.mix1*d2S1eta1
                                                                                                                       -c.copula2.be1*dS1eta1^2-c.copula.be1*d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLI*c( mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be1.mix1*dS1eta1-c.copula.be1*dS1eta1)  )*VC$X1)*der2.par1 ) ) )
        )
        
        be2.be2 <-  -( 
                
                crossprod(VC$weights*VC$indLI*c(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(dS2eta2-c.copula.be2*dS2eta2)^2 
                                                + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2)*dS2eta2^2 
                                                + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(1-c.copula.be2)*d2S2eta2 )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLI*c( mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(1-c.copula.be2)*dS2eta2 
                        )*VC$X2)*der2.par2 ) ) )+
                        
                        crossprod(VC$weights*VC$indLI*c(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.mix1*dS2eta2.2-dS2eta2.2)^2
                                                        + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2.mix1)*dS2eta2.2^2 
                                                        +  mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2.mix1-1)*d2S2eta2.2 )*dereta2derb2.2, dereta2derb2.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indLI*c(mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2.mix1-1)*dS2eta2.2
                        )*VC$X2.2)*der2.par2 ) ) )+ 
                        
                        
                        crossprod(VC$weights*VC$indLI*c(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*((dS2eta2-c.copula.be2*dS2eta2)*(c.copula.be2.mix1*dS2eta2.2-dS2eta2.2)) )*dereta2derb2, dereta2derb2.2)+
                        crossprod(VC$weights*VC$indLI*c(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*((dS2eta2-c.copula.be2*dS2eta2)*(c.copula.be2.mix1*dS2eta2.2-dS2eta2.2)) )*dereta2derb2.2, dereta2derb2)
                )
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indLI*c((-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.mix1-c.copula.be1)*(1-c.copula.be2)
                                                 +mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula2.be1be2) )*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indLI*c((-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.mix1-c.copula.be1)*(c.copula.be2.mix1-1)
                                                         +mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(c.copula2.be1be2.mix1))*dS1eta1*dS2eta2.2)*dereta1derb1, dereta2derb2.2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -( 
                
                VC$weights*VC$indLI*( -mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.thet.mix1*derteta.derteta.st-c.copula.thet*derteta.derteta.st)^2 
                                      + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(bit1.th2ATE.mix1-bit1.th2ATE)*derteta.derteta.st^2 
                                      + rotConst*mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.thet.mix1-c.copula.thet)*der2teta.derteta.stteta.st )
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -(
                
                crossprod(VC$weights*VC$indLI*c(rotConst*(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.mix1-c.copula.be1)*(c.copula.thet.mix1-c.copula.thet) 
                                                          + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t.mix1-c.copula2.be1t))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)
        )
        
        be2.rho <- -(
                
                crossprod(VC$weights*VC$indLI*c(rotConst*(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be2)*(c.copula.thet.mix1-c.copula.thet) 
                                                          + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be2t))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)+
                        
                crossprod(VC$weights*VC$indLI*c(rotConst*(-mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.mix1-1)*(c.copula.thet.mix1-c.copula.thet) 
                                                                  + mm(p2-p2.2+p00.mix1-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t.mix1))*dS2eta2.2*derteta.derteta.st)*dereta2derb2.2, X3)
                
        )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
}
#IL
if(sum(VC$indIL)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indIL*log(mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indIL*(c(mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(dS1eta1-c.copula.be1*dS1eta1)*dereta1derb1
                +c(c.copula.be1.mix2*dS1eta1.2-dS1eta1.2)*dereta1derb1.2)
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indIL*(c(mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula.be2.mix2*dS2eta2)*dereta2derb2
                -c(c.copula.be2*dS2eta2)*dereta2derb2)
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*( VC$indIL*(c(mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c.copula.theta.mix2
                -c.copula.theta)
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                
                crossprod(VC$weights*VC$indIL*c(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(dS1eta1-c.copula.be1*dS1eta1)^2 
                                                + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(d2S1eta1-c.copula2.be1*dS1eta1^2-c.copula.be1*d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIL*c( mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(dS1eta1-c.copula.be1*dS1eta1)  )*VC$X1)*der2.par1 ) ) )+
                        
                        
                crossprod(VC$weights*VC$indIL*c(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-dS1eta1.2+c.copula.be1.mix2*dS1eta1.2)^2 
                                                + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-d2S1eta1.2+c.copula2.be1.mix2*dS1eta1.2^2+c.copula.be1.mix2*d2S1eta1.2))*dereta1derb1.2, dereta1derb1.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIL*c( mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-dS1eta1.2+c.copula.be1.mix2*dS1eta1.2)  )*VC$X1.2)*der2.par1 ) ) )+
                        
                        
                        crossprod(VC$weights*VC$indIL*c(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-dS1eta1.2+c.copula.be1.mix2*dS1eta1.2)*(dS1eta1-c.copula.be1*dS1eta1)
                        )*dereta1derb1, dereta1derb1.2)+
                        
                        crossprod(VC$weights*VC$indIL*c(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-dS1eta1.2+c.copula.be1.mix2*dS1eta1.2)*(dS1eta1-c.copula.be1*dS1eta1)
                        )*dereta1derb1.2, dereta1derb1)
        )
        
        be2.be2 <-  -(
                
                crossprod(VC$weights*VC$indIL*c(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.mix2*dS2eta2-c.copula.be2*dS2eta2)^2
                                                + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2.mix2-c.copula2.be2)*dS2eta2^2 
                                                + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2.mix2-c.copula.be2)*d2S2eta2 )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIL*c( mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2.mix2-c.copula.be2)*dS2eta2 
                        )*VC$X2)*der2.par2 ) ) )
        )
        
        be1.be2 <- -(
                
                crossprod(VC$weights*VC$indIL*c((-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be1)*(c.copula.be2.mix2-c.copula.be2)+
                                                         mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula2.be1be2) )*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indIL*c((-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.mix2-1)*(c.copula.be2.mix2-c.copula.be2)+
                                                                 mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(c.copula2.be1be2.mix2))*dS1eta1.2*dS2eta2)*dereta1derb1.2, dereta2derb2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -( 
                
                VC$weights*VC$indIL*( -mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.thet.mix2*derteta.derteta.st-c.copula.thet*derteta.derteta.st)^2 
                                      + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(bit1.th2ATE.mix2-bit1.th2ATE)*derteta.derteta.st^2 
                                      + rotConst*mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.thet.mix2-c.copula.thet)*der2teta.derteta.stteta.st )
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        
        be1.rho <- -( 
                
                crossprod(VC$weights*VC$indIL*c(rotConst*(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(1-c.copula.be1)*(c.copula.thet.mix2-c.copula.thet) 
                                                          + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1t))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)+
                        
                        crossprod(VC$weights*VC$indIL*c(rotConst*(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.mix2-1)*(c.copula.thet.mix2-c.copula.thet)
                                                                  + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t.mix2) )*dS1eta1.2*derteta.derteta.st)*dereta1derb1.2, X3)
                )
        
        be2.rho <- -( 
                
                crossprod(VC$weights*VC$indIL*c(rotConst*(-mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.mix2-c.copula.be2)*(c.copula.thet.mix2-c.copula.thet)  
                                                          + mm(p1-p1.2+p00.mix2-p00, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t.mix2-c.copula2.be2t))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)
                )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
}
#II
if(sum(VC$indII)>1){
        #Likelihood
        l.par <- VC$weights*(VC$indII*log( mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr) ))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        #Gradient
        dl.dbe1 <- -VC$weights*(VC$indII*(c(mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                (c(c.copula.be1*dS1eta1)*dereta1derb1)
                -(c(c.copula.be1.mix1*dS1eta1) * dereta1derb1)
                -(c(c.copula.be1.mix2*dS1eta1.2)*dereta1derb1.2)
                +(c(c.copula.be1.2*dS1eta1.2)*dereta1derb1.2))
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indII*(c(mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula.be2*dS2eta2)*dereta2derb2
                -c(c.copula.be2.mix1*dS2eta2.2)*dereta2derb2.2
                -c(c.copula.be2.mix2*dS2eta2)*dereta2derb2
                +c(c.copula.be2.2*dS2eta2.2)*dereta2derb2.2)
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indII*(c(mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula.theta
                  -(c.copula.theta.mix1)
                  -(c.copula.theta.mix2)
                  +c.copula.theta.2))
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1*dS1eta1-c.copula.be1.mix1*dS1eta1)^2 
                                                + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1*dS1eta1^2+c.copula.be1*d2S1eta1-c.copula2.be1.mix1*dS1eta1^2-c.copula.be1.mix1*d2S1eta1))*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indII*c( mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be1*dS1eta1-c.copula.be1.mix1*dS1eta1)  )*VC$X1)*der2.par1 ) ) ) +
                        
                        crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula.be1.mix2*dS1eta1.2+c.copula.be1.2*dS1eta1.2)^2 + 
                                                                mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1.2*dS1eta1.2^2+c.copula.be1.2*d2S1eta1.2
                                                                                                                                            -c.copula2.be1.mix2*dS1eta1.2^2-c.copula.be1.mix2*d2S1eta1.2))*dereta1derb1.2, dereta1derb1.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indII*c( mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be1.2*dS1eta1.2-c.copula.be1.mix2*dS1eta1.2)  )*VC$X1.2)*der2.par1 ) ) ) +
                        
                        
                        crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-2)*(c.copula.be1*dS1eta1-c.copula.be1.mix1*dS1eta1)*(-c.copula.be1.mix2*dS1eta1.2+c.copula.be1.2*dS1eta1.2) 
                        )*dereta1derb1, dereta1derb1.2)+
                        crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-2)*(c.copula.be1*dS1eta1-c.copula.be1.mix1*dS1eta1)*(-c.copula.be1.mix2*dS1eta1.2+c.copula.be1.2*dS1eta1.2) 
                        )*dereta1derb1.2, dereta1derb1) 
                
        )
        
        be2.be2 <-  -(
                
               
                crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2*dS2eta2-c.copula.be2.mix2*dS2eta2)^2 
                                                + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2-c.copula2.be2.mix2)*dS2eta2^2
                                                +mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2-c.copula.be2.mix2)*d2S2eta2   )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indII*c( mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2-c.copula.be2.mix2)*dS2eta2 
                        )*VC$X2)*der2.par2 ) ) ) + 
                        
                        
                        crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.2*dS2eta2.2-c.copula.be2.mix1*dS2eta2.2)^2 
                                                        + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2.2-c.copula2.be2.mix1)*dS2eta2.2^2
                                                        +mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2.2-c.copula.be2.mix1)*d2S2eta2.2  )*dereta2derb2.2, dereta2derb2.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indII*c( mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.be2.2-c.copula.be2.mix1)*dS2eta2.2  )*VC$X2.2)*der2.par2 ) ) )+
                       
                        crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*((c.copula.be2*dS2eta2-c.copula.be2.mix2*dS2eta2)*(c.copula.be2.2*dS2eta2.2-c.copula.be2.mix1*dS2eta2.2)) )*dereta2derb2, dereta2derb2.2)+
                        crossprod(VC$weights*VC$indII*c(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*((c.copula.be2*dS2eta2-c.copula.be2.mix2*dS2eta2)*(c.copula.be2.2*dS2eta2.2-c.copula.be2.mix1*dS2eta2.2)) )*dereta2derb2.2, dereta2derb2)
        )
        be1.be2 <- -(
                #II
                #be1 be2
                crossprod(VC$weights*VC$indII*c((-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2-c.copula.be2.mix2)*(c.copula.be1-c.copula.be1.mix1)+ 
                                                         mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1be2)  )*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)+
                        #be1 b2.2 
                        crossprod(VC$weights*VC$indII*c((-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.2-c.copula.be2.mix1)*(c.copula.be1-c.copula.be1.mix1)+
                                                                 mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula2.be1be2.mix1) )*dS1eta1*dS2eta2.2)*dereta1derb1, dereta2derb2.2)+
                        #be1.2 be2
                        crossprod(VC$weights*VC$indII*c((-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2-c.copula.be2.mix2)*(c.copula.be1.2-c.copula.be1.mix2)+
                                                                 mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-c.copula2.be1be2.mix2) )*dS1eta1.2*dS2eta2)*dereta1derb1.2, dereta2derb2)+
                        #be1.2 #be2.2
                        crossprod(VC$weights*VC$indII*c((-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.2-c.copula.be2.mix1)*(c.copula.be1.2-c.copula.be1.mix2)+
                                                                 mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(c.copula2.be1be2.2) )*dS1eta1.2*dS2eta2.2)*dereta1derb1.2, dereta2derb2.2)
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -(
                #II
                VC$weights*VC$indII*( -mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.thet*derteta.derteta.st-c.copula.thet.mix1*derteta.derteta.st-c.copula.thet.mix2*derteta.derteta.st+c.copula.thet.2*derteta.derteta.st)^2 
                                      + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(bit1.th2ATE-bit1.th2ATE.mix1-bit1.th2ATE.mix2+bit1.th2ATE.2)*derteta.derteta.st^2 
                                      + rotConst*mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula.thet-c.copula.thet.mix1-c.copula.thet.mix2+c.copula.thet.2)*der2teta.derteta.stteta.st )
        )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        be1.rho <- -( 
                #II
                crossprod(VC$weights*VC$indII*c(rotConst*(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1-c.copula.be1.mix1)*(c.copula.thet-c.copula.thet.mix1-c.copula.thet.mix2+c.copula.thet.2) 
                                                          + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t-c.copula2.be1t.mix1))*dS1eta1*derteta.derteta.st)*dereta1derb1, X3)+
                        
                        crossprod(VC$weights*VC$indII*c(rotConst*(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be1.2-c.copula.be1.mix2)*(c.copula.thet-c.copula.thet.mix1-c.copula.thet.mix2+c.copula.thet.2) 
                                                                  + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t.2-c.copula2.be1t.mix2))*dS1eta1.2*derteta.derteta.st)*dereta1derb1.2, X3)
                )
        be2.rho <- -( 
                #II
                crossprod(VC$weights*VC$indII*c(rotConst*(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2-c.copula.be2.mix2)*(c.copula.thet-c.copula.thet.mix1-c.copula.thet.mix2+c.copula.thet.2) 
                                                          + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t-c.copula2.be2t.mix2))*dS2eta2*derteta.derteta.st)*dereta2derb2, X3)+
                        
                        crossprod(VC$weights*VC$indII*c(rotConst*(-mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula.be2.2-c.copula.be2.mix1)*(c.copula.thet-c.copula.thet.mix1-c.copula.thet.mix2+c.copula.thet.2) 
                                                                  + mm(p00-p00.mix1-p00.mix2+p00.2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t.2-c.copula2.be2t.mix1))*dS2eta2.2*derteta.derteta.st)*dereta2derb2.2, X3)
                )
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
        
        
       
}
#UI

if(sum(VC$indUI)>1){
        l.par <- VC$weights*( VC$indUI*( log( mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr) ) + log(-dS1eta1)+ log(Xd1P)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        #grad
        dl.dbe1 <- -VC$weights*(VC$indUI*(c(mm(c.copula.be1 - c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)) *
                                                  c( (c.copula2.be1 - c.copula2.be1.mix1) * dS1eta1) * dereta1derb1 +
                                                  c((dS1eta1)^(-1)*d2S1eta1) * dereta1derb1
                                          +c(Xd1P)^(-1)*der2eta1dery1b1
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indUI*( mm(c.copula.be1 - c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1) *
                                                   ( c(c.copula2.be1be2 * dS2eta2) * dereta2derb2 - c(c.copula2.be1be2.mix1 * dS2eta2.2) * dereta2derb2.2 )
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indUI*(c(mm(c.copula.be1 - c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula2.be1th - c.copula2.be1th.mix1))
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                #UI
                crossprod(VC$weights*VC$indUI*c(-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be1*dS1eta1-c.copula2.be1.mix1*dS1eta1)^2 
                                                + mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der3C.derp1p1p1-der3C.derp1p1p1.mix1)*dS1eta1^2 
                                                + mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1-c.copula2.be1.mix1)*d2S1eta1 
                                                -dS1eta1^-2*d2S1eta1^2 + dS1eta1^-1*d3S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUI*c(mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1-c.copula2.be1.mix1)*dS1eta1 + dS1eta1^-1*d2S1eta1)*VC$X1)*der2.par1 ) ) ) +
                        
                        crossprod(VC$weights*VC$indUI*c(-Xd1P^-2)*der2eta1dery1b1, der2eta1dery1b1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUI*c(Xd1P^-1)*VC$Xd1)*der2.par1 ) ) )
        )
        be2.be2 <-  -(
                #UI
                
                #no mix
                crossprod(VC$weights*VC$indUI*c(-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be1be2*dS2eta2)^2   +  
                                                        mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der2h.derp1p2)*dS2eta2^2 +
                                                        mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*d2S2eta2 )*dereta2derb2, dereta2derb2)+
                        
                        
                        
                        diag( colSums( t( t(VC$weights*VC$indUI*c( mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1be2)*dS2eta2 )*VC$X2)*der2.par2 ) ) )+
                        #mix
                        crossprod(VC$weights*VC$indUI*c(-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix1*dS2eta2.2)^2   +  
                                                                mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*-der2h.derp1p2.mix1*dS2eta2.2^2 +
                                                                mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*-c.copula2.be1be2.mix1*d2S2eta2.2 )*dereta2derb2.2, dereta2derb2.2) +  
                        
                        diag( colSums( t( t(VC$weights*VC$indUI*c(mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*-c.copula2.be1be2.mix1*dS2eta2.2)*VC$X2.2)*der2.par2 ) ) )+
                        #doppio
                        
                        
                        crossprod(VC$weights*VC$indUI*c(-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix1*c.copula2.be1be2)*dS2eta2*dS2eta2.2 )*dereta2derb2, dereta2derb2.2)+
                        crossprod(VC$weights*VC$indUI*c(-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix1*c.copula2.be1be2)*dS2eta2*dS2eta2.2 )*dereta2derb2.2, dereta2derb2)
                
        )
        be1.be2 <- -(
                #UI
                crossprod(VC$weights*VC$indUI*c((-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be1-c.copula2.be1.mix1)*(c.copula2.be1be2)
                                                 +mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(der2h.derp1p1) )*dS1eta1*dS2eta2)*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indUI*c((-mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be1-c.copula2.be1.mix1)*(-c.copula2.be1be2.mix1)
                                                         +mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)*(-der2h.derp1p1.mix1))*dS1eta1*dS2eta2.2)*dereta1derb1, dereta2derb2.2) 
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        
        d2l.rho.rho <- -( 
                #UI
                VC$weights*VC$indUI*( -mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be1t*derteta.derteta.st-c.copula2.be1t.mix1*derteta.derteta.st)^2 
                                      + mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der3C.derp1tetateta-der3C.derp1tetateta.mix1)*derteta.derteta.st^2 
                                      + rotConst*mm(c.copula.be1-c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be1t-c.copula2.be1t.mix1)*der2teta.derteta.stteta.st )
                )
        
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)
        
    be1.rho <- -(crossprod(VC$weights * VC$indUI * c( 
                                                       (rotConst *-mm(c.copula.be1 - c.copula.be1.mix1, min.pr = VC$min.pr, 
                                                            max.pr = VC$max.pr)^-2 * (c.copula2.be1 - c.copula2.be1.mix1) * 
                                                          (c.copula2.be1t - c.copula2.be1t.mix1) + mm(c.copula.be1 - 
                                                                                                        c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1 * 
                                                          (der3C.p1p1teta - der3C.p1p1teta.mix1)) * dS1eta1 * 
                                                       derteta.derteta.st) * dereta1derb1, X3))
    be2.rho <- -(crossprod(VC$weights * VC$indUI * c( 
                                                       (rotConst *-mm(c.copula.be1 - c.copula.be1.mix1, min.pr = VC$min.pr, 
                                                            max.pr = VC$max.pr)^-2 * (c.copula2.be1be2) * 
                                                          (c.copula2.be1t - c.copula2.be1t.mix1) + mm(c.copula.be1 - 
                                                                                                        c.copula.be1.mix1, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1 * 
                                                          (der2h.derp1teta)) * dS2eta2 * derteta.derteta.st) * 
                             dereta2derb2, X3) + crossprod(VC$weights * VC$indUI * 
                                                             c( (rotConst *-mm(c.copula.be1 - c.copula.be1.mix1, 
                                                                               min.pr = VC$min.pr, max.pr = VC$max.pr)^-2 * 
                                                                             (-c.copula2.be1be2.mix1) * (c.copula2.be1t - 
                                                                                                           c.copula2.be1t.mix1) + mm(c.copula.be1 - c.copula.be1.mix1, 
                                                                                                                                     min.pr = VC$min.pr, max.pr = VC$max.pr)^-1 * 
                                                                             (-der2h.derp1teta.mix1)) * dS2eta2.2 * derteta.derteta.st) * 
                                                             dereta2derb2.2, X3))    
        
        
        
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
       
        
}


#IU
if(sum(VC$indIU)>1){
        l.par <- VC$weights*(VC$indIU*( log( mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr) ) + log(-dS2eta2)+log(Xd2P)))
        res <- -sum(l.par)
        likelihood<-likelihood+ res
        #grad
        dl.dbe1 <- -VC$weights*(VC$indIU*( mm(c.copula.be2 - c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1) *
                                                   ( c(c.copula2.be1be2 * dS1eta1) * dereta1derb1 - c(c.copula2.be1be2.mix2 * dS1eta1.2) * dereta1derb1.2 )
        ))
        dl.dbe1 <- colSums(dl.dbe1)
        
        dl.dbe2 <- -VC$weights*(VC$indIU*(c(mm(c.copula.be2 - c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1)) *
                                                  ( c( (c.copula2.be2 - c.copula2.be2.mix2) * dS2eta2 ) * dereta2derb2 )
                                          + c((dS2eta2)^(-1)*d2S2eta2) * dereta2derb2
                                          + c(Xd2P)^(-1)*der2eta2dery2b2
        ))
        dl.dbe2 <- colSums(dl.dbe2)
        
        dl.dteta.st    <- -VC$weights*(VC$indIU*(c(mm(c.copula.be2 - c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^(-1))*(
                c(c.copula2.be2th - c.copula2.be2th.mix2))
        ))*X3
        dl.dteta.st <- colSums( dl.dteta.st)
        G <-G+ c( dl.dbe1, dl.dbe2, dl.dteta.st )
        
        #Hessian
        be1.be1 <- -(
                
                #IU
                crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be1be2*dS1eta1)^2
                                                + mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der2h.derp1p1*dS1eta1^2 )
                                                + mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*d2S1eta1)*dereta1derb1, dereta1derb1)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIU*c( mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*c.copula2.be1be2*dS1eta1  )*VC$X1)*der2.par1 ) ) )+
                        
                        
                        #mix
                        crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix2*dS1eta1.2)^2
                                                        + mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-der2h.derp1p1.mix2*dS1eta1.2^2 )
                                                        + mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-c.copula2.be1be2.mix2*d2S1eta1.2))*dereta1derb1.2, dereta1derb1.2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIU*c( (mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1)*(-c.copula2.be1be2.mix2*dS1eta1.2)  )*VC$X1.2)*der2.par1 ) ) )+
                        #doppio
                        
                        crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix2*c.copula2.be1be2)*dS1eta1*dS1eta1.2 )*dereta1derb1, dereta1derb1.2)+
                        crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix2*c.copula2.be1be2)*dS1eta1*dS1eta1.2 )*dereta1derb1.2, dereta1derb1)
        )
        be2.be2 <-  -(
                #IU
                
                
                crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be2-c.copula2.be2.mix2)^2*dS2eta2^2 
                                                -dS2eta2^-2*(d2S2eta2^2)+dS2eta2^-1*(d3S2eta2)+
                                                        + mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der2h.derp2p2-der2h.derp2p2.mix2)*dS2eta2^2 + 
                                                        mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2-c.copula2.be2.mix2)*d2S2eta2 )*dereta2derb2, dereta2derb2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIU*c( mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2-c.copula2.be2.mix2)*dS2eta2+
                                                                           dS2eta2^-1*d2S2eta2      )*VC$X2)*der2.par2 ) ) )+
                        
                        crossprod(VC$weights*VC$indIU*c(-Xd2P^-2)*der2eta2dery2b2, der2eta2dery2b2)  +  
                        
                        diag( colSums( t( t(VC$weights*VC$indIU*c(Xd2P^-1)*VC$Xd2)*der2.par2 ) ) )  
        )
        be1.be2 <- -(
                #IU
                
                crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*c.copula2.be1be2*dS1eta1*(c.copula2.be2*dS2eta2-c.copula2.be2.mix2*dS2eta2))*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indIU*c(mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der2h.derp1p2*dS1eta1*dS2eta2))*dereta1derb1, dereta2derb2)+
                        
                        crossprod(VC$weights*VC$indIU*c(-mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(-c.copula2.be1be2.mix2*dS1eta1.2)*(c.copula2.be2*dS2eta2-c.copula2.be2.mix2*dS2eta2))*dereta1derb1.2, dereta2derb2)+
                        crossprod(VC$weights*VC$indIU*c(mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(-der2h.derp1p2.mix2*dS1eta1.2*dS2eta2))*dereta1derb1.2, dereta2derb2)
                
                
        )
        if(VC$BivD %in% c("GAL180","C180","J180","G180","GAL90","C90","J90","G90","GAL270","C270","J270","G270") ) rotConst <- -1
        if(VC$BivD %in% VC$BivD2) rotConst <- VC$my.env$signind
        d2l.rho.rho <- -( 
                #IU
                VC$weights*VC$indIU*( -mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-2*(c.copula2.be2t*derteta.derteta.st-c.copula2.be2t.mix2*derteta.derteta.st)^2 
                                      + mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(der2h.derteta.teta.st -der2h.derteta.teta.st.mix2 )*derteta.derteta.st^2
                                      + rotConst*mm(c.copula.be2-c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1*(c.copula2.be2t-c.copula2.be2t.mix2)*der2teta.derteta.stteta.st ) 
                )
        rho.rho <- crossprod(X3*c(d2l.rho.rho), X3)

    be1.rho <- -(crossprod(VC$weights * VC$indIU * c( 
                                                       (rotConst *-mm(c.copula.be2 - c.copula.be2.mix2, min.pr = VC$min.pr, 
                                                            max.pr = VC$max.pr)^-2 * (c.copula2.be1be2) * 
                                                          (c.copula2.be2t - c.copula2.be2t.mix2) + mm(c.copula.be2 - 
                                                                                                        c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1 * 
                                                          (der2h.derp1teta)) * dS1eta1 * derteta.derteta.st) * 
                             dereta1derb1, X3) + crossprod(VC$weights * VC$indIU * 
                                                             c( (rotConst *-mm(c.copula.be2 - c.copula.be2.mix2, 
                                                                               min.pr = VC$min.pr, max.pr = VC$max.pr)^-2 * 
                                                                             (-c.copula2.be1be2.mix2) * (c.copula2.be2t - 
                                                                                                           c.copula2.be2t.mix2) + mm(c.copula.be2 - c.copula.be2.mix2, 
                                                                                                                                     min.pr = VC$min.pr, max.pr = VC$max.pr)^-1 * 
                                                                             (-der2h.derp1teta.mix2)) * dS1eta1.2 * derteta.derteta.st) * 
                                                             dereta1derb1.2, X3))
    be2.rho <- -(crossprod(VC$weights * VC$indIU * c( 
                                                       (rotConst *-mm(c.copula.be2 - c.copula.be2.mix2, min.pr = VC$min.pr, 
                                                            max.pr = VC$max.pr)^-2 * (c.copula2.be2 - c.copula2.be2.mix2) * 
                                                          (c.copula2.be2t - c.copula2.be2t.mix2) + mm(c.copula.be2 - 
                                                                                                        c.copula.be2.mix2, min.pr = VC$min.pr, max.pr = VC$max.pr)^-1 * 
                                                          (der2h.derp2teta - der2h.derp2teta.mix2)) * 
                                                       dS2eta2 * derteta.derteta.st) * dereta2derb2, X3))
    
        H <- H+ rbind( cbind( be1.be1        ,   be1.be2      ,      be1.rho    ), 
                       cbind( t(be1.be2)     ,   be2.be2      ,      be2.rho    ), 
                       cbind( t(be1.rho)     ,   t(be2.rho)   ,      rho.rho    ) )
       
}



########################################################################

if(VC$extra.regI == "pC") H <- regH(H, type = 1)

S.h  <- ps$S.h + monP2                                # hess
S.h1 <- 0.5*crossprod(params, ps$S.h)%*%params + monP # lik
S.h2 <- S.h%*%params + monP1                          # grad   

S.res <- likelihood # res
res   <- S.res + S.h1
G     <- G + S.h2
H     <- H + S.h  




if(VC$extra.regI == "sED") H <- regH(H, type = 2)   


list(value=res, gradient=G, hessian=H, S.h=S.h, S.h1=S.h1, S.h2=S.h2, 
     l=S.res, l.ln = l.ln, l.par=l.par,ps = ps, 
     eta1=eta1, eta2=eta2, etad=etad, etas1 = 1, etas2 = 1, 
     BivD=VC$BivD,               p1 = p1, p2 = p2, pdf1 = -dS1eta1, pdf2 = -dS2eta2,          
     c.copula.be2 = c.copula.be2,
     c.copula.be1 = c.copula.be1,
     c.copula2.be1be2 = c.copula2.be1be2, 
     dl.dbe1          = NULL,       
     dl.dbe2          = NULL,       
     dl.dteta.st      = NULL, 
     teta.ind2 = teta.ind2, teta.ind1 = teta.ind1,
     Cop1 = Cop1, Cop2 = Cop2, teta1 = teta1, teta2 = teta2,
     indNeq1 = indNeq1, indNeq2 = indNeq2,
     Veq1 = Veq1, Veq2 = Veq2, 
     k1 = VC$my.env$k1, k2 = VC$my.env$k2, monP2 = monP2) 

}



