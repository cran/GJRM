inform.setup <- function(gam1, gam2, inform.cov, start.v1, start.v2, lsgam1, lsgam2){ 
  
scv <- pcv <- stermsgam1 <- stermsgam2 <- pv1 <- pv2 <- pos1 <- pos2 <- inds2 <- poss2 <- NULL
j   <- 1
Lposs2    <- list()
Lposs2PAR <- list()
poss2PAR  <- NA

ptermsgam1 <- attr(gam1$pterms, "term.labels")
ptermsgam2 <- attr(gam2$pterms, "term.labels")


pcv <- intersect(inform.cov, ptermsgam1) # check how many and if any info covs in the set of parametric ones
                                         # I do this only for ptermsgam1 because they are the same in both equations

pcv <- ptermsgam1[which(ptermsgam1 %in% pcv)] # order here depends on how list of
                                              # inform covs is provided and this is not good
                                              # hence the solution above; the order must be that of ptermsgam1


lp1 <- lp2 <- length(pcv) # total no. of covariates to consider

if(lp1 == 0) pcv <- NULL 




if(lp1 != 0){ 

    for(i in 1:lp1){ Lposs2PAR[[i]] <- grep(pcv[i], names(start.v1))  
                     poss2PAR[i] <- Lposs2PAR[[i]][1]
                   } # I do not add another one as equations are the same

                   start.v1[unlist(Lposs2PAR)] <- start.v2[unlist(Lposs2PAR)] <- (start.v1[unlist(Lposs2PAR)] + start.v2[unlist(Lposs2PAR)])/2
                 }







for(i in 1:lsgam1) stermsgam1 <- c(stermsgam1, gam1$smooth[[i]]$vn) # same as above
for(i in 1:lsgam2) stermsgam2 <- c(stermsgam2, gam2$smooth[[i]]$vn)

scv <- intersect(inform.cov, stermsgam1)

if(length(scv) == 0) scv <- NULL

# scv is not affected by above problem as scv is never actually used in the calculations





if(length(scv) != 0){

for(i in 1:lsgam1){ if(gam1$smooth[[i]]$vn %in% scv) pos1 <- c(pos1, gam1$smooth[[i]]$first.para:gam1$smooth[[i]]$last.para) }


for(i in 1:lsgam2){ if(gam2$smooth[[i]]$vn %in% scv) {pos2 <- c(pos2, gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para) 
                                                      
                                                       poss2[j]   <- gam2$smooth[[i]]$first.para
                                                      Lposs2[[j]] <- gam2$smooth[[i]]$first.para:gam2$smooth[[i]]$last.para

                                                      inds2[j] <- i
                                                      j <- j + 1
                                                      
 
                                                      }}
                                                      
                                                      
                                                      

start.v1[pos1] <- start.v2[pos2] <- (start.v1[pos1] + start.v2[pos2])/2 

                 }


start.v2 <- start.v2[-c(unlist(Lposs2PAR), pos2)] 


list(start.v1 = start.v1, start.v2 = start.v2, 
     lp1 = lp1, pcv = pcv, 
     par.pos1 = unlist(Lposs2PAR), par.pos2 = unlist(Lposs2PAR), 
     smo.pos1 = pos1, smo.pos2 = pos2, 
     poss2 = poss2, Lposs2 = Lposs2, Lposs2PAR = Lposs2PAR, poss2PAR = poss2PAR, 
     inds2 = inds2, pcv = pcv, scv = scv)

}

