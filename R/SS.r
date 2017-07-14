SS <- function(offtemp, Sx){

          totnM <- length(unique(offtemp))
          off2 <- as.data.frame(cbind(1:length(offtemp),offtemp))
          off2 <- split(off2, off2$offtemp)
          St <- list()
          for(i in 1:totnM) St[[i]] <- Reduce("+", Sx[ off2[[i]][,-2] ])
          St
          
}
