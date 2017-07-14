adjCov <- function(x, id){

scores <- aCov(x)

scores <- aggregate.data.frame(scores,by=list(id),FUN=sum)[,-1]
nclusters <- dim(scores)[1]
meat   <- (nclusters-1)*var(scores)
covsan <- x$Vb %*% meat %*% x$Vb

x$Vb <- covsan

rm(scores, nclusters, meat, covsan)

x

}


