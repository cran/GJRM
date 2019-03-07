prediction.CLM <- function(x, object, newdata, type){

#object$knots <- NULL              #
#object$drop.unused.levels <- TRUE # TO BE MADE DYNAMYC

K1 <- object$VC$K1

if(!("(weights)" %in% names(newdata))){
	weights <- rep(1, dim(newdata)[1]) 
	newdata$weights <- weights
	names(newdata)[length(names(newdata))] <- "(weights)"
} else { 
	weights <- newdata[,"(weights)"]    
}


##### Setting-up quantities for prediction #####

if (!(c(object$formula[[1]][[2]]) %in% names(newdata))) {
	newdata$dependentVariable_FD <- rep(1)
		names(newdata)[which(names(newdata) == "dependentVariable_FD")] <- c(object$formula[[1]][[2]])
}

gam1.false <- eval(substitute(gam(object$formula[[1]], gamma = 1, weights = weights,
	data = newdata, knots = object$knots, fit = FALSE, drop.unused.levels = object$drop.unused.levels), 
	list(weights = weights))) 

X1 <- gam1.false$X ; X1 <- as.matrix(X1[, -1]) ; n <- dim(X1)[1]


##### Some checks #####

w.empty <- which(colnames(X1) == "")
	if (length(w.empty) != 0) n.check <- colnames(X1)[-w.empty] else n.check <- colnames(X1)

if (length(x$smooth) != 0) {
	n.coeff <- names(x$coefficients)[-(x$smooth[[1]]$first.para : length(names(x$coefficients)))]
} else {
	n.coeff <- names(x$coefficients)
}

w.coeff <- which(!(n.coeff %in% n.check))

if (length(w.coeff) != 0) stop("The covariates in newdata do not correspond to those of the fitted model.\nThis may happen, for example, when categorical covariates are included in the model specification but not all the levels are observed.")

###


#X1 <- model.matrix(x)
#	w.n1 <- which(names(x$model) %in% colnames(X1))
#	w.n2 <- which(colnames(X1) %in% names(x$model))
#	X1[, w.n2] <- as.matrix(x$model[, w.n1])

	eta1 <- X1 %*% x$coefficients
#	X1 <- as.matrix(X1) ; n <- dim(X1)[1]

cut  <- t(matrix(nrow = K1 - 1, ncol = n, object$coefficients[1 : (object$VC$K1 - 1)])) 
eta1 <- matrix(nrow = n, ncol = K1 - 1, eta1)
	lp1 <- cbind(cut - eta1)

pk <- probm(lp1, object$VC$margins[1], only.pr = FALSE, bc = TRUE)$pr

p1 <- matrix(nrow = n, ncol = object$VC$K1, 0)
	p1[, 1] <- pk[, 1]
	p1[, K1] <- 1 - pk[, K1 - 1]

if (K1 > 2) {
	for (k in 2 : (K1 - 1)) {
		p1[, k] <- pk[, k] - pk[, k - 1]

	}
}

pk <- cbind(pk, rep(1))


if (type == "response") {
	res <- list(p1.mass = p1,
                    p1.cum  = pk )
}

if (type == "lpmatrix") {
#	infty <- 1e+25
#
#	y1 <- gam1.false$y
#	sel <- model.matrix(~ as.factor(y1) - 1)
#	lp1 <- cbind(lp1, infty) ; lp1.sel <- rowSums(lp1 * sel)
#
#	Xs <- as.matrix(lp1.sel) %*% x$coefficients # IS THAT CORRECT?!?
#
#	res <- as.matrix(Xs)
#
#	res <- list(sel = sel,
#                   X1  = X1  )

	res <- list(X1 = X1)
}


return(res)

}

