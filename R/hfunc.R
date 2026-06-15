hfunc <- function(p1, p2, theta, copula = "N", dof = NULL){

	res.h <- copgHsAT(p1, p2, teta = theta, BivD = copula, Ln = TRUE, par2 = dof, 
	                  min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.9999999)

	list(h1 = res.h$c.copula.be1, h2 = res.h$c.copula.be2)

}


