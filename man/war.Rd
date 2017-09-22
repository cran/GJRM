\name{war}
\alias{war}
\docType{data}


\title{Civil war data}

\description{Civil war data from Fearon and Laitin (2003).}

\usage{
data(war)
}



\format{ \code{war} is a 6326 row data frame with the following columns 

\describe{
           
\item{onset}{equal to 1 for all country-years in which a civil war started.}
\item{instab}{equal to 1 if unstable government.}
\item{oil}{equal to 1 for oil exporter country.}
\item{warl}{equal to 1 if the country had a distinct civil war ongoing in the previous year.}
\item{gdpenl}{GDP per capita (measured as thousands of 1985 U.S. dollars) lagged one year.}
\item{ncontig}{equal to 1 for non-contiguous state.}
\item{nwstate}{equal to 1 for new state.}
\item{lpopl}{log(population size).}
\item{lmtnest}{log(mountainous).}
\item{ethfrac}{measure of ethnic fractionalization (calculated as the probability that two randomly drawn individuals from a country are not from the same ethnicity).}
\item{relfrac}{measure of religious fractionalization.}
\item{polity2l}{measure of political democracy (ranges from -10 to 10) lagged one year.}

}


} 




\source{

Data are from:

Fearon J.D., Laitin D.D. (2003), Ethnicity, Insurgency, and Civil War. \emph{The American Political Science Review}, 97, 75-90.

}



\examples{

\dontrun{ 

#########################################################
#########################################################

library("GJRM")

data("war", package = "GJRM")        

###################################################
# Bivariate brobit model with partial observability
###################################################  

reb.eq <- onset ~ instab + oil + warl + lpopl + lmtnest + ethfrac +
                  polity2l + s(gdpenl) + s(relfrac) 
gov.eq <- onset ~ instab + oil + warl + ncontig + nwstate + s(gdpenl)   

bpo <- gjrm(list(reb.eq, gov.eq), data = war, Model = "BPO",
            margins = c("probit", "probit"))
conv.check(bpo)

# perhaps model is to complex

set.seed(1)
sbpo <- summary(bpo)
sbpo$theta; sbpo$CItheta

# let's exclude the correlation parameter in fitting

bpo0 <- gjrm(list(reb.eq, gov.eq), data = war, Model = "BPO0",
             margins = c("probit", "probit"))
conv.check(bpo0)

summary(bpo0)


war.eq <- onset ~ instab + oil + warl +  ncontig + nwstate + lpopl +
                  lmtnest + ethfrac + polity2l + s(gdpenl) + s(relfrac)    
Probit <- gam(war.eq, family = binomial(link = "probit"), data = war)                 
summary(Probit)


coef(Probit)[(which(names(coef(Probit)) == "s(gdpenl).9"))]

coef(bpo0)[(which(names(coef(bpo)) == "s(gdpenl).9"))]


probitW <- bpoW <- bpoReb <- bpoGov <- NA
gdp.grid <- seq(0, 8)

median.values <- data.frame(t(apply(war, 2, FUN = median)))

for (i in 1:length(gdp.grid)){

newd <- median.values; newd$gdpenl <- gdp.grid[i]
eta1 <- predict(bpo0, eq = 1, newd)
eta2 <- predict(bpo0, eq = 2, newd)
probitW[i] <- predict(Probit, newd, type = "response") 
bpoW[i]    <- pnorm(eta1)*pnorm(eta2) 
bpoReb[i]  <- pnorm(eta1) 
bpoGov[i]  <- pnorm(eta2) 

}


plot(gdp.grid, probitW, type = "l", ylim = c(0, 0.55), lwd = 2, 
     col = "grey", xlab = "GDP per Capita (in thousands)", 
     ylab = "Pr(Outcome)", main = "Probabilities for All Outcomes", 
     cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.3)
lines(gdp.grid, bpoW,   lwd = 2)
lines(gdp.grid, bpoReb, lwd = 2, lty = 2)
lines(gdp.grid, bpoGov, lwd = 2, lty = 3)

#dev.copy(postscript, "probWAR.eps", width = 8)
#dev.off() 

}

#

}
