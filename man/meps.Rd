\name{meps}
\alias{meps}
\docType{data}


\title{MEPS data}

\description{2008 MEPS data.}

\usage{
data(meps)
}



\format{ \code{meps} is a 18592 row data frame with the following columns 

\describe{
        
\item{bmi}{body mass index.}
\item{age}{age in years.}
\item{gender}{equal to 1 if male.}
\item{race}{levels: 2 white, 3 black, 4 native American, 5 others.}
\item{education}{years of education.}
\item{health}{levels: 5 excellent, 6 very good, 7 good, 8 fair, 9 poor.}
\item{limitation}{equal to 1 if health limits physical activity.}
\item{region}{levels: 2 northeast, 3 mid-west, 4 south, 5 west.}
\item{private}{equal to 1 if individual has private health insurance.}
\item{visits.hosp}{equal to 1 if at least one visit to hospital outpatient departments.}
\item{diabetes}{equal to 1 if diabetic.}
\item{hypertension}{equal to 1 if hypertensive.}
\item{hyperlipidemia}{equal to 1 if hyperlipidemic.}
\item{income}{income (000's).}

}


} 




\source{

The data have been obtained from http://www.meps.ahrq.gov/. 


}




\examples{

\dontrun{ 

###################################################
###################################################

library("GJRM")
data("meps", package = "GJRM") 

###################################################
# Bivariate brobit models with endogenous treatment
###################################################

treat.eq <- private ~ s(bmi) + s(income) + s(age) + s(education) +
                      as.factor(health) + as.factor(race) +
                      as.factor(limitation) + as.factor(region) + 
                      gender  + hypertension + hyperlipidemia + diabetes
out.eq <- visits.hosp ~ private + s(bmi) + s(income) + s(age) + 
                        s(education) + as.factor(health) + 
                        as.factor(race) + as.factor(limitation) + 
                        as.factor(region) + gender + hypertension + 
                        hyperlipidemia + diabetes

f.list <- list(treat.eq, out.eq) 
mr     <- c("probit", "probit")
bpN    <- gjrm(f.list, data = meps, margins = mr, Model = "B")
bpF    <- gjrm(f.list, data = meps, margins = mr, BivD = "F", Model = "B")
bpC0   <- gjrm(f.list, data = meps, margins = mr, BivD = "C0", Model = "B")
bpC180 <- gjrm(f.list, data = meps, margins = mr, BivD = "C180", Model = "B")
bpJ0   <- gjrm(f.list, data = meps, margins = mr, BivD = "J0", Model = "B")
bpJ180 <- gjrm(f.list, data = meps, margins = mr, BivD = "J180", Model = "B")
bpG0   <- gjrm(f.list, data = meps, margins = mr, BivD = "G0", Model = "B")
bpG180 <- gjrm(f.list, data = meps, margins = mr, BivD = "G180", Model = "B")

conv.check(bpJ0)

AIC(bpN, bpF, bpC0, bpC180, bpJ0, bpJ180, bpG0, bpG180) 

set.seed(1)
summary(bpJ0, cex.axis = 1.6, 
        cex.lab = 1.6, cex.main = 1.7)

#dev.copy(postscript, "contplot.eps")
#dev.off()

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 2), 
    cex.axis = 1.6, cex.lab = 1.6)
plot(bpJ0, eq = 1, seWithMean = TRUE, scale = 0, shade = TRUE, 
     pages = 1, jit = TRUE)

#dev.copy(postscript, "spline1.eps")
#dev.off() 

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2, 2), 
    cex.axis = 1.6, cex.lab = 1.6)
plot(bpJ0, eq = 2, seWithMean = TRUE, scale = 0, shade = TRUE, 
     pages = 1, jit = TRUE)

#dev.copy(postscript, "spline2.eps")
#dev.off() 

set.seed(1)
AT(bpJ0, nm.end = "private", hd.plot = TRUE, cex.axis = 1.5, 
   cex.lab = 1.5, cex.main = 1.6)

#dev.copy(postscript, "hd.plotAT.eps")
#dev.off()

AT(bpJ0, nm.end = "private", type = "univariate")

AT(bpJ0, nm.end = "private", type = "naive")

}

#

}







