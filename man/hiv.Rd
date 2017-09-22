\name{hiv}
\alias{hiv}
\alias{hiv.polys}
\docType{data}


\title{HIV Zambian data}

\description{HIV Zambian data by region, together with polygons describing the regions' shapes.}

\usage{
data(hiv)
data(hiv.polys)
}



\format{ \code{hiv} is a 6416 row data frame with the following columns 



\describe{


\item{hivconsent}{binary variable indicating consent to test for HIV.}
\item{hiv}{binary variable indicating whether an individual is HIV positive.}
\item{age}{age in years.}
\item{education}{years of education.}
\item{region}{code identifying region, and matching \code{names(hiv.polys)}. It can take nine possible values: 1 central, 2 copperbelt, 3 eastern, 
             4 luapula, 5 lusaka, 6 northwestern, 7 northern, 8 southern, 9 western.}
\item{marital}{never married, currently married, formerly married.}
\item{std}{had a sexually transmitted disease.}
\item{highhiv}{had high risk sex.}
\item{condom}{used condom during last intercourse.}
\item{aidscare}{equal to 1 if would care for an HIV-infected relative.}
\item{knowsdiedofaids}{equal to 1 if know someone who died of HIV.}
\item{evertestedHIV}{equal to 1 if previously tested for HIV.}           
\item{smoke}{smoker.}         
\item{ethnicity}{bemba, lunda (luapula), lala, ushi, lamba, tonga, luvale,
lunda (northwestern), mbunda, kaonde, lozi, chewa, nsenga,
ngoni, mambwe, namwanga, tumbuka, other.}              
\item{language}{English, Bemba, Lozi, Nyanja, Tonga, other.}                     
\item{interviewerID}{interviewer identifier.}                        
\item{sw}{survey weights.}  

}



\code{hiv.polys} contains the polygons defining the areas in the format described below.
} 


\details{The data frame \code{hiv} relates to the regions whose boundaries are coded in \code{hiv.polys}.
\code{hiv.polys[[i]]} is a 2 column matrix, containing the vertices of the polygons defining the boundary of the ith 
region. \code{names(hiv.polys)} matches \code{hiv$region} (order unimportant).

}

\source{

The data have been produced as described in: 

McGovern M.E., Barnighausen T., Marra G. and Radice R. (2015), On the Assumption of Joint Normality in Selection Models: A Copula Approach Applied to Estimating HIV Prevalence. \emph{Epidemiology}, 26(2), 229-237.


}


\references{

Marra G., Radice R., Barnighausen T., Wood S.N. and McGovern M.E. (2017), A Simultaneous Equation Approach to Estimating HIV Prevalence with Non-Ignorable Missing Responses. \emph{Journal of the American Statistical Association}, 112(518), 484-496.


}



\examples{

\dontrun{ 

#########################################################
#########################################################

library("GJRM")

data("hiv", package = "GJRM")        
data("hiv.polys", package = "GJRM")  

#########################################################
#########################################################
## The stuff below is useful if the user wishes to employ  
## a Markov Random Field (MRF) smoother. It provides
## the instructions to set up polygons automatically
## and the dataset variable needed to fit a model with 
## MRF.
#########################################################
#########################################################
#
# ## hiv.polys was already created and
# ## made available via the call 
# ## data("hiv.polys", package = "GJRM") 
# ## hiv.polys was created using the code below
#
# obj <- readRDS("ZMB_adm1.rds") 
# ## RDS Zambian Level 1 file obtained from 
# ## http://www.gadm.org. 
#
# pol <- polys.setup(obj)
#
# hiv.polys <- pol$polys   
# name <- cbind(names(hiv.polys), pol$names1)
# name
#
## last step was to create a factor variable with range
## range(name[,1]) where the numerical values were linked   
## to the regions in name[, 2]. This is what was done in 
## the hiv dataset; see hiv$region. Specifically,
## the procedure used was
##
# reg <- NULL
# 
# for(i in 1:dim(hiv)[1]){
# 
# if(hiv$region[i] == "Central")       reg[i] <- 1
# if(hiv$region[i] == "Copperbelt")    reg[i] <- 2
# if(hiv$region[i] == "Eastern")       reg[i] <- 3
# if(hiv$region[i] == "Luapula")       reg[i] <- 4
# if(hiv$region[i] == "Lusaka")        reg[i] <- 5
# if(hiv$region[i] == "North-Western") reg[i] <- 6
# if(hiv$region[i] == "Northern")      reg[i] <- 7
# if(hiv$region[i] == "Southern")      reg[i] <- 8
# if(hiv$region[i] == "Western")       reg[i] <- 9
# 
# }
# 
# hiv$region <- as.factor(reg)
# 
# 
#########################################################
#########################################################

xt <- list(polys = hiv.polys) 

# neighbourhood structure info for MRF  
# to use in model specification

#########################################################
# Bivariate brobit model with non-random sample selection
#########################################################          
     
sel.eq <- hivconsent ~ s(age) + s(education) + s(wealth) + 
                       s(region, bs = "mrf", xt = xt, k = 7) + 
                       marital + std + age1sex_cat + highhiv + 
                       partner + condom + aidscare + 
                       knowsdiedofaids + evertestedHIV + 
                       smoke + religion + ethnicity + 
                       language + s(interviewerID, bs = "re")
 
out.eq <- hiv        ~ s(age) + s(education) + s(wealth) + 
                       s(region, bs = "mrf", xt = xt, k = 7) + 
                       marital + std + age1sex_cat + highhiv + 
                       partner + condom + aidscare + 
                       knowsdiedofaids + evertestedHIV + 
                       smoke + religion + ethnicity + 
                       language      

theta.eq <-          ~ s(region, bs = "mrf", xt = xt, k = 7)                       
    
fl <- list(sel.eq, out.eq, theta.eq)    
     
# the above model specification is fairly
# complex and it serves to illustrate the 
# flexibility of the modelling approach
     
bss <- gjrm(fl, data = hiv, BivD = "J90", Model = "BSS",
            margins = c("probit", "probit"))

conv.check(bss)

set.seed(1)
sb <- summary(bss)
sb

plot(bss, eq = 1, seWithMean = TRUE, scheme = 1,   
     scale = 0, pages = 1, jit = TRUE)
                    
plot(bss, eq = 2, seWithMean = TRUE, scheme = 1,
     scale = 0, pages = 1, jit = TRUE)


prev(bss, sw = hiv$sw, type = "naive") 

set.seed(1)
prev(bss, sw = hiv$sw, type = "univariate") 

prev(bss, sw = hiv$sw) 


lr <- length(hiv.polys) 
prevBYreg  <- matrix(NA, lr, 2)
thetaBYreg <- NA

for(i in 1:lr) {
prevBYreg[i,1] <- prev(bss, sw = hiv$sw, ind = hiv$region==i, 
                       type = "univariate")$res[2]
prevBYreg[i,2] <- prev(bss, sw = hiv$sw, ind = hiv$region==i)$res[2]
thetaBYreg[i]  <- bss$theta[hiv$region==i][1]
}


zlim <- range(prevBYreg)  # to establish a common prevalence range

par(mfrow = c(1, 3), cex.axis = 1.3)

polys.map(hiv.polys, prevBYreg[,1], zlim = zlim, lab = "",  
          cex.lab = 1.5, cex.main = 1.5, 
          main = "HIV - Imputation Model")
          
polys.map(hiv.polys, prevBYreg[,2], zlim = zlim, cex.main = 1.5, 
          main = "HIV - Selection Model")
          
polys.map(hiv.polys, thetaBYreg, rev.col = FALSE, cex.main = 1.7, 
          main = expression(paste("Copula parameter (",hat(theta),")")))
          
sb$CItheta[1,]
}

#

}





