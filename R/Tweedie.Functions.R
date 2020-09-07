########################################################################################
############ Functions for Tweedie CDF and derivatives
########################################################################################

###
# Gradient and hessian wrt mu, phi and p
derivPGamma_lag <- function(y, la, a, g, k, log = TRUE, deriv = 0){
  
  lga <- pgamma(y, -k * a, scale = g, log.p = TRUE)
  lpo <- dpois(k, la, log = TRUE)
  
  lP <- lga + lpo
  if( !log ){ lP <- exp(lP) }
  
  grad <- NULL
  Hess <- NULL
  if( deriv ){
    
    funcD <- function(.x) pgamma(y, -k * .x, scale = g)
    # nde <- numgh(funcD, a)
    
    nde <- list()
    nde$fi <- jacobian(function(x) funcD(x), x = a)
    nde$se <- jacobian(function(x) jacobian(function(x) funcD(x), x = x), x = a)
    
    # Deriv wrt lambda, alpha and gamma
    dpo_la <- dpois(k, la) * (k / la - 1)
    dga_ga <- ifelse(y > 0, - y / g^2 * exp( log(y / g) * (- k * a - 1) - y / g - lgamma(-k*a)), 0)
    grad <- c(dpo_la * exp( lga ),
              nde$fi * exp( lpo ), 
              dga_ga * exp( lpo )) 
    
    if( any(!is.finite(grad)) ) stop("nooo")
    
    if( deriv > 1 ){
      
      Hess <- matrix(c(exp( lga ) * (dpo_la * (k / la - 1) - k * exp( lpo ) / la^2), 
                       dpo_la * nde$fi, 
                       dpo_la * dga_ga,
                       dpo_la * nde$fi, 
                       nde$se * exp( lpo ), 
                       ifelse(y*k > 0, k * dga_ga * (digamma(-k*a) - log(y/g)) * exp( lpo ), 0), 
                       dpo_la * dga_ga,
                       ifelse(y*k > 0, k * dga_ga * (digamma(-k*a) - log(y/g)) * exp( lpo ), 0), 
                       ifelse(y*k > 0, y * exp(-y/g - (k*a + 1) * log(y/g) - lgamma(-k*a)) * (g * (1 - k*a) - y) / g^4 * exp( lpo ) , 0)
      ), 3, 3, byrow = TRUE)
    }
    
  }
  
  out <- list("p" = lP, "grad" = grad, "Hess" = Hess)
  
  return( out )
  
}

###
# Gradient and hessian wrt mu, phi and p
derivPGamma_mphip  <- function(y, mu, phi, p, k, log = TRUE, deriv = 0){
  
  la <- mu^(2-p) / ( phi * (2-p) )
  ga <- phi * (p-1) * mu^( p-1 )
  al <- (2-p) / ( 1-p )
  
  out <- derivPGamma_lag(y = y, la = la, a = al, g = ga, k = k, log = log, deriv = deriv)
  
  J11 <- la * (2 - p) / mu
  J12 <- J22 <- 0
  J13 <- (p - 1) * ga / mu
  J21 <- - la / phi
  J23 <- ga / phi
  J31 <- mu^(2-p) * (1 + (p-2) * log(mu)) / ((p-2)^2 * phi)
  J32 <- 1 / ((p-1)^2)
  J33 <- ga * (log(mu) + 1/(p-1))
  
  J <- matrix(c(J11, J12, J13, J21, J22, J23, J31, J32, J33), 3, 3, byrow = TRUE)
  
  if( deriv ){
    
    grad0 <- out$grad
    out$grad <- drop( J %*% grad0 )
    
    if( deriv > 1 ){
      
      Jmu <- matrix(c(J11 * (1 - p) / mu,
                      0,
                      J13 * (p - 2) / mu,
                      - J11 / phi,
                      0,
                      J13 / phi,
                      J31 * ( (2 - p)/mu ) + mu ^ (2 - p) / (p - 2) / phi / mu,
                      0,
                      ga / mu + J13 * (log(mu) + 1 / (p - 1))),
                    3, 3, byrow = TRUE)
      
      Jphi <- matrix(c(J21 * (2 - p) / mu,
                       0,
                       J23 * (p - 1) / mu,
                       - J21 * 2 / phi,
                       0,
                       0,
                       - J31 / phi,
                       0,
                       J23 * (log(mu) + 1 / (p-1))),
                     3, 3, byrow = TRUE)
      
      Jp <- matrix(c((J31 * (2 - p) - la) / mu,
                     0,
                     (J33 * (p - 1) + ga)/ mu,
                     - J31 / phi,
                     0,
                     J33 / phi,
                     - J31 * (log(mu) + 2/(p-2)) + mu^(2-p) * log(mu) / (phi*(p-2)^2),
                     - 2 / ( p - 1 )^3,
                     - ga / ((p - 1)^2) + J33 * (log(mu) + 1 / (p - 1))),
                   3, 3, byrow = TRUE)
      
      out$Hess <- J %*% out$Hess %*% t(J) + cbind(Jmu %*% grad0, Jphi %*% grad0, Jp %*% grad0)
      
    }
    
  }
  
  return(out)
  
}

##############################
# CDF of tweedie distribution
##############################
pTweed <- function(y, mu, phi, p, log = FALSE, eps = 39, deriv = 0){
  
  la <- mu^(2-p) / ( phi * (2-p) )
  ga <- phi * (p-1) * mu^( p-1 )
  al <- (2-p) / ( 1-p )
  
  if( y == 0 ){
    lpPT <- dpois(0, la, log = TRUE)
    lpPT <- if( !log ){ exp(lpPT) } 
    d1 <- d2 <- NULL
    if( deriv >= 1 ){ d1 <- - exp(-la) * c( (2-p)/mu*la, -la/phi, mu^(2-p)*(1+(p-2)*log(mu))/((p-2)^2*phi) ) }
    if( deriv >= 2 ){ d2 <- d1 %*% t(d1) / exp(-la) + 
                            matrix(c(d1[1]*(1-p)/mu, -d1[1]/phi, d1[3]*(2-p)/mu - mu^(2-p)/(p-2)/phi/mu*exp(-la),
                                     d1[2]*(2-p)/mu, -2*d1[2]/phi, -d1[3]/phi,
                                     d1[3]*(2-p)/mu+la/mu*exp(-la), -d1[3]/phi, -d1[3]*(log(mu)+2/(p-2)) - mu^(2-p)*log(mu)/((p-2)^2*phi)*exp(-la)),
                                                3, 3) }
    return( list("d0" = lpPT, "d1" = d1, "d2" = d2) )
  }
  
  # Mode of Poisson is lambda and get log density at mode
  N <- N0 <- max(floor(la), 1)
  mxlpP <- lP <- dpois(N, la, log = TRUE)
  
  # Get log-probability contribution at Poisson mode, and use it to avoid underflow
  co0 <- mxlpP + pgamma(y, shape = - N * al, scale = ga, log.p = TRUE)
  lpPT <- 0
  x <- NA # New log-probability value will be stored here
  
  # Get derivative distribution at Poisson mode
  d1 <- d2 <- NULL
  if( deriv ){
    
    der <- derivPGamma_mphip(y = y, mu = mu, phi = phi, p = p, k = N, deriv = deriv) # Storage for gradient
    d1 <- der$grad
    d2 <- der$Hess # Storage for Hessian
    
  }
  
  # Sum from mode until N_max (latter is unknown)
  while ( lP > mxlpP - eps  ){
    
    N <- N + 1
    
    lP <- dpois(N, la, log = TRUE)
    x <- lP + pgamma(y, shape = - N * al, scale = ga, log.p = TRUE)
    
    # # We have a new max log-prob, so we reset the constant used to avoid underflow
    if ( x > co0 ){
      lpPT <- lpPT + co0 - x
      co0 <- x
    }
    
    lpPT <- log( exp(lpPT) + exp(x - co0) )
    
    if( deriv ){
      
      der <- derivPGamma_mphip(y = y, mu = mu, phi = phi, p = p, k = N, deriv = deriv)
      d1 <- d1 + der$grad
      d2 <- d2 + der$Hess
      
    }
    
  }
  Nmax <- N
  
  N <- N0
  lP <- mxlpP 
  
  # Sum from mode until N_min (latter is unknown)
  while ( lP > mxlpP - eps && N > 0  ){
    
    N <- N - 1
    
    lP <- dpois(N, la, log = TRUE)
    x <- lP + pgamma(y, shape = - N * al, scale = ga, log.p = TRUE)
    
    # We have a new max log-prob, so we reset the constant used to avoid underflow
    if ( x > co0 ){
      lpPT <- lpPT + co0 - x
      co0 <- x
    }
    
    lpPT <- log( exp(lpPT) + exp(x - co0) )
    
    if( deriv ){
      
      der <- derivPGamma_mphip(y = y, mu = mu, phi = phi, p = p, k = N, deriv = deriv)
      d1 <- d1 + der$grad
      d2 <- d2 + der$Hess
      
    }
    
  }
  Nmin <- N
  
  lpPT <- lpPT + co0
  
  if( !log ) { lpPT <- exp( lpPT ) }
  
  return( list("d0" = lpPT, "d1" = d1, "d2" = d2) )
  
}
