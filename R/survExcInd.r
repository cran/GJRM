survExcInd <- function(n, cens, type.cens, hrate, d.lchrate, d.rchrate, d.lchrate.td, d.rchrate.td){ 
  
  ####################################
  ####################################
  # Survival indices for excess hazard 
  ####################################
  ####################################

  # Three things to look out for
  
  ### 1. none at all -- if there is at least one excess hazard constant it means we are in a survival setting, so we end up in this function. 
  ###    Check that the user provide all of the necessary quantities before checking whether their length is ok.
  
  ### 2. too little -- the user should input survival indices only for the correspondingly censored observations (if not what is the meaning of inserted quantity?) 
  
  ### 3. too many -- the user might provide excess hazard quantities which are not needed, this too is covered for
  
  # ---------------------------------- NON MIXED CENSORING ---------------------------------- #
  
  if(type.cens != 'mixed'){
  
  # Letf-truncated observations are not admitted in non-mixed censoring setting
  if( !is.null(d.lchrate.td) || !is.null(d.rchrate.td) ) stop( 'Left truncated observations are admitted only when type.cens = \'mixed\' so d.lchrate.td and d.rchrate.td are not needed here.')
    
  # For uncensored observations (in any type of non mixed censoring)  
  if( length(hrate) != sum(cens == 1))  stop(  ifelse(length(hrate) == 0, 'In an excess hazard setting with uncensored observations you must provide hrate.', paste('Length of excess hazard vector doesn\'t match number of uncensored observations in dataset.\n There are', sum(cens == 1), 'such observations in the dataset and', length(hrate), 'excess hazard values were provided.'  ))  )  
  
  
  # For left censored observations only
  if( type.cens == 'L' && length(d.lchrate) != sum(cens == 0) ) stop( ifelse(length(d.lchrate) == 0, 'In an excess hazard setting with only uncensored and left censored observations you must provide d.lchrate.', paste('Length of cumulative excess hazard vector evaluated at left censoring time doesn\'t match number of left censored observations in dataset.\n There are', sum(cens == 0), 'such observations in the dataset and', length(d.lchrate), 'cumulative excess hazard values were provided.'  )) )
  if( type.cens == 'L' && !is.null(d.rchrate) ) stop( 'In an excess hazard setting with only uncensored and right censored observations only hrate need be provided.' )
  
  # For right censored observations only
  if( type.cens == 'R' && (!is.null(d.lchrate) || !is.null(d.rchrate)) ) stop( 'In an excess hazard setting with only uncensored and right censored observations only hrate need be provided.' ) 

  
  # For interval censored observations only
  if( type.cens == 'I' && (length(d.lchrate) != sum(cens == 0) || length(d.rchrate) != sum(cens == 0) ) ){
    
    if( length(d.lchrate) == 0 && length(d.rchrate) == 0 ) stop('In an excess hazard setting with only uncensored and interval censored observations you must provide both d.lchrate and d.rchrate.')
    if( length(d.lchrate) == 0 ) stop('In an excess hazard setting with only uncensored and interval censored observations you must provide d.lchrate in addition to d.rchrate.')
    if( length(d.rchrate) == 0 ) stop('In an excess hazard setting with only uncensored and interval censored observations you must provide d.rchrate in addition to d.lchrate.')
    
    if( length(d.lchrate) != sum(cens == 0) ) stop( paste('Length of cumulative excess hazard vector evaluated at left censoring time doesn\'t match number of interval censored observations in dataset.\n There are', sum(cens == 0), 'such observations in the dataset and', length(d.lchrate), 'cumulative excess hazard values were provided.'  ) )
    if( length(d.rchrate) != sum(cens == 0) ) stop( paste('Length of cumulative excess hazard vector evaluated at right censoring time doesn\'t match number of interval censored observations in dataset.\n There are', sum(cens == 0), 'such observations in the dataset and', length(d.rchrate), 'cumulative excess hazard values were provided.'  ) )
    
  }
  
  
  hrate.t <- rep(0, n) 
  if(!is.null(hrate)) hrate.t[cens == 1] <- hrate 
  hrate <- hrate.t
  rm(hrate.t)
  
  d.lchrate.t <- rep(0, n) 
  if( !is.null(d.lchrate) && type.cens %in% c('L', 'I') ) d.lchrate.t[cens == 0] <- d.lchrate 
  d.lchrate <- d.lchrate.t
  rm(d.lchrate.t)
  
  d.rchrate.t <- rep(0, n) 
  if(!is.null(d.rchrate) && type.cens == 'I') d.rchrate.t[cens == 0] <- d.rchrate 
  d.rchrate <- d.rchrate.t
  rm(d.rchrate.t)
  
  d.lchrate.td <- rep(0, n) 
  d.rchrate.td <- rep(0, n) 
  
  
}
    
  # --------------------------------- MIXED CENSORING --------------------------------- #
  
  if(type.cens == 'mixed'){
    
    if( length(hrate) != sum(cens %in% c("U", "UT"))) stop( ifelse(length(hrate) == 0, 'In an excess hazard setting with uncensored observations you must provide hrate.', paste('Length of excess hazard vector doesn\'t match number of uncensored and/or left truncated uncensored observations in dataset.\n There are', sum(cens %in% c("U", "UT")), 'such observations in the dataset and', length(hrate), 'excess hazard values were provided.'  )) )
    if( length(d.lchrate) != sum(cens %in% c("L", "I")) ) stop( ifelse(length(d.lchrate) == 0, 'In an excess hazard setting with left censored and/or interval censored observations you must provide d.lchrate.', paste('Length of cumulative excess hazard vector evaluated at left censoring time doesn\'t match number of left censored and/or interval censored observations in dataset.\n There are', sum(cens %in% c("L", "I")), 'such observations in the dataset and', length(d.lchrate), 'cumulative excess hazard values were provided.'  )) )
    if( length(d.lchrate.td) != sum(cens %in% c("LT", "IT")) ) stop(ifelse(length(d.lchrate.td) == 0, 'In an excess hazard setting with left truncated left censored and/or left truncated interval censored observations you must provide d.lchrate.td.', paste('Length of cumulative excess hazard vector evaluated at left truncated left censoring time \n doesn\'t match number of left truncated left censored and/or left truncated interval censored observations in dataset.\n There are', sum(cens %in% c("LT", "IT")), 'such observations in the dataset and', length(d.lchrate.td), 'cumulative excess hazard \n values were provided.'  )) )
    if( length(d.rchrate) != sum(cens == "I") ) stop( ifelse(length(d.rchrate) == 0, 'In an excess hazard setting with interval censored observations you must provide d.rchrate.', paste('Length of cumulative excess hazard vector evaluated at right censoring time \n doesn\'t match number of interval censored observations in dataset.\n There are', sum(cens == "I"), 'such observations in the dataset and', length(d.rchrate), 'cumulative excess hazard \n values were provided.'  )) )
    if( length(d.rchrate.td) != sum(cens == "IT") ) stop( ifelse(length(d.rchrate.td) == 0, 'In an excess hazard setting with left truncated interval censored observations you must provide d.rchrate.td.', paste('Length of cumulative excess hazard vector evaluated at left truncated right censoring time \n doesn\'t match number of left truncated interval censored observations in dataset.\n There are', sum(cens == "IT"), 'such observations in the dataset and', length(d.rchrate.td), 'cumulative excess hazard \n values were provided.'  )) )
    
    
    hrate.t <- rep(0, n) 
    if(!is.null(hrate)) hrate.t[cens %in% c("U", "UT")] <- hrate
    hrate <- hrate.t
    rm(hrate.t)
    
    d.lchrate.t <- rep(0, n)
    if(!is.null(d.lchrate)) d.lchrate.t[cens %in% c("L", "I")] <- d.lchrate 
    d.lchrate <- d.lchrate.t
    rm(d.lchrate.t)
    
    d.rchrate.t <- rep(0, n)
    if(!is.null(d.rchrate)) d.rchrate.t[cens == "I"] <- d.rchrate 
    d.rchrate <- d.rchrate.t
    rm(d.rchrate.t)
    
    d.lchrate.td.t <- rep(0, n)
    if(!is.null(d.lchrate.td)) d.lchrate.td.t[cens %in% c("LT", "IT")] <- d.lchrate.td 
    d.lchrate.td <- d.lchrate.td.t
    rm(d.lchrate.td.t)
    
    d.rchrate.td.t <- rep(0, n)
    if(!is.null(d.lchrate.td)) d.rchrate.td.t[cens == "IT"] <- d.rchrate.td 
    d.rchrate.td <- d.rchrate.td.t
    rm(d.rchrate.td.t)
    
    
  }
  
  
  # ----------------------------------------------------------------------------------------------------------------------------------- #
  
  L <- list(hrate = hrate, d.lchrate = d.lchrate, d.rchrate = d.rchrate, d.lchrate.td = d.lchrate.td, d.rchrate.td = d.rchrate.td)
  
  L
  
}

