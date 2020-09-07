mm <- function(ob, min.pr, max.pr){

  res <- ifelse(ob > max.pr, max.pr, ob) 
  res <- ifelse(res < min.pr, min.pr, res) 

  res
      
}


