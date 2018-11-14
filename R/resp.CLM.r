resp.CLM <- function(res){

m <- min(res) ; M <- max(res)

seq.test <- seq(m, M)
seq.real <- sort(unique(res), decreasing = FALSE)

if (length(seq.test) != length(seq.real)) stop("The ordinal response has one (or more) missing level(s).")
if (!(all(seq.test == seq.real) == TRUE)) stop("The ordinal response has one (or more) missing level(s).")

K <- length(seq.test)


return(K)

}

