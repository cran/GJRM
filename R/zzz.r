
.onAttach <- function(...) { 

  library(help=GJRM)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  hello <- paste("\nThis is GJRM ",version,".\nFor overview type 'help(\"GJRM-package\")'.\n",sep="")
  packageStartupMessage(hello)
  
}






