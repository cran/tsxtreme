## Copyright (C) 2017 Thomas Lugrin
## tic and toc function, as in Matlab
##################################################

#####################################
## CLOCK FUNCTIONS
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
    tic <- proc.time()[type]         
    assign(".tic", tic, envir=baseenv())
    invisible(tic)
}

toc <- function(silent=FALSE)
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  if(!silent) print(toc - tic)
  invisible(toc -tic)
}