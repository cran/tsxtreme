## Copyright (C) 2017 Thomas Lugrin
## tic and toc functions are provided by tictoc package,
## but need a wrapper around the latter
##################################################

#####################################
## CLOCK FUNCTIONS
toc <- function(silent = FALSE)
{
  tt <- tictoc::toc(quiet = TRUE)#list of 3
  ttic <- tt$tic
  ttoc <- tt$toc
  if(!silent) print(ttoc - ttic)
  invisible(ttoc - ttic)
}