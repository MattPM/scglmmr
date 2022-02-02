# source: https://github.com/MattPM/scglmmr
# author: Matt Mulè
# email: mattmule@gmail.com

scglmmr_message <- function()
{
  mesg <-c(paste0(
" loaded package 'scglmmr' "
))
  return(mesg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- scglmmr_message()
  if(!interactive())
    msg[1] <- paste("https://github.com/MattPM/scglmmr")
  packageStartupMessage(msg)
  invisible()
}
