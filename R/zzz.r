dsbmessage <- function()
{
  mesg <-c(paste0(
" loaded package 'scglmmr' "
))
  return(mesg)
}

.onAttach <- function(lib, pkg)
{
  # startup message
  msg <- dsbmessage()
  if(!interactive())
    msg[1] <- paste("https://github.com/MattPM/scglmmr")
  packageStartupMessage(msg)
  invisible()
}
