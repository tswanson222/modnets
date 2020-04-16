netTable <- function(fits, nodes = FALSE, d = 4, alpha = .05, 
                     names = NULL, order = TRUE, s = "res"){
  args <- as.list(match.call()[-1])
  FUN <- ifelse("ggm" %in% names(attributes(fits[[1]])), "modTable", "SURtable")
  if(FUN == "modTable"){args["s"] <- NULL}
  do.call(match.fun(FUN), args)
}
