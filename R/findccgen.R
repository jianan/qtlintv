#' @export
findccgen <- function(ped){
  gen <- ped[ped[,"do"]==0,"gen"]
  x <- table(gen)/2
  names(x) <- as.numeric(names(x)) - 3
  x <- x[as.character((max(gen)-3):1)]
  x <- c(x[1], diff(x))
  x <- x[x!=0]
  x <- x[as.character(sort(as.numeric(names(x))))]
  return(x)
}
