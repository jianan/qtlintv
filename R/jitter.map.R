#' @export
jitter.map <- function(map, eps=1e-6){
  if(!any(class(map) == "map"))
      stop("Input must be a map")
  for(i in 1:length(map)){
    x <- map[[i]]
    while(1){
      ind <- which(diff(x) < eps*0.9)
      if(length(ind)>0) x[ind+1] <- x[ind+1] + eps else break
    }
    map[[i]] <- x
  }
  return(map)
}


## ## test
## eps <- 1e-7
## map.whole <- gen_map(n.mar.cM=1)
## chrs <- c(1:19,"X")
## map.whole <- map.whole[chrs]
## map.jj.my <- jitter.map(map.whole, eps)
## n.same <- numeric(length(chrs))
## names(n.same) <- chrs
## for(i.chr in chrs){
##   n.same[i.chr] <- (max(map.jj.my[[i.chr]] - map.whole[[i.chr]]))/eps + 1
## }
## n.same2 <- c(unlist(lapply(map.whole, function(x) max(table(x)))))
## all(abs(n.same - n.same2) < 0.5)

