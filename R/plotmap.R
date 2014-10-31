#' plot map
#'
#' plot marker postion and pesudo-marker position of a map
#' 
#' @export
plotmap <- function(map, sep=10, main=""){
  
  ylim <- c(0, sep+1)
  xlim <- c(-1, 11)
  cuttoff <- seq(0, max(map), length=sep+1)
  
  plot(0,0,xlim=xlim, ylim=ylim, xaxt="n",type="n", yaxt="n", 
       xlab="Chromosome Position (cM)", ylab="index", main=main)
  at <- seq(xlim[1], xlim[2], length=13)
  axis(1, at=at, labels=round(at, 1), tick=FALSE)
  
  for(i in 1:10){
    snp <- names(map)[map >= cuttoff[i] & map <= cuttoff[i+1]]
    ind <- grep("loc-*[0-9]+(\\.[0-9]+)*$", snp)
    pseudomarker <- map[snp[ind]] - cuttoff[i]
    origmarker <- map[snp[-ind]] - cuttoff[i]
    abline(h=i)
    segments(x0=origmarker, x1=origmarker, y0=i-0.25, y1=i, col="black")
    segments(x0=pseudomarker, x1=pseudomarker, y0=i, y1=i+0.25, col="red")
  }

}
