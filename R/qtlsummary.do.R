#' summarize of do simulation results
#'
#' take LOD score and marker infomation, return coverage and interval
#' width for LOD support interval, Bayes interval, lod.drop and so on. 
#' 
#' @export
#' 
qtlsummary.do <- function(LOD, pos, snp, qtl.pos, lod.thr=3, drop=1.5, out.qtl){
  ## to-do: write functions to calc bayesint and lodint for 100 phenos at one time.
  
  require(qtl)

  stopifnot(nrow(LOD) == length(pos))
  stopifnot(length(snp) == length(pos))
  stopifnot(nrow(out.qtl) == length(pos))
  stopifnot(length(qtl.pos) == 1)
  names(pos) <- snp
  n.pheno <- ncol(LOD)
  
  max.lod <- apply(LOD, 2, max, na.rm = TRUE)
  LOD <- LOD[, max.lod >= lod.thr, drop=FALSE]
  max.lod <- max.lod[max.lod >= lod.thr]
  max.pos <- pos[apply(LOD, 2, which.max)]
  nr <- ncol(LOD)
  
  intv <- matrix(NA, 8, nr)
  for(i.simu in 1:ncol(LOD)){
    out.qtl$lod <- LOD[, i.simu]
    li <- qtl::lodint(out.qtl, drop=drop)
    bi <- qtl::bayesint(out.qtl)
    li.ep <- qtl::lodint(out.qtl, drop=drop, expandtomarkers=TRUE)
    bi.ep <- qtl::bayesint(out.qtl, expandtomarkers=TRUE)
    intv[, i.simu] <- c(bi[c(1,nrow(bi)), "pos"],
                        li[c(1,nrow(li)), "pos"],
                        bi.ep[c(1,nrow(bi.ep)), "pos"],
                        li.ep[c(1,nrow(li.ep)), "pos"])
  }
  intv <- t(intv)
  intv <- cbind(intv, max.lod, max.pos)
  colnames(intv) <- c("bi.left", "bi.right",
                      "li.left", "li.right",
                      "bi.ep.left", "bi.ep.right",
                      "li.ep.left", "li.ep.right",
                      "maxlod","maxpos")

  width <- intv[,c(2,4,6,8)] - intv[,c(1,3,5,7)]
  cover <- (intv[,c(1,3,5,7)] - qtl.pos) <= 0 & (intv[,c(2,4,6,8)] - qtl.pos) >= 0
  res <- c(apply(cover, 2, mean), apply(cover, 2, sd)/sqrt(nr),
           apply(width, 2, mean), apply(width, 2, sd)/sqrt(nr))
  names(res) <- paste(rep(c("bi","li","bi.ep","li.ep"), 4),
                      rep(c("cover","cover.sd", "width", "width.sd"),each=4),
                      sep=".")
  
  res.loddrop <- lod.drop(pos, snp, LOD, qtl.pos, probs=0.95, lod.thr=lod.thr)
  res <- c(res, res.loddrop)
  res["power"] <- nr/n.pheno
  res["maxlod"] <- median(max.lod)
  res["maxpos"] <- median(max.pos)
  attr(res, "intv") <- intv
  return(res)
}
