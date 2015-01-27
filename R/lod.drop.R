#' LOD-drop for lod support interval to reach {95\%} coverage rate.
#'
#' @export
lod.drop <- function(pos, snp, LOD, qtl.pos, probs=0.95, lod.thr=3, eps=1e-8){

  n.phenos <- ncol(LOD)
  colnames(LOD) <- 1:n.phenos
  rownames(LOD) <- snp

  ## LOD score at QTL position
  wh <- which.min(abs(pos-qtl.pos))
  qtl.lod <- LOD[wh, ]
  max.lod <- apply(LOD, 2, max, na.rm = TRUE)
  max.pos <- pos[apply(LOD, 2, which.max)]

  names(pos) <- snp
  pos.marker <- pos[-grep("loc-*[0-9]+(\\.[0-9]+)*$", snp)]
  candi.pos <- pos.marker[pos.marker - qtl.pos < -eps]
  mn.left <- names(candi.pos)[which.max(candi.pos)]
  candi.pos <- pos.marker[pos.marker - qtl.pos > +eps]
  mn.right <- names(candi.pos)[which.min(candi.pos)]

  qtl.lod.marker <- numeric(n.phenos)
  for(j in 1:n.phenos){
    if(max.pos[j] <= pos[mn.left]){
      qtl.lod.marker[j] <- LOD[mn.left, j]
    }else if(max.pos[j] >= pos[mn.right]){
      qtl.lod.marker[j] <- LOD[mn.right, j]
    }else{
      qtl.lod.marker[j] <- max.lod[j]
    }
  }

  lod.drop1 <- max.lod - qtl.lod
  lod.drop2 <- max.lod - qtl.lod.marker
  lod.drop <- apply(cbind(lod.drop1, lod.drop2), 1, min)

  ## only use those with bigger lod score.
  lod.drop1 <- lod.drop1[max.lod >= lod.thr]
  lod.drop <- lod.drop[max.lod >= lod.thr]

  ## use bootstrap to estimate sd(standard errors) of quantile of x.
  bs <- function(x, qt=0.95, n.bs=1000){
    stopifnot(length(qt)==1)
    n.x <- length(x)
    res <- numeric(n.bs)
    for(i in 1:n.bs){
      res[i] <- quantile(sample(x, n.x, replace=TRUE), qt)
    }
    return(sd(res))
  }

  drop <- quantile(lod.drop1, probs)
  drop.sd <- bs(lod.drop1, qt=probs, n.bs=1000)
  drop.ep <- quantile(lod.drop, probs)
  drop.ep.sd <- bs(lod.drop, qt=probs, n.bs=1000)
  if(length(probs)==1){
    res <- c(drop, drop.sd, drop.ep, drop.ep.sd)
    names(res) <- c("drop", "drop.sd", "drop.ep", "drop.ep.sd")
  }else{
    res <- cbind(drop, drop.sd, drop.ep, drop.ep.sd)
    colnames(res) <- c("drop", "drop.sd", "drop.ep", "drop.ep.sd")
    rownames(res) <- probs
  }
  return(res)
}
