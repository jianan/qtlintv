#' @export
map_df2list <- function(snps){
  stopifnot(!any(is.na(match(c("snp", "chr", "dist"), names(snps) ))))
  
  chrs <- unique(snps$chr)
  ## chrs <- c(1:19, "X")
  map <- NULL
  for(i in 1:length(chrs)){
    id <- which(snps$chr == chrs[i])
    map[[i]] <- structure(snps[id, "dist"], .Names=snps[id, "snp"])
    attr(map[[i]], "class") <- ifelse(chrs[i]!="X", "A", "X")
  }
  names(map) <- chrs
  class(map) <- "map"
  return(map)
}

#' @export
map_list2df <- function(map){
  res <- data.frame(snp=unlist(lapply(map, names)),
                    chr=rep(names(map), nmar(map)),
                    dist=unlist(map),
                    stringsAsFactors=FALSE)
  return(res)
}

#' @export
scanone_rel2qtl <- function(dat){
  out <- data.frame(chr=dat$chr, pos=dat$dist, lod=dat$p/(2*log(10)),
                    stringsAsFactors=FALSE)
  class(out) <- c("scanone", "data.frame")
  return(out)
}

#' @export
scanone_do2qtl <- function(qtl, include.x=FALSE){
  chr <- qtl$lod$A[,2]
  pos <- qtl$lod$A[,4]
  lod <- qtl$lod$A$neg.log10.p * 2
  snp.id <- qtl$lod$A[,1]
  if(include.x){
    chr <- c(chr, qtl$lod$X[,2])
    pos <- c(pos, qtl$lod$X[,4])
    lod <- c(lod, qtl$lod$X$neg.log10.p * 2)
    snp.id <- c(snp.id, qtl$lod$X[,1])
  }
  out.qtl <- data.frame(chr=chr, pos=pos, lod=lod,
                        stringsAsFactors=FALSE)
  rownames(out.qtl) <- snp.id
  class(out.qtl) <- c("scanone", "data.frame")
  out.qtl
}

##' generate a mice genetic map with one gape region on CHR  
##' @param n.mar.cM Number of markers per cM
gen_map <- function(n.mar.cM, gap.chr=1, gap.start=40, gap.end=60,
                    anchor.tel=TRUE, include.x=FALSE, sex.sp=FALSE, eq.spacing=FALSE){

  ## Chr length from snps in megamuga
  chr.len <- c(100, 104, 83, 89, 91, 80, 90, 77, 76, 78, 89, 64, 
               68, 67, 60, 58, 62, 60, 57, 80) 
  names(chr.len) <- c(1:19, "X")
  if(!include.x) chr.len <- chr.len[1:19]
  chr.n.mar <- chr.len * n.mar.cM + 1
  map <- sim.map(len=chr.len, n.mar=chr.n.mar, anchor.tel=anchor.tel,
                 include.x=include.x, sex.sp=sex.sp, eq.spacing=eq.spacing)
  map <- addblank(map, gap.chr, gap.start, gap.end)
  return(map)
}

##' add a gap region into a chr by removing markers.
addblank <- function(map, chr, gap.start, gap.end){
  x <- map[[chr]]
  class.save <- class(x)
  index <- which(x <= gap.start | x >= gap.end)
  x <- x[index]
  class(x) <- class.save
  map[[chr]] <- x
  return(map)
}

