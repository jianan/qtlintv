get_qtl_geno <- function(xodat, id, chr, pos){
  if(missing(id)) id <- 1:length(xodat[[1]])
  g <- get_geno(xodat[[chr]][id], pos)
  g <- (g[,1] + g[,2] - 1)
  g
}

#' generate some minor qtls from other Chrs, not X and the QTL chr.
get_minor_qtl <- function(xodat, id, map, n.mqtl, qtl.chr){
  chrs <- names(map)
  chrs <- chrs[chrs!=qtl.chr & chrs!="X"]
  mqtl.chr <- sample(chrs, size=n.mqtl, replace=TRUE)
  mqtl.pos <- sapply(map[mqtl.chr], sample, size=1)
  mqtl.geno <- matrix(NA, length(id), n.mqtl)
  for(i in 1:n.mqtl){
    mqtl.geno[, i] <- get_qtl_geno(xodat, id, mqtl.chr[i], mqtl.pos[i])
  }
  return(mqtl.geno)
}

#' generate pheno for ail data
#'
#' @export
gen_pheno_ail<- function(xodat, map, id, qtl.chr, qtl.pos,
                         h.qtl=0.3, h.kin=0.1, n.mqtl=10){

  ## get genotype at QTL
  qtlgeno <- get_qtl_geno(xodat, id, qtl.chr, qtl.pos)
  mqtl.geno <- get_minor_qtl(xodat, id, map, n.mqtl, qtl.chr)

  mqtl.eff <- sample(c(-1, 1), n.mqtl, replace=TRUE)
  mqtl.eff <- mqtl.eff * runif(n=n.mqtl)

  ## generate phenotype
  n <- length(id)
  y.qtl <- qtlgeno
  s.qtl <- sqrt(h.qtl/(1-h.qtl-h.kin)/var(y.qtl))
  y.kin <- c(mqtl.geno %*% mqtl.eff)
  s.kin <- sqrt(h.kin/(1-h.qtl-h.kin)/var(y.kin))
  pheno <- s.qtl * y.qtl + s.kin * y.kin + rnorm(n)

  pheno
}
