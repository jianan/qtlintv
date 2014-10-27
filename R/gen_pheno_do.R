get_qtl_geno_do <- function(xodat, id, chr, pos, allele){
  if(missing(id)) id <- 1:length(xodat[[1]])
  g <- get_geno(xodat[[chr]][id], pos)
  g[g > 8] <- g[g > 8] - 8
  g <- (g[,1] %in% allele) + (g[,2] %in% allele)
  g
}

#' generate some minor qtls from other Chrs, not X and the QTL chr.
get_minor_qtl_do <- function(xodat, id, map, n.mqtl, qtl.chr){
  chrs <- names(map)
  chrs <- chrs[chrs!=qtl.chr & chrs!="X"]
  mqtl.chr <- sample(chrs, size=n.mqtl, replace=TRUE)
  mqtl.pos <- sapply(map[mqtl.chr], sample, size=1)
  mqtl.geno <- matrix(NA, length(id), n.mqtl)
  for(i in 1:n.mqtl){
    allele <- sample(1:8, size=sample(8,1))
    mqtl.geno[, i] <- get_qtl_geno_do(xodat, id, mqtl.chr[i], mqtl.pos[i], allele)
  }
  return(mqtl.geno)
}

#' generate pheno from pedigree and xo data.
#'
#' Add main qtl effect for a fixed qtl and multiple small effect QTLs
#' randomly selected from chr 2-19.
#' 
#' @export
gen_pheno_do <- function(xodat, map, id, qtl.chr, qtl.pos, allele,
                         h.qtl=0.3, h.kin=0.1, n.mqtl=10){
  
  ## get genotype at QTL
  qtlgeno <- get_qtl_geno_do(xodat, id, qtl.chr, qtl.pos, allele)
  mqtl.geno <- get_minor_qtl_do(xodat, id, map, n.mqtl, qtl.chr)

  mqtl.eff <- sample(c(-1, 1), n.mqtl, replace=TRUE)
  mqtl.eff <- mqtl.eff * runif(n=n.mqtl)

  ## generate phenotype
  n <- length(id)
  y.qtl <- qtlgeno
  s.qtl <- sqrt(h.qtl/(1-h.qtl-h.kin)/var(y.qtl)) 
  y.kin <- c(mqtl.geno %*% mqtl.eff)
  s.kin <- sqrt(h.kin/(1-h.qtl-h.kin)/var(y.kin))
  pheno <- s.qtl * y.qtl + s.kin * y.kin + rnorm(n)
  pheno <- matrix(pheno, nrow=n)
  rownames(pheno) <- paste0("id.", id)

  pheno
}

#' generate pheno for all the qtlpositons
#'
#' returns a matrix of size n x \code{n_qtl*n.simu}
#'
#' @export
gen_pheno_do_allqtl <- function(xodat, map, id, qtl.chr, qtl.allpos, allele, 
                                h.qtl=0.3, h.kin=0.1, n.mqtl=10, n.simu=1){

  ## get genotype at QTL
  n.qtl <- length(qtl.allpos)
  qtlgeno <- matrix(NA, length(id), n.qtl)
  for(i.qtl in 1:n.qtl){
    qtlgeno[, i.qtl] <- get_qtl_geno_do(xodat, id, qtl.chr,
                                        qtl.allpos[i.qtl], allele)
  }

  n <- length(id)
  pheno <- matrix(NA, n, n.qtl*n.simu)
  for(i.qtl in 1:n.qtl){
    y.qtl <- qtlgeno[, i.qtl]
    s.qtl <- sqrt(h.qtl/(1-h.qtl-h.kin)/var(y.qtl))
    for(i.simu in 1:n.simu){
      mqtl.geno <- get_minor_qtl_do(xodat, id, map, n.mqtl, qtl.chr)
      mqtl.eff <- sample(c(-1, 1), n.mqtl, replace=TRUE)
      mqtl.eff <- mqtl.eff * runif(n.mqtl)
      y.kin <- c(mqtl.geno %*% mqtl.eff)
      s.kin <- sqrt(h.kin/(1-h.qtl-h.kin)/var(y.kin))
      pheno[, (i.qtl-1)*n.simu+i.simu] <- s.qtl * y.qtl + s.kin * y.kin + rnorm(n)
    }
  }
  rownames(pheno) <- paste0("id.", id)
  if(n.simu==1)  colnames(pheno) <- paste0("pheno.", 1:n.qtl)
  else  colnames(pheno) <- paste0("pheno.", rep(1:n.qtl, each=n.simu), ".simu.", 1:n.simu)
  pheno
}
