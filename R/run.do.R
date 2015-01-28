#' run simulations for DO type experiments
#'
#' For a fix set of parameters, generate DO pedigree and simulate
#' crossovers, randomly generate phenotypes with main QTL and minor
#' QTLs, then run DOQTL::scanone to do QTL mapping. save LOD score as
#' a matrix and marker informations to result.dir.
#'
#' @param method design method for the pedigree, choose from 'sub2' or 'last2'.
#' @param para data.frame that saves all the parameter settings for the simulation
#' @param i.para which setting to use for the current run
#' @param j.simu which simulation is the currtent run at
#' @param n.simu number of phenotypes for each qtl position. suggested
#' to be 1, otherwise the result will be correlated.
#' @param qtl.chr chromosome for the qtl
#' @param ochr chromosomes other than qtl.chr, used to estimate kinship matrix
#' @param n.mqtl number of minor qtls for the polygenic effects
#' @param map.whole genetic map for all chromosomes
#' @param qtl.allpos named vector of all candidate qtl positions.
#' @param f.geno.chr founder genotype for QTL chromosome.
#' @param snps SNP information saved in a data.frame: name, chr and
#' position. Only for SNPs on the QTL chromosome.
#' @param design 'nosib' or 'random'
#' @param write.gp36 logic value
#' @param cleanup logic value
#'
#' @export
run.do <- function(method=c("sub2", "last2"), para, i.para, j.simu, output.dir="DO.output/", result.dir="./",
                   ## n.gen, n.kids, n.sample, h.qtl, h.kin, allele.freq, seed,
                   n.simu=1, qtl.chr=1, ochr=2:19, n.mqtl=10, n.ccgen=15, npairs_small=30, npairs_big=300,
                   map.whole, qtl.allpos, f.geno.chr, snps, design=c("nosib", "random"),
                   write.gp36=FALSE, cleanup=TRUE){

  method <- match.arg(method)
  design <- match.arg(design)

  ## check that the QTLs themself are marker/pseudomarkers.
  stopifnot(all(qtl.allpos %in% snps$dist))

  if(!file.exists(output.dir))   dir.create(output.dir)
  if(!file.exists(result.dir))   dir.create(result.dir)

  ## load("para.do.real.RData")
  n.gen <- para$n.gen[i.para]
  n.kids <- para$n.kids[i.para]
  n.sample <- para$n.sample[i.para]
  h.qtl <- para$h.qtl[i.para]
  h.kin <- para$h.kin[i.para]
  allele.freq <- para$allele.freq[i.para]
  seed <- para$seed[i.para]

  set.seed(seed)
  n.seeds <- ifelse(j.simu <= 100, 100, 1e+8)
  seeds <- runif(n.seeds, 0, 1e+8)
  set.seed(seeds[j.simu])
  allele <- 1:allele.freq

  tic <- proc.time() ## starting time

  cat("Simulating DO Pedigree... \n")
  ped <- simcross::sim_do_pedigree_fix_n(ngen=n.gen, nkids=n.kids, nccgen=n.ccgen,
                                         nsample=n.sample, npairs_small=npairs_small, npairs_big=npairs_big,
                                         method=method, design=design)
  id <- attr(ped, "last.gen.id")
  attr(ped, "last.gen.id") <- NULL

  while(1){
    cat("Simulating genotypes...  \n")
    xodat <- sim_from_pedigree_wholechr(pedigree=ped, map.whole)

    cat("Generating phenotypes...  \n")
    ## for each qtl position, generate multiple phenotypes, result is a matrix of n.sample x (n.simu*n.qtlpos)
    pheno <- gen_pheno_do_allqtl(xodat, map.whole, id,
                                 qtl.chr, qtl.allpos, allele,
                                 h.qtl, h.kin, n.mqtl, n.simu)

    min.sd <- min(apply(pheno, 2, sd))
    if(is.na(pheno[1,1]) | is.na(min.sd) | min.sd<0){
      cat("QTL has single genotype, re-do simulation of genotypes... \n")
    } else break
  }

  cat("Preparing for DOQTL::calc.genoprob...  \n")
  result <- prep_data_founders_chr(ped, xodat, map.whole, id, f.geno.chr, qtl.chr)

  cat("Calculating Kinship matrix...  \n")
  nm.ochr <- sum(unlist(lapply(map.whole[ochr], length)))
  f.genoAH <- matrix(LETTERS[1:8], nrow=nm.ochr, 8, byrow=TRUE)
  genoAH <- convert2geno_wholechr(xodat[ochr], map.whole[ochr], id, f.genoAH)
  K <- DOQTL::kinship.alleles(t(genoAH))
  dimnames(K) <- list(paste0("id.", id), paste0("id.", id))
  save(pheno, K, snps, file=paste0(output.dir, "pheno.K.snps.RData"))

  rm(xodat, f.geno.chr, f.genoAH, genoAH)

  ########################################################################
  cc.gen <- findccgen(ped)
  cat("Calculating genoprob...  \n")
  DOQTL::calc.genoprob(data=result$data, chr = qtl.chr, output.dir = output.dir,
                       plot = FALSE, array="other", sampletype="DO",
                       method="allele", snps=snps, founders=result$founders,
                       write.gp36=write.gp36, cc.gen=cc.gen)
  remove(result)

  cat("QTL mapping...  \n")
  load(paste0(output.dir, "Chr", qtl.chr, ".founder.probs.Rdata"))
  n.pheno <- ncol(pheno)
  qtl <- DOQTL::scanone(pheno = pheno, pheno.col=1:n.pheno,
                        probs = model.probs, K = K, snps = snps)

  out <- scanone_do2qtl(qtl[[1]])
  pos <- out$pos
  snp <- rownames(out)

  LOD <- matrix(NA, nrow(out), n.pheno)
  for(i.qtl in 1:n.pheno){
    LOD[,i.qtl] <- scanone_do2qtl(qtl[[i.qtl]])[, "lod"]
  }
  colnames(LOD) <- colnames(pheno)

  toc <- proc.time() - tic ## run time

  file.result <- paste0(result.dir, "result.para.", i.para, ".simu.", j.simu, ".RData")
  save(qtl.allpos, pos, snp, LOD, out, toc, file=file.result)

  if(cleanup){
    files <- list.files(output.dir)
    if(length(files)>0)  file.remove(paste0(output.dir, "/", files))
    file.remove(output.dir)
  }

  return(file.result)
}
