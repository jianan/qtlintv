#' run simulations for DO type experiments
#'
#' For a fix set of parameters, generate DO pedigree and simulate
#' crossovers, randomly generate phenotypes with main QTL and minor
#' QTLs, then run DOQTL::scanone to do QTL mapping. save LOD score as
#' a matrix and marker informations to result.dir.
#'
#' @param method design method for the pedigree, choose from 'sub2' or
#' 'last2' or 'fixcc'.
#' @param para data.frame that saves all the parameter settings for the simulation
#' @param i.para which setting to use for the current run
#' @param i.simu which simulation is the currtent run at
#' @param n.simu_qtlpos number of phenotypes for each qtl position. suggested
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
#' @param realpb8 will the real 8 allele probabilities be used while
#' mapping.
#' @param noK do we control for Kinship matrix while mapping.
#' @param i.qtlpos which candidate qtlpos will be used.
#'
#' @export
run.do.single <- function(
    method=c("sub2", "last2", "fixcc"), para, i.para, i.simu, i.qtlpos,
    output.dir="DO.output/", result.dir="./", n.simu_qtlpos=1, n.simu=10000,
    qtl.chr=1, ochr=2:19, n.mqtl=10, n.ccgen=15,
    map.whole, qtl.allpos, f.geno.chr, snps, design=c("nosib", "random"),
    selc.method=c("byfamily","byindiv"),
    write.gp36=FALSE, cleanup=TRUE, realpb8=FALSE, noK=FALSE){

  method <- match.arg(method)
  design <- match.arg(design)

  result.dir <- add.slash(result.dir)
  if(!file.exists(result.dir))   dir.create(result.dir)

  output.dir <- add.slash(output.dir)
  if(file.exists(output.dir) & cleanup){
    warning(paste0("output.dir \"", output.dir,
                   "\" already exists, will not cleanup"),
            immediate.=TRUE)
    cleanup <- FALSE
  }
  if(!file.exists(output.dir))   dir.create(output.dir)

  n.gen <- para$n.gen[i.para]
  n.kids <- para$n.kids[i.para]
  n.sample <- para$n.sample[i.para]
  h.qtl <- para$h.qtl[i.para]
  h.kin <- para$h.kin[i.para]
  allele.freq <- para$allele.freq[i.para]
  seed <- para$seed[i.para]

  ## check that the QTLs themself are marker/pseudomarkers.
  stopifnot(all(qtl.allpos %in% snps$dist))
  n.qtlpos <- length(qtl.allpos)
  stopifnot(all(i.qtlpos <= n.qtlpos))
  stopifnot(i.simu <= n.simu)

  set.seed(seed)
  seeds <- matrix(runif(n.qtlpos*n.simu, 0, 1e+8), n.qtlpos, n.simu)
  set.seed(seeds[i.qtlpos, i.simu])

  allele <- 1:allele.freq

  tic <- proc.time() ## starting time

  cat("Simulating DO Pedigree... \n")
  ped <- simcross::sim_do_pedigree_fix_n(
      ngen=n.gen, nkids=n.kids, nccgen=n.ccgen,
      nsample=n.sample,
      method=method, design=design, selc.method=selc.method)
  id <- attr(ped, "last.gen.id")
  attr(ped, "last.gen.id") <- NULL

  while(1){
    cat("Simulating genotypes...  \n")
    xodat <- sim_from_pedigree_wholechr(pedigree=ped, map.whole)

    cat("Generating phenotypes...  \n")
    pheno <- gen_pheno_do_allqtl(xodat, map.whole, id,
                                 qtl.chr, qtl.allpos[i.qtlpos], allele,
                                 h.qtl, h.kin, n.mqtl, n.simu_qtlpos)

    min.sd <- min(apply(pheno, 2, sd))
    if(is.na(pheno[1,1]) | is.na(min.sd) | min.sd<0){
      cat("QTL has single genotype, re-do simulation of genotypes... \n")
    } else break
  }

  if(!noK){
    cat("Calculating Kinship matrix...  \n")
    K <- calcK(xodat, id, map.whole, ochr)
  }

  ########################################################################
  if(!realpb8){
    cat("Preparing for DOQTL::calc.genoprob...  \n")
    result <- prep_data_founders_chr(ped, xodat, map.whole, id, f.geno.chr, qtl.chr)

    cc.gen <- findccgen(ped)
    cat("Calculating genoprob...  \n")
    DOQTL::calc.genoprob(data=result$data, chr = qtl.chr, output.dir = output.dir,
                         plot = FALSE, array="other", sampletype="DO",
                         method="allele", snps=snps, founders=result$founders,
                         write.gp36=write.gp36, cc.gen=cc.gen)
    load(paste0(output.dir, "Chr", qtl.chr, ".founder.probs.Rdata"))
  }else{
    model.probs <- realprob8(map.whole, qtl.chr, xodat, id)
  }

  cat("QTL mapping...  \n")
  n.pheno <- ncol(pheno)
  if(noK){
    qtl <- DOQTL::scanone(pheno = pheno, pheno.col=1:n.pheno,
                          probs = model.probs, snps = snps)
  }else{
    qtl <- DOQTL::scanone(pheno = pheno, pheno.col=1:n.pheno,
                          probs = model.probs, K = K, snps = snps)
  }

  if(n.pheno == 1) qtl <- list(qtl)
  out <- scanone_do2qtl(qtl[[1]])
  pos <- out$pos
  snp <- rownames(out)

  LOD <- matrix(NA, nrow(out), n.pheno)
  for(i.qtl in 1:n.pheno){
    LOD[,i.qtl] <- scanone_do2qtl(qtl[[i.qtl]])[, "lod"]
  }
  colnames(LOD) <- colnames(pheno)

  toc <- proc.time() - tic ## run time

  str.qtlpos <- ifelse(length(i.qtlpos)==1, i.qtlpos, paste0(min(i.qtlpos), "-", max(i.qtlpos)))
  str.noK <- ifelse(noK, "noK","hasK")
  str.realP <- ifelse(realpb8, "realP", "calcP")
  file.result <- paste0(add.slash(result.dir), "result.", str.noK, ".", str.realP, ".",method,
                        ".para.", i.para, ".simu.", i.simu, ".qtl.", str.qtlpos, ".RData")
  save(LOD, toc, file=file.result)

  if(cleanup){
    files <- list.files(output.dir)
    if(length(files)>0)  file.remove(paste0(output.dir, files))
    file.remove(output.dir)
  }

  return(file.result)
}
