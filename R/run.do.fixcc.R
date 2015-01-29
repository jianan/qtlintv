#' run simulations for DO type experiments, with the same pre-CC
#' structure as real experiment.
#'
#' For a fix set of parameters, generate DO pedigree and simulate
#' crossovers, randomly generate phenotypes with main QTL and minor
#' QTLs, then run DOQTL::scanone to do QTL mapping. save LOD score as
#' a matrix and marker informations to result.dir.
#'
#' @inheritParams run.do
#'
#' @export
run.do.fixcc <- function(para, i.para, j.simu, output.dir="DO.output/", result.dir="./",
                   n.simu=1, qtl.chr=1, ochr=2:19, n.mqtl=10,
                   map.whole, qtl.allpos, f.geno.chr, snps, design=c("nosib", "random"),
                   write.gp36=FALSE, cleanup=TRUE){

  design <- match.arg(design)

  ## check that the QTLs themself are marker/pseudomarkers.
  stopifnot(all(qtl.allpos %in% snps$dist))

  output.dir <- add.slash(output.dir)
  result.dir <- add.slash(result.dir)

  if(file.exists(output.dir) & cleanup){
    warning(paste0("output.dir \"", output.dir,
                   "\" already exists, will not cleanup"))
    cleanup <- FALSE
  }

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

  ## use fixed ccgen value
  alpha <- c(21, 64, 24, 10, 5, 9, 5, 3, 3)
  names(alpha) <- 4:12
  npairs <- sum(alpha)
  ccgen <- rep(as.numeric(names(alpha)), alpha)
  ped <- sim_do_pedigree(ngen = n.gen, npairs = npairs, ccgen=ccgen,
                         nkids_per = n.kids, design = design)
  id <- ped[ped[, "gen"] == n.gen & ped[, "do"] == 1, "id"]
  id <- sample(id, n.sample)

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
    if(length(files)>0)  file.remove(paste0(output.dir, files))
    file.remove(output.dir)
  }

  return(file.result)
}
