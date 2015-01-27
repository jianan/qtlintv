#' run simulations for AIL
#'
#' For a fix set of parameters, generate AIL pedigree and simulate
#' crossovers, randomly generate phenotypes with main QTL and minor
#' QTLs, then run QTLRel::scanOne to do QTL mapping. save LOD score as
#' a matrix and marker infomations to result.dir.
#'
#' @param n.simu number of phenotypes for the current run
#' @inheritParams run.do
#'
#' @export
run.ail <- function(method=c("sub2", "last2"), para, i.para, n.simu, n.mar.cM=2, result.dir="./",
                    qtl.chr=1, ochr=2:19, n.mqtl=10, design=c("nosib", "random")){

  method <- match.arg(method)
  design <- match.arg(design)

  if(!file.exists(result.dir))   dir.create(result.dir)

  stopifnot(colnames(para) == c("n.gen","n.kids","n.sample","h.qtl","h.kin","qtl.pos","seed"))
  ## load("para.do.real.RData")
  n.gen <- para$n.gen[i.para]
  n.kids <- para$n.kids[i.para]
  n.sample <- para$n.sample[i.para]
  h.qtl <- para$h.qtl[i.para]
  h.kin <- para$h.kin[i.para]
  gap.len <- para$gap.len[i.para]
  qtl.pos <- para$qtl.pos[i.para]
  seed <- para$seed[i.para]

  gap.start <- gap.middle - gap.len/2
  gap.end <- gap.middle + gap.len/2

  set.seed(seed)

  tic <- proc.time() ## starting time

  for(i.simu in 1:n.simu){

    cat(i.simu, "/", n.simu, "\n")
    ## Simulating AIL Pedigree...
    ped <- simcross::sim_ail_pedigree_fix_n(ngen=n.gen, nkids=n.kids,
                                            nsample=n.sample, npairs_small=30, npairs_big=300,
                                            method=method, design=design)
    id <- attr(ped, "last.gen.id")
    attr(ped, "last.gen.id") <- NULL

    ## Generate marker map...
    map <- gen_map(n.mar.cM=n.mar.cM, gap.chr=gap.chr, gap.start=gap.start, gap.end=gap.end,
                   eq.spacing=eq.spacing)
    map.df <- map_list2df(map)

    ## Simulating genotypes...
    xodat <- sim_from_pedigree_wholechr(pedigree=ped, map)

    ## convert xodata to marker genotypes, only for the last generation
    geno.chr <- convert2geno_wholechr(xodat[qtl.chr], map[qtl.chr], id)
    mapdf.chr <- subset(map.df, chr==qtl.chr)

    ## G matrix
    gdat <- convert2geno_wholechr(xodat[ochr], map[ochr], id)
    gm <- genMatrix(gdat)

    ## generate positions for mapping
    pseudo.pos <- seq(from=max(0, gap.start-gap.len/2),
                      to=min(max(map[[qtl.chr]]), gap.end+gap.len/2), by=step)
    if(!(qtl.pos %in% mapdf.chr$dist | qtl.pos %in% pseudo.pos)){
      ## makesure that the QTLs themself are marker/pseudomarkers.
      pseudo.pos <- c(pseduo.pos, qtl.allpos)
    }
    pos <- data.frame(snp=paste0("c", qtl.chr, ".loc", 1:length(pseudo.pos)),
                      chr=qtl.chr, dist=pseudo.pos)
    index <- pos$dist %in% mapdf.chr$dist   ## remove duplicate pesudo-marker.
    pos <- rbind(pos[!index, ], mapdf.chr)
    pos <- pos[order(pos$dist), ]

    ## init LOD matrix
    if(i.simu == 1) LOD <- matrix(NA, nrow(pos), n.simu)

    prDat <- genoProb(gdat=geno.chr, gmap=mapdf.chr, gr=n.gen, pos=pos)

    pheno <- gen_pheno_ail(xodat, map, id, qtl.chr, qtl.pos, h.qtl, h.kin, n.mqtl)
    vc <- QTLRel::estVC(y=pheno, v=list(AA=gm$AA, DD=NULL, HH=NULL, AD=NULL, MH=NULL,
                                     EE=diag(length(pheno))))
    out.rel <- scanOne(y=pheno, prd=prDat, vc=vc)
    out.qtl <- scanone_rel2qtl(out.rel)
    LOD[, i.simu] <- out.qtl$lod
  }

  attr(LOD, "pos") <- pos$dist
  attr(LOD, "snp") <- as.character(pos$snp)
  toc <- proc.time() - tic ## run time

  file.result <- paste0(result.dir, "result.para.", i.para, ".RData")
  save(LOD, out.qtl, toc, file=file.result)

  return(file.result)
}
