#' create fake snps for DO founders, while keeping the same allele
#' pattern as the real snps.
#' @param n.snp Number of snps needed
#' @param f.geno founder genotype at all snps
#' @param p.missing proportion of missing data in the result
#' @return A n.snp x 8 matrix with values "A", "T", "G", "C" for
#' alleles and "N" for missing data.
#'
#' @export
simu_snp_do_founder <- function(n.snp, p.missing=0.03, f.geno){

  if(missing(p.missing)) p.missing <- mean(f.geno=="N")

  pattern <- matrix(NA, nrow(f.geno), 8)
  pattern <- apply(f.geno, 1, function(x){paste0(as.numeric(x==x[1]), collapse="")})
  x <- table(pattern)
  pattern.name <- names(x)
  x.pn <- x/sum(x)

  rand <- rmultinom(n=n.snp, size=1, prob=x.pn)
  res <- apply(rand==1, 2, which)
  res <- pattern.name[res]
  IND <- matrix(unlist(strsplit(res, ""))=="1", ncol=n.snp)

  alleles <- c("A", "T", "G", "C")
  mat <- matrix(0, 4, 4) ## transition prob matrix for snp
  rownames(mat) <- colnames(mat) <- alleles
  mat["A", "G"] <- 0.8
  mat["A", "C"] <- 0.2
  mat["T", "G"] <- 0.2
  mat["T", "C"] <- 0.8
  mat["G", "A"] <- 0.8
  mat["G", "T"] <- 0.2
  mat["C", "A"] <- 0.2
  mat["C", "T"] <- 0.8

  f.geno.new <- matrix(NA, n.snp, 8)
  f.geno.new[,1] <- sample(alleles, n.snp, replace=TRUE) ## founder 1
  for(i in 1:n.snp){
    g1 <- f.geno.new[i, 1]
    g2 <- sample(alleles, size=1, prob=mat[g1, ])
    ind.1 <- IND[,i]
    f.geno.new[i, ind.1] <- g1
    f.geno.new[i, !ind.1] <- g2
  }

  ## add missing values "N"
  f.geno.new[sample(1:(n.snp*8), size=n.snp*8*p.missing)] <- "N"
  if(n.snp == nrow(f.geno)) rownames(f.geno.new) <- rownames(f.geno)

  return(f.geno.new)
}
