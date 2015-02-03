#' @export
calcK <- function(xodat, id, map.whole, ochr){
  nm.ochr <- sum(unlist(lapply(map.whole[ochr], length)))
  f.genoAH <- matrix(LETTERS[1:8], nrow = nm.ochr, 8, byrow = TRUE)
  genoAH <- convert2geno_wholechr(xodat[ochr], map.whole[ochr],
                                  id, f.genoAH)
  K <- DOQTL::kinship.alleles(t(genoAH))
  dimnames(K) <- list(paste0("id.", id), paste0("id.", id))
  return(K)
}
