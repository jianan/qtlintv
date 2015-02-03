#' @export
realprob8 <- function(map.whole, qtl.chr, xodat, id){
  nm <- length(map.whole[[qtl.chr]])
  f.genoAH <- matrix(LETTERS[1:8], nrow = nm, 8, byrow = TRUE)
  genoAH <- convert2geno_wholechr(xodat[qtl.chr], map.whole[qtl.chr], id, f.genoAH)
  n <- dim(genoAH)[1]
  m <- dim(genoAH)[2]
  g1 <- substr(genoAH, 1, 1)
  g2 <- substr(genoAH, 2, 2)
  prob <- array(NA, dim=c(n,m,8))
  for(i in 1:8){
    prob[, , i] <- (g1 == LETTERS[i]) + (g2 == LETTERS[i])
  }
  prob <- prob/2
  prob <- aperm(prob, c(1,3,2))
  if(is.numeric(id)) id <- paste0("id.", id)
  dimnames(prob) <- list(id, LETTERS[1:8], colnames(genoAH))
  return(prob)
}
