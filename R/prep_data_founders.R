#' prepare for DOQTL::calc.genoprob
#'
#' generate list 'data' and 'founders' from simulated DO pedigree and
#' crossovers.
#' 
#' @export
prep_data_founders <- function(ped, xodat, map, id, f.geno){
  
  geno <- convert2geno_wholechr(xodat, map, id, f.geno)
  sex <- ifelse(ped[id, "sex"]==1, "M", "F")
  gen <- paste0(ifelse(ped[id,"do"]==1, "DO", "CC"), ped[id, "gen"])
  gen[substr(gen,1,2)=="CC"] <- "CC"
  rownames(geno) <- names(sex) <- names(gen) <- paste0("id.", ped[id, "id"])
  
  founder.ped <- ped[ped[, "gen"] %in% c(0, 1) & ped[, "do"]==0, ]
  founder.sex <- ifelse(founder.ped[, "sex"]==1, "M", "F")
  founder.id <- founder.ped[, "id"]
  founder.geno <- convert2geno_wholechr(xodat, map, founder.id, f.geno)
  names(founder.sex) <- rownames(founder.geno) <- paste0("id.", founder.id)

  ## founder.code
  F1.ped <- ped[ped[, "gen"] == 1 & ped[, "do"]==0, ]
  mother.letter <- LETTERS[F1.ped[,"mom"]]
  father.letter <- LETTERS[F1.ped[,"dad"]-8]
  F1.code <- apply(cbind(mother.letter, father.letter), 1,
                   function(x){paste(sort(x), collapse="")})
  founder.code <- rep(NA, nrow(founder.ped))
  names(founder.code) <- founder.id
  founder.code[as.character(F1.ped[,"id"])] <- F1.code
  names(founder.code) <- paste0("id.", founder.id)
  F0.code <- rep(paste0(LETTERS[1:8], LETTERS[1:8]), 2)
  founder.code[1:16] <- F0.code

  data <- list(geno = geno, sex = sex, gen = gen)
  founders <- list(geno = founder.geno, sex = founder.sex, code = founder.code)
  return(list(data=data, founders=founders))
}

#' prepare for DOQTL::calc.genoprob
#'
#' generate list 'data' and 'founders' from simulated DO pedigree and
#' crossovers for a chromosome.
#' 
#' @export
prep_data_founders_chr <- function(ped, xodat, map, id, f.geno, chr){
  f.geno.chr <- f.geno[names(map[[chr]]),]
  result <- prep_data_founders(ped, xodat[chr], map[chr], id, f.geno.chr)
  return(result)
}
