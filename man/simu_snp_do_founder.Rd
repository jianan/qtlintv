% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/simu_snp_do_founder.R
\name{simu_snp_do_founder}
\alias{simu_snp_do_founder}
\title{create fake snps for DO founders, while keeping the same allele
pattern as the real snps.}
\usage{
simu_snp_do_founder(n.snp, p.missing = 0.03, f.geno)
}
\arguments{
\item{n.snp}{Number of snps needed}

\item{p.missing}{proportion of missing data in the result}

\item{f.geno}{founder genotype at all snps}
}
\value{
A n.snp x 8 matrix with values "A", "T", "G", "C" for
alleles and "N" for missing data.
}
\description{
create fake snps for DO founders, while keeping the same allele
pattern as the real snps.
}

