#' Returns gene names with a given phenotype
#'
#' \code{makePhenotypeSearch} returns a character vector of gene names. 
#'
#' Of interest maybe groups of genes that have the same phenotype 
#'
#' @family output
#'
#' @seealso \code{\link{makePhenotypes}}
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param term Character vector representing phenotype of interest 
#'
#' @param phen List produced by \code{\link{makePhenotypes}}
#'
#' @param saveFile A logical flag indicating whether a csv file ('PhenotypeSearch_[term].csv') should be saved in the current directory
#'
#' @examples
#' makePhenotypes(myMgl, saveFile = FALSE) -> myPhenotypes
#' makePhenotypeSearch(myMgl, term = 'Mean corpuscular hemoglobin', 
#'     phen = myPhenotypes, saveFile = FALSE) -> myPhenotypeSearch
#'
#'@export
#'@importFrom utils write.csv
 
makePhenotypeSearch <- function(mgl, term = c(''), phen, saveFile = FALSE){

# use to display genes with certain GOID
as.character(phen[[3]][which(phen[[3]][,1] == term),2]) -> pnames

# Saving if flag
if (saveFile == TRUE)
write.csv(pnames, paste('PhenotypeSearch_', stringr::str_replace_all(term, " ", "_"), '.csv', sep = ""))

#message('\nPlease cite the NHGRI-EBI GWAS Catalog or GRASP as appropriate. \n \n NHGRI-EBI GWAS Catalog: \n Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H.The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006 \n \n GRASP: \n (1) Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390 genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP Build 2.0.0.0\n (2) Carey V (2016). grasp2db: grasp2db, sqlite wrap of GRASP 2.0. R package version 0.1.14.\n')

return(pnames)

}

