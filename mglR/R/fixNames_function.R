#' Fix genes with missing names 
#'
#' \code{fixNames} fixes gene name element for those that were not found in ensembl
#'
#' This gives basic information on gene nomenclature for those genes that were missed by the \code{\link{buildFromNames}} function.  Use of ENSG identifier eliminates ambiguity while pairing with external gene name enhances human readability.  ENSG identifiers are required for many of the remaining functions in the package.
#'
#' @seealso \code{\link{missNames}}
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param ensgs A character vector of ensgs for missing genes. Note: These should be in the same order as the output from missNames!
#'
#' @examples
#' \dontrun{buildFromNames(c('CETP', 'APBO', '1ABCA')) -> myMgl}
#' \dontrun{missNames(myMgl)}
#' \dontrun{fixNames(myMgl, c('ENSG00000084674', 'ENSG00000165029'))}
#'
#' @section Details:
#' If a gene is not found, a warning message will appear when the function geneNames is run: 'Gene names missing: ...'. This must be corrected or no other elements can be filled in as the remaining elements all build off of the disambigous ENSG identifier.  There are two strategies to fix this.  The first is to check gene names for typos or use of less common colloquial names.  The second is to use the missNames and fixNames functions in this pacakge to fill in the missing ENSG identifiers.  Note: googling the colloquial gene name and 'gene cards' is an excellent way to find an ENSG id.  Genecards does an exceptional job of cataloging alternative colloquial names. 
#'
#'@export
 
fixNames <- function(mgl, ensgs){
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')	
	x = length(mgl)
	tmp <- list(); for (i in 1:x) {tmp[[i]] <- mgl[[i]][[1]]}
	lapply(tmp, class) -> check
	for (i in 1:length(which(check == 'character'))) {
		mgl[[which(check == 'character')[i]]][[1]] <- biomaRt::getBM(attributes = c('external_gene_id', 'ensembl_gene_id', 'description'), filters = 'ensembl_gene_id', values = ensgs[i], mart = mart)}
	return(mgl)
}

