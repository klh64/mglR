#' Idenify genes with missing names 
#'
#' \code{missNames} returns gene names that were not found in ensembl
#'
#' This gives a vector of gene names where basic information on gene nomenclature in particular identification of the unambiguous ENSG identifier was NOT found.
#'
#' @family elements
#' @seealso \code{\link{fixNames}}
#'
#' @param mgl Object of class 'mgl'; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' \dontrun{buildFromNames(c('CETP', 'APBO'), elements = 
#'	c('name', 'enst', 'location', 'antisense')) -> myMgl}
#' \dontrun{missNames(myMgl)}
#'
#' @section Details:
#' If a gene is not found, a warning message will appear: 'Gene names missing: ...'.  This must be corrected or no other elements can be filled in as the remaining elements all build off of the disambigous ENSG identifier.  There are two strategies to fix this.  The first is to check gene names for typos or use of less common colloquial names.  The second is to use the missNames and fixNames functions in this pacakge to fill in the missing ENSG identifiers.  Note: googling the colloquial gene name and 'gene cards' is an excellent way to find an ENSG id.  Genecards does an exceptional job of cataloging alternative colloquial names. 
#'
#'@export
 
missNames <- function(mgl){
	x = length(mgl)
	tmp <- list(); for (i in 1:x) {tmp[[i]] <- mgl[[i]][[1]]}
	lapply(tmp, class) -> check
	names(mgl)[which(check == 'character')] -> missing
	return(missing)
}

