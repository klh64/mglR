#' Add gene location information to list
#'
#' \code{addLoc} returns an 'mgl' list with the third element as a dataframe with position information for the gene of interest: chromosome name, start position, stop position, and strand for each gene.
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addLoc(myMgl)}
#'
#'@export
 
addLoc <- function(mgl){
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', 
		path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')
	for (i in 1:length(mgl)) {
		mgl[[i]][[3]] <- biomaRt::getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'strand'), filters = 'ensembl_gene_id', values = mgl[[i]][[1]][1,2], mart = mart)}
	return(mgl)
}
