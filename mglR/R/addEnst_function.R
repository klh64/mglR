#' Add transcript information to list
#'
#' \code{addEnst} returns an 'mgl' list with the second element containing information about the trascripts for the gene of interest including external transcript name, ensembl transcript name, and biotype (an Ensembl defined classification of transcript type: e.g. non-coding, protein coding; see \url{http://www.ensembl.org/Help/Faq?id=468})  
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addEnst(myMgl)}
#'
#'@export
 
addEnst <- function(mgl){
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', 
		path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')
	for (i in 1:length(mgl)) {
		mgl[[i]][[2]] <- biomaRt::getBM(attributes = c('external_transcript_id', 
			'ensembl_transcript_id', 'transcript_biotype'), filters = 'ensembl_gene_id', values = 
			mgl[[i]][[1]][1,2], mart = mart)}
	return(mgl)
}

