#' Add antisense transcript information to list
#'
#' \code{addAntisense} returns a list with the fourth element as a dataframe containing external gene name, ensembl transcript name, biotype (an Ensembl defined classification of transcript type: e.g. non-coding, protein coding; see \url{http://www.ensembl.org/Help/Faq?id=468}), chromosome name, start position, stop position, and strand for those transcripts antisense to the gene of interest.  
#'
#' This gives basic information on transcripts that are antisense to the gene of interest.  It works by defining the gene region, selecting the opposite strand, and quering biomaRt for transcripts.  It will thus detect any transcript in a position antisense to the gene of interest regardless of its name or coding status.  
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addAntisense(myMgl)}
#'
#'@export
#'@import biomaRt
 
addAntisense <- function(mgl){
	
# stop process if location (element 3) has not been filled in
if (unique(unlist(lapply(mgl, function(x) class(x[[3]])))) == 'integer') stop('location data not available. see function addLoc') 

# load biomart
mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')

# add to list
for (i in 1:length(mgl)) {
	mgl[[i]][[4]] <- biomaRt::getBM(attributes = c('external_gene_id', 'ensembl_transcript_id', 
	'transcript_biotype','chromosome_name','start_position', 'end_position', 'strand'), 
	filters = 'chromosomal_region', values = paste(unique(mgl[[i]][[3]][,1]), ':', 
	min(mgl[[i]][[3]][,2]), ":", max(mgl[[i]][[3]][,3]), ":", 
	(unique(mgl[[i]][[3]][,4])*-1), sep = ""), mart = mart)
}
	return(mgl)
}

