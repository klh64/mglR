#' Build empty list using ENSG gene ids.
#'
#' \code{buildFromEnsgs} returns an empty 'mgl' list given ensembl ENSG gene ids.
#'
#' This is one of three functions that can be used to set up the list structure.  It starts
#'   with ENSG gene ids.  
#'
#' @family Build list
#'
#' @param ensg A character vector with ENSG gene ids
#' 
#' @examples
#' \dontrun{buildFromEnsgs(ensg = 
#'	c('ENSG00000087237', 'ENSG00000165029', 'ENSG00000084674')) -> myMgl}
#'
#'@export


buildFromEnsgs <- function(ensg){	
	x = length(ensg)
	mgl <- as.list(1:x)
	for(i in 1:x){
 		mgl[[i]] <- as.list(1:20)
 		names(mgl[[i]]) <- c('name', 'enst', 'location', 'antisense', 'go', 'pubmed', 'gtex.normalized', 'gtex.gene.counts', 'gtex.transcript.counts', 'gtex.gene.rpkm', 'gtex.transcript.rpkm', 'dnase','transEqtls', 'cisEqtls', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'gwasCatalog', 'grasp', 'aei')}
	names(mgl) <- ensg
	
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', 
		path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')
		
	bm <- biomaRt::getBM(attributes = c('external_gene_id','ensembl_gene_id','description'), filters = 'ensembl_gene_id', values = 
			names(mgl), mart = mart)
			

	bm <- bm[match(ensg, bm[,2]),]
	
	for(x in 1:length(ensg)){
		mgl[[x]][[1]] <- bm[which(bm[,2] == ensg[x]),]
		}
	
	names(mgl) <- bm[,1]

		
	return(mgl)
}


