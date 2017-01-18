#' Build empty list using gene names.
#'
#' \code{buildFromNames} returns an empty 'mgl' list given gene names.
#'
#' This is one of three functions that can be used to set up the list structure.  It starts with gene names.  
#'
#' @family Build list
#'
#' @param nm A character vector with gene names
#' 
#' @examples
#' \dontrun{buildFromNames(nm = c('CETP', 'ABCA1', 'APOB')) -> myMgl}
#'
#' @section Warning:
#' If a gene is not found, a warning message will appear: 'Gene names missing: ...'.  This must be corrected or no other elements can be filled in as the remaining elements all build off of the disambigous ENSG identifier.  There are two strategies to fix this.  The first is to check gene names for typos or use of less common colloquial names.  The second is to use the missNames and fixNames functions in this pacakge to fill in the missing ENSG identifiers.  Note: googling the colloquial gene name and 'gene cards' is an excellent way to find an ENSG id.  Genecards does an exceptional job of cataloging alternative colloquial names. 
#'
#'@export

buildFromNames <- function(nm){	
	x = length(nm)
	mgl <- as.list(1:x)
	for(i in 1:x){
 		mgl[[i]] <- as.list(1:20)
 		names(mgl[[i]]) <- c('name', 'enst', 'location', 'antisense', 'go', 'pubmed', 'gtex.normalized', 'gtex.gene.counts', 'gtex.transcript.counts', 'gtex.gene.rpkm', 'gtex.transcript.rpkm', 'dnase','transEqtls', 'cisEqtls', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'gwasCatalog', 'grasp', 'aei')}
	names(mgl) <- nm
	
mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')
		
bm <- biomaRt::getBM(attributes = c('uniprot_genename','external_gene_id','ensembl_gene_id', 'description'), filters = 'uniprot_genename', values = names(mgl), mart = mart)
			
	m <- which(names(mgl) %in% bm[,1]) 
	k <- 1:length(mgl)
	n <- k[-which(k %in% m)]
	
	if(length(n) != 0){
	warning('Gene names missing: ', paste(names(mgl)[n], collapse = ", "))
	for(x in 1:length(n)){
		mgl[[n[x]]][[1]] <- 'NA'
		}
	}
	
	for(x in 1:length(m)){
		mgl[[m[x]]][[1]] <- bm[which(bm[,1] == names(mgl)[m[x]]),2:4] 
		}
		
	return(mgl)
}


