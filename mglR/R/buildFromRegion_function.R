#' Build empty list using genomic region.
#'
#' \code{buildFromRegion} returns an empty 'mgl' list given a genomic region.
#'
#' This is one of three functions that can be used to set up the list structure.  It starts
#'   with a position range.
#'
#' @family Build list
#'
#' @param chr Chromosome
#' @param start Start position
#' @param stop	Stop position
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#'
#'@export

buildFromRegion <- function(chr, start, stop){
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', 
		path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')	
	
	nm <- biomaRt::getBM(attributes = c('external_gene_id', 'ensembl_gene_id'), filters = 
		c('chromosome_name', 'start', 'end'), values = list(chr, start, stop), mart = mart)
	
	x = length(nm[,1])
	mgl <- as.list(1:x)
	for(i in 1:x){
 		mgl[[i]] <- as.list(1:20)
 		names(mgl[[i]]) <- c('name', 'enst', 'location', 'antisense', 'go', 'pubmed', 'gtex.normalized', 'gtex.gene.counts', 'gtex.transcript.counts', 'gtex.gene.rpkm', 'gtex.transcript.rpkm', 'dnase','transEqtls', 'cisEqtls', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'gwasCatalog', 'grasp', 'aei')}
	names(mgl) <- nm[,1]
	
	bm <- biomaRt::getBM(attributes = c('external_gene_id','ensembl_gene_id', 
			'description'), filters = 'ensembl_gene_id', values = nm[,2], mart = mart)
			
	m <- which(nm[,2] %in% bm[,2]) 
	k <- 1:length(mgl)
	n <- k[-which(k %in% m)]
	
	if(length(n) > 0){
	warning('Gene names missing: ', paste(names(mgl)[n], collapse = ", "))
	for(x in 1:length(n)){
		mgl[[n[x]]][[1]] <- 'NA'
		}
	}
	
	for(x in 1:length(m)){
		mgl[[x]][[1]] <- bm[which(nm[x,2] == bm[,2]),] 
		}
		
	return(mgl)
}


