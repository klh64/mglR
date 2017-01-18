#' Add gene ontology to list
#'
#' \code{addGo} returns gene ontlology (GO) terms for the gene using biomaRt as the fifth element in the 'mgl' list.  
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addGo(myMgl)}
#'
#'@export
 
addGo <- function(mgl){
	mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', host = 'feb2014.archive.ensembl.org', 
		path = '/biomart/martservice', dataset = 'hsapiens_gene_ensembl')
	
	if (length(mgl) > 5){		
	as.numeric(stringr::word((length(mgl)/5), 1, sep = stringr::fixed('.'))) -> y
	for(x in 0:(y-1)){
		z <- (5*x)+1
		p <- z+4
		for (i in z:p) {mgl[[i]][[5]] <- biomaRt::getBM(attributes = c('name_1006'), filters = 
			'ensembl_gene_id', values = mgl[[i]][[1]][1,2], mart = mart)}}
		for(i in ((5*y)+1):((5*y)+length(mgl) - (y*5))) {mgl[[i]][[5]] <- biomaRt::getBM(attributes = 
			c('name_1006'), filters = 'ensembl_gene_id', values = mgl[[i]][[1]][1,2], mart = mart)}}
			
	if (length(mgl) <= 5){
		for(i in 1:length(mgl)){
		mgl[[i]][[5]] <- biomaRt::getBM(attributes = c('name_1006'), filters = 
			'ensembl_gene_id', values = mgl[[i]][[1]][1,2], mart = mart)}
	}
	return(mgl)
}


