#' Add trans eQTL data to list
#'
#' \code{addTransEqtl} returns an 'mgl' list with the thirteenth element as a list with each element being a dataframe with trans eQTL data from GTEx for a given tissue for the gene of interest.    
#'
#' This gives basic information on tissue specific eQTLs as reported by GTEx for the gene of interest.  It pulls eQTL data based on the ENSG identifier.  Data is downloaded from GTEx \url{gtexportal.org/}.  
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param download A logical vector indicating if the data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_transEqtls.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' myMgl <- addTransEqtl(myMgl)
#'
#'@export

addTransEqtl <- function(mgl, download = T, saveDownload = F, fpsource = "./"){

# Getting data
if (download == TRUE){	
	# Download data
	download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_trans_eQTLs.xlsx', 'GTEx_Analysis_v6p_trans_eQTLs.xlsx')
#	download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_intrachromosomal_eQTLs.xlsx', 'GTEx_Analysis_v6p_intrachromosomal_eQTLs.xlsx')	
	# Read into R
	trans <- as.list(1:4)	
	names(trans) <- c('genome-wide', 'LD-pruned', 'cis-eVariants', 'trait-associated')
	for(x in 1:4){
		trans[[x]] <- gdata::read.xls("GTEx_Analysis_v6p_trans_eQTLs.xlsx", sheet = x, header = TRUE)
		trans[[x]] <- cbind(names(trans)[x], trans[[x]])
	}
#	gdata::read.xls("GTEx_Analysis_v6p_intrachromosomal_eQTLs.xlsx", sheet = 1, header = TRUE) -> t1
#	gdata::read.xls("GTEx_Analysis_v6p_intrachromosomal_eQTLs.xlsx", sheet = 2, header = TRUE) -> t2
	
	trans[[1]] <- cbind(trans[[1]], NA, NA, NA, NA, NA)
	trans[[2]] <- cbind(trans[[2]], NA, NA, NA, NA, NA)
	trans[[3]] <- cbind(trans[[3]][,c(1:7, 10:12, 8:9)], NA, NA, NA)
	trans[[4]] <- cbind(trans[[4]][, c(1:10)], NA, NA, trans[[4]][,c(11:13)])
	for(x in 1:4){
		colnames(trans[[x]]) <- c('DiscoveryMode', 'rs_id_dbSNP142_GRCh37p13', 'variant_id', 'snp_chr', 'snp_pos', 'gene_id', 'gene_name', 'tissue_site_detail', 'pvalue', 'FDR', 'cis.eGene_id', 'cis.eGene_name', 'GWAS_pvalue', 'trait', 'PMID')
	}
	do.call(rbind, trans) -> trans
	rownames(trans) <- NULL
	tissues <- as.character(unique(trans[,8]))
	lapply(tissues, function(x) trans[which(trans[,8] == x),]) -> trans
	names(trans) <- tissues
	
	# save
	if(saveDownload == TRUE){
		save(trans, file = 'RawData_transEqtls.RData')
	}
	# Deleting files
	unlink("GTEx_Analysis_v6p_trans_eQTLs.xlsx")
}

if (download == FALSE){	
	load(paste(fpsource, 'RawData_transEqtls.RData', sep = ""))	
}

## add to mgl
for (i in 1:length(mgl)) {
	mgl[[i]][[13]] <- list()
	for(x in 1:length(trans)){
		tryCatch({
			mgl[[i]][[13]][[x]] <- trans[[x]][grep(mgl[[i]][[1]][2][1,], trans[[x]][,6]),]}, error = function(e){})
}
names(mgl[[i]][[13]]) <- names(trans)
}

return(mgl)

}

