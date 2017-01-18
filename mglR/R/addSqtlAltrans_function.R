#' Add results from Altrans splicingQTL algorithm to list
#'
#' \code{addSqtlAltrans} returns an 'mgl' list with the sixteenth element as a dataframe containing results from Altrans splicingQTL algorithm from GTEx; see \url{http://www.gtexportal.org}  
#'
#' This gives the Altrans algorithm results for each gene.  Data is downloaded from the GTEx website \url{http://www.gtexportal.org/static/datasets/gtex_analysis_pilot_v3/splicing_qtls_sqtls/Altrans_FDR05_bestPerLink.tgz}.
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param download A logical vector indicating if the data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_sqtlAltrans.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' myMgl <- addSqtlAltrans(myMgl)
#'
#'@export
 
addSqtlAltrans <- function(mgl, download = TRUE, saveDownload = FALSE, fpsource = "./"){

# Getting data
if (download == TRUE){	
	# Download data
	download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_pilot_v3/splicing_qtls_sqtls/Altrans_FDR05_bestPerLink.tgz', 'Altrans_FDR05_bestPerLink.tgz')
	untar('Altrans_FDR05_bestPerLink.tgz')	
	files <- list.files(paste(getwd(), '/Altrans_FDR05_bestPerLink', sep = ""), full.names = T)
	
	# Read into R
	altrans <- list()
	for(x in 1:length(files)){
		altrans[[x]] <- read.table(files[x], header = T, stringsAsFactors = F, comment.char = "")
	}
	names(altrans) <- stringr::word(stringr::word(files, 2, sep = stringr::fixed('bestPerLink/')), 1, sep = stringr::fixed('.Altrans'))	
	
	# save
	if(saveDownload == TRUE){
		save(altrans, file = 'RawData_sqtlAltrans.RData')
	}
	# Deleting files
	unlink("Altrans_FDR05_bestPerLink.tgz")
	unlink("Altrans_FDR05_bestPerLink", recursive = T)
}

if (download == FALSE){	
	load(paste(fpsource, 'RawData_sqtlAltrans.RData', sep = ""))	
}

# add to mgl
	for (i in 1:length(mgl)) {mgl[[i]][[16]] <- lapply(as.list(1:9), function(a) 
		altrans[[a]][grep(mgl[[i]][[1]][2][1,], altrans[[a]][,2]),c(1,6)])}

# with names
	for (i in 1:length(mgl)) {names(mgl[[i]][[16]]) <- names(altrans)}
	
	return(mgl)
}


