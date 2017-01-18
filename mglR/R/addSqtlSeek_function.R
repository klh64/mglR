#' Add results from SqtlSeek splicingQTL algorithm to list
#'
#' \code{addSqtlSeek} returns an 'mgl' list with the fifteenth element as a dataframe containing results from SqtlSeek splicingQTL algorithm from GTEx; see \url{http://www.gtexportal.org}.  
#'
#' This gives the SqtlSeek algorithm results for each gene.  Data is downloaded from the GTEx website \url{http://www.gtexportal.org/static/datasets/gtex_analysis_pilot_v3/splicing_qtls_sqtls/sQTLs-sQTLseeker-merged.tgz}. 
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param download A logical vector indicating if the data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_sqtlSeek.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' myMgl <- addSqtlSeek(myMgl)
#'
#'@export
 
addSqtlSeek <- function(mgl, download = TRUE, saveDownload = FALSE, fpsource = "./"){

# Getting data
if (download == TRUE){	
	# Download data
	download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_pilot_v3/splicing_qtls_sqtls/sQTLs-sQTLseeker-merged.tgz', 'sQTLs-sQTLseeker-merged.tgz')
	untar('sQTLs-sQTLseeker-merged.tgz')	
	files <- list.files(getwd(), pattern = 'FDR05', full.names = T)
	
	# Read into R
	seek <- list()
	for(x in 1:length(files)){
		seek[[x]] <- read.table(files[x], header = T, stringsAsFactors = F,  sep = '\t', quote = "")
	}
	names(seek) <- stringr::word(stringr::word(files, 2, sep = stringr::fixed('sQTLs-')), 1, sep = stringr::fixed('-sQTLseekeR'))	
	
	# save
	if(saveDownload == TRUE){
		save(seek, file = 'RawData_sqtlSeek.RData')
	}
	# Deleting files
	unlink("sQTLs-sQTLseeker-merged.tgz")
	lapply(files,unlink)
}

if (download == FALSE){	
	load(paste(fpsource, 'RawData_sqtlSeek.RData', sep = ""))	
}

# add to mgl
x = length(mgl)

for (i in 1:x) {mgl[[i]][[15]] <- lapply(as.list(1:9), function(a) seek[[a]][grep(mgl[[i]][[1]][2][1,], seek[[a]][,1]),c(2,3,7,8)])}

# add names
for (i in 1:x) {names(mgl[[i]][[15]]) <- names(seek)}
	
	return(mgl)
}


