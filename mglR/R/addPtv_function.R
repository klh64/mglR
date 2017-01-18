#' Add results from protein truncating variant algorithm made available by GTEx
#'
#' \code{addPtv} returns an 'mgl' list with the seventeenth element as a dataframe containing results from psiQTL algorithm; see \url{http://www.gtexportal.org/home/}.  
#'
#' This gives the protein truncating variant results for each gene.  Data was originally downloaded from the GTEx website \url{http://www.gtexportal.org/static/datasets/gtex_analysis_pilot_v3/ptv_data/gtex_psiqtls.zip}.
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param range Number indicating distance upstream of start and downstream of stop that should be used for filtering SNPs.
#'
#' @param download A logical vector indicating if the data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_Ptv.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' myMgl <- addPtv(myMgl, range = 0)
#'
#'@export
#'@importFrom utils unzip
 
addPtv <- function(mgl, range = 0, download = T, saveDownload = F, fpsource = "./"){
	
# stop process if location (element 3) has not been filled in
if (unique(unlist(lapply(mgl, function(x) class(x[[3]])))) == 'integer') stop('location data not available. see function addLoc') 

# Getting data
if (download == TRUE){	
	# Download data
	download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_pilot_v3/ptv_data/gtex_psiqtls.zip', 'gtex_psiqtls.zip')
	unzip('gtex_psiqtls.zip')	
	files <- list.files(getwd(), pattern = "Assoc-total", full.names = T)
	
	# Read into R
	ptv <- list()
	for(x in 1:length(files)){
		ptv[[x]] <- read.table(files[x], header = F, stringsAsFactors = F, comment.char = "")
		colnames(ptv[[x]]) <- c('identifier of snp', 'alternative identifier of SNP', 'name of exon that snp belongs to', 'chr_start_end', 'chromosome of SNP', 'chromosome of exon', 'location of SNP', 'location of middle point of exon', 'distance between 8 and 9', 'spearman correlation', 'pvalue', '-log10(pvalue)')
	}
	names(ptv) <- stringr::word(stringr::word(files, 2, sep = stringr::fixed('Assoc-total.')), 1, sep = stringr::fixed('.txt'))	
	
	# save
	if(saveDownload == TRUE){
		save(ptv, file = 'RawData_Ptv.RData')
	}
	# Deleting files
	unlink("gtex_psiqtls.zip")
	unlink("__MACOSX", recursive = T)
	lapply(files,unlink)
}

if (download == FALSE){	
	load(paste(fpsource, 'RawData_Ptv.RData', sep = ""))	
}

# add to mgl
for (i in 1:length(mgl)) {
	mgl[[i]][[17]] <- 
# gene position: "CHR_ID" & "CHR_POS"
lapply(ptv, function(a) a[which(a[,5] == mgl[[i]][[3]][1,1] & a[,7] <= (mgl[[i]][[3]][1,3] + range) & a[,7] >= (mgl[[i]][[3]][1,2] - range)),])
	}
# with names
for (i in 1:length(mgl)) {names(mgl[[i]][[17]]) <- names(ptv)}
	
return(mgl)
}


