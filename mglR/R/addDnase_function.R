#' Add data from Maurano et. al. Nat Genet. 2015 regarding allele specific DNAse hypersensitivty.
#'
#' \code{addDnase} returns an 'mgl' list with the twelfth element as a dataframe containing results from Maurano et. al. Nat Genet. 2015  
#'
#' This gives the results for each gene, subsetting the data to reflect those SNPs that fall within the gene region.  Data was originally downloaded from \url{http://www.nature.com/ng/journal/v47/n12/extref/ng.3432-S5.txt}.  
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param range Number indicating distance upstream of start and downstream of stop that should be used for filtering SNPs.
#'
#' @param download A logical vector indicating if the data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_Dnase.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' myMgl <- addDnase(myMgl, range = 0)
#'
#' @section Warning:
#' All results will be added.  To view only SNPs with a significant imbalance use \code{\link{makeDnaseSig}}
#'
#'@export

addDnase <- function(mgl, range = 0, download = TRUE, saveDownload = FALSE, fpsource = "./"){
	
# stop process if location (element 3) has not been filled in
if (unique(unlist(lapply(mgl, function(x) class(x[[3]])))) == 'integer') stop('location data not available. see function addLoc') 

# Getting data
if (download == TRUE){	
	# Download data
	temp <- tempfile()
	download.file('http://www.nature.com/ng/journal/v47/n12/extref/ng.3432-S5.txt', temp)	
	# Read into R
	read.table(temp, stringsAsFactors =F, header = T) -> dnase
	if(saveDownload == TRUE){
		save(dnase, file = 'RawData_Dnase.RData')
	}
unlink(temp)
}

if (download == FALSE){	
	load(paste(fpsource, 'RawData_Dnase.RData', sep = ""))	
}

# Adding to mgl
for(i in 1:length(mgl)) {
mgl[[i]][[12]] <- rbind(dnase[which(dnase[,1] == paste('chr', mgl[[i]][[3]][1,1], sep = '') & dnase[,2] >= (mgl[[i]][[3]][1,2] - range) & dnase[,3] <= (mgl[[i]][[3]][1,3] + range)),])}
		
for(i in 1:length(mgl)) {rownames(mgl[[i]][[12]]) <- NULL}
	
	return(mgl)
}


