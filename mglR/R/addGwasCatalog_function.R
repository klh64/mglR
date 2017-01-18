#' Add GWAS data from NHGRI-EBI GWAS Catalog to list
#'
#' \code{addGwasCatalog} returns a list with the eighteenth element as a dataframe with GWAS data from the NHGRI-EBI GWAS Catalog.    
#'
#' This gives basic information on trait associated variants as reported by the NHGRI-EBI GWAS Catalog for the gene of interest.  It pulls data based on the gene name and position (i.e. any SNP that falls between the start and stop position).  Note a wider range can be taken using the \emph{range} flag.   
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param range Number indicating distance upstream of start and downstream of stop that should be used for filtering SNPs.
#'
#' @param download A logical vector indicating if most recent version of the NHGRI-EBI GWAS Catalog should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_GwasCatalog.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addGwasCatalog(myMgl, range = 0, download = TRUE, 
#'	saveDownload = TRUE, fpsource = "./")}
#'
#'@export
#'@importFrom utils download.file read.table
 

addGwasCatalog <- function(mgl, range = 0, download = TRUE, saveDownload = FALSE, fpsource = "./"){

# stop process if location (element 3) has not been filled in
if (unique(unlist(lapply(mgl, function(x) class(x[[3]])))) == 'integer') stop('location data not available. see function addLoc') 

# to download the most recent version
if (download == TRUE){

# Download data
temp <- tempfile()
download.file('http://www.ebi.ac.uk/gwas/api/search/downloads/full', temp)

# Read in file
read.table(temp, header = T, stringsAsFactors = F, sep = '\t', quote = "", comment.char = "", as.is = TRUE) -> gwas

if (saveDownload == TRUE){
# Save files
save(gwas, file = 'RawData_GwasCatalog.RData')
	}
# Delete files
unlink(temp)
}

if (download == FALSE){
load(paste(fpsource, "RawData_GwasCatalog.RData", sep = ""))}	
	

# add to mgl
for (i in 1:length(mgl)) {mgl[[i]][[18]] <- gwas[unique(c(
	# gene name: "REPORTED.GENE.S."
	grep(names(mgl)[i], gwas[,14]),
	# gene name: "MAPPED_GENE" 
	grep(names(mgl)[i], gwas[,15]), 
	# gene name: "UPSTREAM_GENE_ID"
	grep(names(mgl)[i], gwas[,16]),
	# gene name: "DOWNSTREAM_GENE_ID" 
	grep(names(mgl)[i], gwas[,17]),
	# gene position: "CHR_ID" & "CHR_POS"
	which(gwas[,12] == mgl[[i]][[3]][1,1] & gwas[,13] <= (mgl[[i]][[3]][1,3] + range) & gwas[,13] >= (mgl[[i]][[3]][1,2] - range))))
	,]
	rownames(mgl[[i]][[18]]) <- NULL}
	
return(mgl)

}

