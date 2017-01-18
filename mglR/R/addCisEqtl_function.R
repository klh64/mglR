#' Add cis eQTL data to list
#'
#' \code{addCisEqtl} returns an 'mgl' list with the fourteenth element as a list with each element being a dataframe with cis eQTL data from GTEx for a given tissue for the gene of interest.    
#'
#' This gives basic information on tissue specific eQTLs as reported by GTEx for the gene of interest.  It pulls eQTL data based on the ENSG identifier.  Data is downloaded from GTEx \url{gtexportal.org/}.  
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param download A logical vector indicating if the data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'RawData_eQTL.RData'
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addCisEqtl(myMgl, download = TRUE, saveDownload = FALSE)}
#'
#'@export
#'@importFrom utils download.file untar read.table

addCisEqtl <- function(mgl, download = TRUE, saveDownload = FALSE, fpsource = "./"){

# make vector of ensgs with CHR at the start to capture header
ensgs <- unlist(lapply(mgl, function(x) x[[1]][1,2]))

# if NAs exist remove them
if (sum(is.na(ensgs)) != 0) {
	ensgs <- ensgs[-c(which(is.na(ensgs) == T))]
}

####

if (download == TRUE){	
# Download data
temp <- tempfile()
download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_eQTL.tar', temp)
	
# Unzip and untar
untar(temp, exdir = ".")

# compile code
files <- list.files(paste(getwd(), '/GTEx_Analysis_v6p_eQTL', sep = ""), pattern = 'signif_snpgene_pairs')
code <- lapply(files, function(x) paste('zgrep ', ensgs, ' ', getwd(), '/GTEx_Analysis_v6p_eQTL/', x, ' >> ', getwd(), '/', x, sep = ""))

unlink(temp)

}

if (download == FALSE){	
# compile code
	files <- list.files(paste(fpsource, '/GTEx_Analysis_v6p_eQTL', sep = ""), pattern = 'signif_snpgene_pairs')
	code <- lapply(files, function(x) paste('zgrep ', ensgs, ' ', fpsource, 'GTEx_Analysis_v6p_eQTL/', x, ' >> ', getwd(), '/', x, sep = ""))	
}

# run code via unix
code <- unlist(code)
for(x in 1:length(code)){
system(code[[x]])
}

# read in eQTLs
eqtls <- as.list(1:length(files))
names(eqtls) <- stringr::word(files, 1, sep = stringr::fixed('_Analysis'))
for(x in 1:length(files)){
	tryCatch({
	eqtls[[x]] <- read.table(files[x], stringsAsFactors = F, header = F, sep = '\t')}, error = function(e){})
}

# split tissue name to separate column and add colnames
for(x in 1:length(eqtls)){
	tryCatch({
eqtls[[x]] <- cbind(names(eqtls[x]), eqtls[[x]])
colnames(eqtls[[x]]) <- c('tissue','variant_id', 'gene_id', 'tss_distance', 'pval_nominal', 'slope', 'slope_se', 'slope_fpkm', 'slope_fpkm_se', 'pval_nominal_threshold', 'min_pval_nominal', 'pval_beta')}, error = function(e){})
}

## add to mgl
for (i in 1:length(mgl)) {
	mgl[[i]][[14]] <- as.list(1:length(eqtls))
	names(mgl[[i]][[14]]) <- names(eqtls)
	for(x in 1:length(eqtls)){
		tryCatch({
			mgl[[i]][[14]][[x]] <- eqtls[[x]][grep(mgl[[i]][[1]][2][1,], eqtls[[x]][,3]),]}, error = function(e){})
}
}

# Delete files
lapply(files, function(x) unlink(x)) -> remove

if(download == TRUE){
if(saveDownload == FALSE){
	unlink("./GTEx_Analysis_v6p_eQTL", recursive = TRUE)
}	
}

return(mgl)

}

