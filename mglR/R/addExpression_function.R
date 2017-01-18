#' Add RNAsequencing expression data list
#'
#' \code{addExpression} returns a list with the eighth through the eleventh elements as lists containing dataframes with tissue specific expression data from GTEx; see \url{http://www.gtexportal.org}. Element 8 is a list of reads by gene, element 9 of reads by transcript, element 10 of RPKM by gene, and element 11 of RPKM by transcript.
#'
#' This gives tissue specific RNAsequencing data as reported by GTEx for the gene of interest.  It pulls data based on the gene name.  Each element of the list is a separate tissue.
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param download A logical vector indicating if expression data should be downloaded.
#'
#' @param saveDownload A logical vector indicating if the data should be saved 
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @param normalized A logical flag indicating if normalized gene expression data should be incorporated into the list.  Note: this is both the gene expression values and covariates used by GTEx to calculate eQTLs
#'
#' @param xpGeneReads A logical flag indicating if gene read expression data should be incorporated into the list
#'
#' @param xpTranscriptReads A logical flag indicating if transcript read expression data should be incorporated into the list
#'
#' @param xpGeneRpkm A logical flag indicating if gene RPKM expression data should be incorporated into the list
#'
#' @param xpTranscriptRpkm A logical flag indicating if transcript RPKM expression data should be incorporated into the list
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addEnst(myMgl)}
#' \dontrun{myMgl <- addExpression(myMgl, download = TRUE, saveDownload = FALSE, 
#'	fpsource = './', normalized = FALSE, xpGeneReads = TRUE, xpTranscriptReads = FALSE, 
#'	xpGeneRpkm = FALSE, xpTranscriptRpkm = FALSE)}
#'
#'@export
#'@importFrom utils download.file untar read.table 
#'@importFrom stats cov
 
addExpression <- function(mgl, download = TRUE, saveDownload = FALSE, fpsource = './', normalized = TRUE, xpGeneReads = TRUE, xpTranscriptReads = TRUE, xpGeneRpkm = TRUE, xpTranscriptRpkm = TRUE){

# Genes
if(xpGeneReads == TRUE | xpGeneRpkm == TRUE | normalized == TRUE){
# make vector of gene names
#genes <- c('GTEX', names(mgl))
genes <- c('GTEX', do.call(rbind,listSeparate(mgl,1))[,2])
# if NAs exist remove them
if (sum(is.na(genes)) != 0) {
	genes <- genes[-c(which(is.na(genes) == T))]}
}

# Transcripts
if(xpTranscriptReads == TRUE | xpTranscriptRpkm == T){
# stop process if asking for transcript expression data, but transcripts (element 2) have not yet been added
if (unique(unlist(lapply(mgl, function(x) class(x[[2]])))) == 'integer') stop('transcript data not available. see function addEnst') 
transcripts <- list(); for (x in 1:length(mgl)){transcripts[[x]] <- mgl[[x]][[2]][,2]}
transcripts <- c('GTEX', unique(unlist(transcripts)))
# if NAs exist remove them
if (sum(is.na(transcripts)) != 0) {
	transcripts <- transcripts[-c(which(is.na(transcripts) == T))]}
}

###

if(normalized == TRUE){
	
	if(download == TRUE){
		download.file ('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_eQTL_expression_matrices.tar', 'GTEx_Analysis_v6p_eQTL_expression_matrices.tar')
		download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/single_tissue_eqtl_data/GTEx_Analysis_v6p_eQTL_covariates.tar.gz', 'GTEx_Analysis_v6p_eQTL_covariates.tar.gz')
		untar('GTEx_Analysis_v6p_eQTL_expression_matrices.tar')
		untar('GTEx_Analysis_v6p_eQTL_covariates.tar.gz')
		files <- list.files(paste(getwd(), '/GTEx_Analysis_v6p_eQTL_expression_matrices', sep = ""))
		files <- files[-c(grep('tbi', files))]
		code <- lapply(files, function(x) paste('zgrep ', genes, ' ', getwd(), '/GTEx_Analysis_v6p_eQTL_expression_matrices/', x, ' >> ', getwd(), '/', x, sep = ""))
		cfiles <- list.files(paste(getwd(), '/GTEx_Analysis_v6p_eQTL_covariates', sep = ""), full.names = T)
}
		
	if (download == FALSE){
		files <- list.files(paste(fpsource, '/GTEx_Analysis_v6p_eQTL_expression_matrices', sep = ""))
		files <- files[-c(grep('tbi', files))]
		code <- lapply(files, function(x) paste('zgrep ', genes, ' ', fpsource, '/GTEx_Analysis_v6p_eQTL_expression_matrices/', x, ' >> ', getwd(), '/', x, sep = ""))
		cfiles <- list.files(paste(fpsource, '/GTEx_Analysis_v6p_eQTL_covariates', sep = ""), full.names = T)
		}

# run code via unix
code <-unlist(code)
for(x in 1:length(code)){
	system(code[[x]])
}

# read in data
dane <- list()
for(x in 1:length(files)){
	tryCatch({
	dane[[x]] <- read.table(files[x], stringsAsFactors = F, header = T, sep = '\t', comment.char = "")}, error = function(e){})
	rownames(dane[[x]]) <- dane[[x]][,4]
	dane[[x]] <- dane[[x]][,-c(1:4)]
}
names(dane) <- stringr::word(files, 1, sep = stringr::fixed('_Analysis'))

# read in covariates
cov <- list()
for(x in 1:length(cfiles)){
	cov[[x]] <- read.table(cfiles[x], stringsAsFactors = F, header = T, sep = '\t', comment.char = "")
	rownames(cov[[x]]) <- cov[[x]][,1]
	cov[[x]] <- cov[[x]][,-1]
}
names(cov) <- stringr::word(stringr::word(cfiles, 2, sep = stringr::fixed('GTEx_Analysis_v6p_eQTL_covariates/')), 1, sep = stringr::fixed('_Analysis'))

for (i in 1:length(mgl)) {
	mgl[[i]][[7]] <- list()
	for(x in 1:length(dane)){
			mgl[[i]][[7]][[x]] <- rbind(
			dane[[x]][grep(mgl[[i]][[1]][2][1,], rownames(dane[[x]])),],
			cov[[x]][,match(colnames(dane[[x]]), colnames(cov[[x]]))]
			)
			}
names(mgl[[i]][[7]]) <- names(dane)
}

# delete files as appropriate
lapply(files, function(x) unlink(x)) -> remove

}

###

xp <- list()
txp <- list()
# Gene Reads
if(xpGeneReads == TRUE){
if (download == TRUE){
download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz', 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz')	
code <- lapply(genes, function(x) paste('zgrep ', x, ' GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz >> xpGeneReads.txt', sep = "" ))
}
if (download == FALSE){
code <- lapply(genes, function(x) paste('zgrep ', x, ' ', fpsource, 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz >> xpGeneReads.txt', sep = "" ))
}
for(x in 1:length(code)){
system(code[[x]])}
read.table('xpGeneReads.txt', header = T, stringsAsFactors = F) -> xp[[1]]
}
# Transcript Reads
if(xpTranscriptReads == TRUE){
if (download == TRUE){
download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz', 'GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz')	
code <- lapply(transcripts, function(x) paste('zgrep ', x, ' GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz >> xpTranscriptReads.txt', sep = "" ))
}
if (download == FALSE){
code <- lapply(genes, function(x) paste('zgrep ', x, ' ', fpsource, 'GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz >> xpTranscriptReads.txt', sep = "" ))
}
for(x in 1:length(code)){
system(code[[x]])}
read.table('xpTranscriptReads.txt', header = T, stringsAsFactors = F) -> xp[[2]]
}
# Gene RPKM
if(xpGeneRpkm == TRUE){
if (download == TRUE){
download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6p/rna_seq_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz', 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz')	
code <- lapply(genes, function(x) paste('zgrep ', x, ' GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz >> xpGeneRpkm.txt', sep = "" ))
}
if (download == FALSE){
code <- lapply(genes, function(x) paste('zgrep ', x, ' ', fpsource, 'GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz >> xpGeneRpkm.txt', sep = "" ))
}
for(x in 1:length(code)){
system(code[[x]])}
read.table('xpGeneRpkm.txt', header = T, stringsAsFactors = F) -> xp[[3]]
}
# Transcript RPKM
if(xpTranscriptRpkm == TRUE){
if (download == TRUE){
download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt.gz', 'GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt.gz')	
code <- lapply(transcripts, function(x) paste('zgrep ', x, ' GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt.gz >> xpTranscriptRpkm.txt', sep = "" ))
}
if (download == FALSE){
code <- lapply(genes, function(x) paste('zgrep ', x, ' ', fpsource, 'GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt.gz >> xpTranscriptRpkm.txt', sep = "" ))
}
for(x in 1:length(code)){
system(code[[x]])}
read.table('xpTranscriptRpkm.txt', header = T, stringsAsFactors = F) -> xp[[4]]
}

if(xpGeneReads == TRUE | xpGeneRpkm == TRUE | xpTranscriptReads == TRUE | xpTranscriptRpkm == T){
# load key to separate ids into appropriate tissue groups
if (download == TRUE){	
# Download sample attributes
download.file('http://www.gtexportal.org/static/datasets/gtex_analysis_v6/annotations/GTEx_Data_V6_Annotations_SampleAttributesDS.txt', 'GTEx_Data_V6_Annotations_SampleAttributesDS.txt')
}
read.table('GTEx_Data_V6_Annotations_SampleAttributesDS.txt', header = T, stringsAsFactors = F, sep = '\t', na.strings = "", quote = "") -> key
sort(unique(key$SMTSD)) -> tissues
# if NAs exist remove them
if (sum(is.na(tissues)) != 0) {
	tissues <- tissues[-c(which(is.na(tissues) == T))]}
if ('Cells - Leukemia cell line (CML)' %in% tissues){
	tissues <- tissues[-c(which(tissues == 'Cells - Leukemia cell line (CML)'))]}

# select only those xp that are being used
a <- lapply(xp, function(x) dim(x)[1])
a[sapply(a, is.null)] <- NA
which(unlist(a) != 0) -> xs

for (x in xs){
# fix column names
stringr::str_replace_all(colnames(xp[[x]]), '\\.', '-') -> colnames(xp[[x]])
# make list where each element is a different tissue
txp[[x]] <- list()
for(y in 1:length(tissues)){txp[[x]][[y]] <- xp[[x]][,which(colnames(xp[[x]]) %in% key[which(key$SMTSD == tissues[y]),grep('SAMPID', colnames(key))])]}
# add tissue names
names(txp[[x]]) <- tissues
# remove decimal point from the ENSG
for(y in 1:length(tissues)){
	# select only those tissues that are actually data frames
	if (class(txp[[x]][[y]]) == "data.frame"){
		rownames(txp[[x]][[y]]) <- stringr::word(xp[[x]][,1], 1, sep = stringr::fixed('.'))
		}
		}
		}

# add xp names
names(txp) <- c('xpGeneReads', 'xpTranscriptReads', 'xpGeneRPKM', 'xpTranscriptRPKM')[xs]

###

# add to mgl
# element 8 of mgl: xpGeneReads
if(xpGeneReads == TRUE){
for(x in 1:length(mgl)){
	mgl[[x]][[8]] <- list()
	for(y in 1:length(txp[[1]])){
		if (class(txp[[1]][[y]]) == "data.frame"){
			mgl[[x]][[8]][[y]] <- txp[[1]][[y]][which(rownames(txp[[1]][[y]]) %in% mgl[[x]][[1]][,2]),]
			}
			}
			}
for(x in 1:length(mgl)){names(mgl[[x]][[8]]) <- names(txp[[1]])}
unlink('xpGeneReads.txt')
}
# element 9 of mgl: xpTranscriptReads
if(xpTranscriptReads == TRUE){
for(x in 1:length(mgl)){
	mgl[[x]][[9]] <- list()
	for(y in 1:length(txp[[2]])){
		if (class(txp[[2]][[y]]) == "data.frame"){
			mgl[[x]][[9]][[y]] <- txp[[2]][[y]][which(rownames(txp[[2]][[y]]) %in% mgl[[x]][[2]][,2]),]
			}
			}
			}
for(x in 1:length(mgl)){names(mgl[[x]][[9]]) <- names(txp[[2]])}
unlink('xpTranscriptReads.txt')
}
# element 10 of mgl: xpGeneRPKM
if(xpGeneRpkm == TRUE){
for(x in 1:length(mgl)){
	mgl[[x]][[10]] <- list()
	for(y in 1:length(txp[[3]])){
		if (class(txp[[3]][[y]]) == "data.frame"){
			mgl[[x]][[10]][[y]] <- txp[[3]][[y]][which(rownames(txp[[3]][[y]]) %in% mgl[[x]][[1]][,2]),]
			}
			}
			}
for(x in 1:length(mgl)){names(mgl[[x]][[10]]) <- names(txp[[3]])}
unlink('xpGeneRpkm.txt')
}
# element 11 of mgl: xpTranscriptRPKM
if(xpTranscriptRpkm == TRUE){
for(x in 1:length(mgl)){
	mgl[[x]][[11]] <- list()
	for(y in 1:length(txp[[4]])){
		if (class(txp[[4]][[y]]) == "data.frame"){
			mgl[[x]][[11]][[y]] <- txp[[4]][[y]][which(rownames(txp[[4]][[y]]) %in% mgl[[x]][[2]][,2]),]
			}
			}
			}
for(x in 1:length(mgl)){names(mgl[[x]][[11]]) <- names(txp[[4]])}
unlink('xpTranscriptRpkm.txt')
}
}

if( download == TRUE){
	unlink('GTEx_Analysis_v6p_eQTL_covariates.tar.gz', recursive = T)
	unlink('GTEx_Analysis_v6p_eQTL_expression_matrices.tar', recursive = T)
	if(saveDownload == FALSE){
		unlink('GTEx_Data_V6_Annotations_SampleAttributesDS.txt')
		unlink('GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_reads.gct.gz')
		unlink('GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz')
		unlink('GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz')
		unlink('GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt.gz')
		unlink("./GTEx_Analysis_v6p_eQTL_expression_matrices", recursive = T)
		unlink('GTEx_Analysis_v6p_eQTL_covariates', recursive = T)
}
}

return(mgl)

}



