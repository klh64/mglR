#' Convert chr_pos to rs id 
#'
#' \code{fixSnpIds} converts chr_pos identifiers to rs ids for cisEqtl data released from GTEx
#'
#'
#' @seealso \code{\link{addCisEqtl}}
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' \dontrun{fixSnpIds(myMgl) -> myMgl}
#'
#'@export
 
fixSnpIds <- function(mgl){

# stop process if eQTL data not yet added
if(unique(unlist(lapply(mgl, function(x) class(x[[14]])))) == 'integer'){
	stop('Cis eQTL data not available. See function listElements to check which elements are present.  See function addCisEqtl to add pertinent information.') 
}
# stop process if location (element 3) has not been filled in
if (unique(unlist(lapply(mgl, function(x) class(x[[3]])))) == 'integer') stop('location data not available. see function addLoc') 

# stop process if rsIds have already been added
for (i in 1:length(mgl)) {
	for(t in 1:length(mgl[[i]][[14]])){
		if(class(mgl[[i]][[14]][[t]]) != 'integer'){ 
			if(dim(mgl[[i]][[14]][[t]])[1] >= 1){
				if(dim(mgl[[i]][[14]][[t]])[2] == 13){
					stop('rsIds have already been added')
				}}}}}

# Via biomaRt - this is dangerous because of how finicky biomaRt can be and because it downloads the whole chromosome and requires subsetting	
mart = biomaRt::useMart(biomart = 'ENSEMBL_MART_SNP', host = 'feb2014.archive.ensembl.org', path = '/biomart/martservice', dataset = 'hsapiens_snp')	
unique(do.call(rbind,listSeparate(mgl, 3))[,1])  -> chr
minimum <- list()
maximum <- list()
for(x in 1:length(chr)){
	lookHere <- which(do.call(rbind,listSeparate(mgl, 3))[,1] == chr[x])
	minimum[[x]] <- list()
	maximum[[x]] <- list()
	for (i in lookHere){
		a1 <- listSeparate(mgl, 14)[[i]]
		a1 <- a1[which(unlist(lapply(a1, class)) == 'data.frame')]
		if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) != 1 ){
			a1 <- a1[which(unlist(lapply(a1, function(x) dim(x)[1])) > 0)]
			minimum[[x]][[i]] <- min(stringr::word(do.call(rbind,a1)[,2], 2, sep = stringr::fixed('_')))
			maximum[[x]][[i]] <- max(stringr::word(do.call(rbind,a1)[,2], 2, sep = stringr::fixed('_')))
			}
	}
minimum[[x]] <- min(unlist(minimum[[x]]))	
maximum[[x]] <- max(unlist(maximum[[x]]))	
}
bm <- list()
for(x in 1:length(chr)){
bm[[x]] <- biomaRt::getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start'), filters = c('chr_name', 'chrom_start', 'chrom_end'), values = list(chr[x], minimum[[x]], maximum[[x]]), mart = mart)	
}
names(bm) <- chr
for (i in 1:length(mgl)) {
	for(t in 1:length(mgl[[i]][[14]])){
		if(class(mgl[[i]][[14]][[t]]) != 'integer'){ 
			if(dim(mgl[[i]][[14]][[t]])[1] >= 1){
			tmp <- do.call(rbind, stringr::str_split(mgl[[i]][[14]][[t]][,2], ('_')))[,1:2] 
			y <- which(names(bm) == tmp[1,1])
			mgl[[i]][[14]][[t]] <- cbind(mgl[[i]][[14]][[t]], bm[[y]][match(tmp[,2], bm[[y]][,3]),1])
			colnames(mgl[[i]][[14]][[t]])[13] <- 'rsId'
			mgl[[i]][[14]][[t]][,13] <- as.character(mgl[[i]][[14]][[t]][,13])
			mgl[[i]][[14]][[t]][which(is.na(mgl[[i]][[14]][[t]][,13]) == T),13] <- mgl[[i]][[14]][[t]][which(is.na(mgl[[i]][[14]][[t]][,13]) == T),2]	
}}}}
return(mgl)
}