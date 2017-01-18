#' Returns SNPs identified as cis-eQTLs for more than one gene
#'
#' \code{makeMultiEqtl} returns information about SNPs that appear as eQTLs for more than one of the candidate genes.  It can be done across all tissues by selecting tissue = 'All' or for a specific tissue.  Note: multiple tissues can be selected at once.
#'
#' The output is a four element list - with each of the four elements being a list corresponding to the number of tissues selected in the function.  For the first element (\emph{'fullTable'}) the list is comprised of dataframes with three columns: gene, tissue, and snp.  The second element of the output (\emph{'ranking'}) is a list of tables ranking the number of genes each eQTL is assigned to.  The third element (\emph{'topSnps'}) of the output is a list of vectors of eQTLs appearing for the maximum number of genes (note the exact number will vary based on candidate genes and tissues considered).  The fourth element (\emph{'topSummary'}) is a list of the fullTable objects subsetted to include only the top SNPs - this is saved in the csv file if \emph{saveFile = TRUE}.  The plot that is generated (\emph{'MultiEqtl.pdf'}) is a heatmap of eQTLs for the selected tissue displaying SNPs along the x-axis and genes along the y-axis.  SNPs that are eQTLs for a given gene are colored in blue, while those that are eQTLs for another gene are colored in black.  Genes are clustered by the number of shared eQTLs as shown in the dendogram on the left.
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param makePlot A logical flag indicating whether a pdf ('MultiEqtl.pdf') should be saved in the current directory
#'
#' @param saveFile A logical flag indicating whether a csv ('MulitEqtl.csv') should be saved in the current directory. 'MultiEqtl.csv' contains those eQTLs, along with their corresponding genes and tissues, that serve as eQTLs for the greatest number of genes
#'
#' @param cis A logical flag indicating whether cisEqtls should be included
#'
#' @param trans A logical flag indicating whether transEqtls should be included
#'
#' @param tissue Character vector indicating tissues to be considered (note tissue names much match those provided here: c('All','Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Breast_Mammary_Tissue', 'Cells_EBV-transformed_lymphocytes', 'Cells_Transformed_fibroblasts', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Liver', 'Lung', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood'))
#'
#' @examples
#' makeMultiEqtl(myMgl, makePlot = FALSE, saveFile = FALSE, cis = TRUE, trans = TRUE, 
#'     tissue = c('All', 'Thyroid', 'Spleen')) -> results
#'
#'@export
#'@importFrom grDevices pdf dev.off
#'@importFrom graphics par mtext
#'@importFrom utils write.table 
 
makeMultiEqtl <- function(mgl, makePlot = FALSE, saveFile = FALSE, cis = TRUE, trans = TRUE, tissue = c('All','Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Breast_Mammary_Tissue', 'Cells_EBV-transformed_lymphocytes', 'Cells_Transformed_fibroblasts', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Liver', 'Lung', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood')){

# print error message if eQTL data has not yet been filled in 
if(unique(unlist(lapply(mgl, function(x) class(x[[14]])))) == 'integer' & unique(unlist(lapply(mgl, function(x) class(x[[13]])))) == 'integer')
stop('eQTL data not available. see functions addTransEqtl addCisEqtl') 


# STEP 1: pull data
# Pull eQTL data from mgl (14th elemeent)
cSNPs <- list()
tSNPs <- list()
eSNPs <- list()
for (i in 1:length(mgl)) {
# cisEqtls
if(cis == TRUE){
	a1 <- listSeparate(mgl, 14)[[i]]
	a1 <- a1[which(unlist(lapply(a1, class)) == 'data.frame')]
	if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) != 1 ){
		a1 <- a1[which(unlist(lapply(a1, function(x) dim(x)[1])) > 0)]
		a1 <- unique(do.call(rbind,a1))
		if(dim(a1)[2] == 12){
			cSNPs[[i]] <- cbind(names(mgl)[i], stringr::word(rownames(a1), 1, sep = stringr::fixed('.')), a1[,2], 'cis')
			message("Warning: It maybe preferable to convert SNP Ids to rs identifiers where possible.  See fixSnpIds.")}
		if(dim(a1)[2] == 13){
			cSNPs[[i]] <- cbind(names(mgl)[i], stringr::word(rownames(a1), 1, sep = stringr::fixed('.')), a1[,13], 'cis')}
}
	if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) == 1){
		cSNPs[[i]] <- a1[[1]]
	}
}
# transEqtls
if(trans == TRUE){
	tSNPs[[i]] <- do.call(rbind, mgl[[i]][[13]])
	# 
	if(dim(tSNPs[[i]])[1] != 0){
	# reformat to add the gene name and tissue
	tSNPs[[i]] <- cbind(names(mgl)[i],stringr::word(rownames(tSNPs[[i]]), 1, sep = stringr::fixed('.')), as.character(tSNPs[[i]][,2]), 'trans')}
}
# combine
	colnames(cSNPs[[i]]) <- NULL
	colnames(tSNPs[[i]]) <- NULL
	eSNPs[[i]] <- rbind(cSNPs[[i]], tSNPs[[i]])
}

# Remove those genes that do not have any eQTLs
tmp <- which(unlist(lapply(eSNPs, function(x) dim(x)[2])) != 4)
if (length(tmp) != 0){
	eSNPs <- eSNPs[-c(tmp)]
}

# put into one dataframe
do.call(rbind, eSNPs) -> tmp
# remove rows where the SNP has . instead of an id
if (length(grep('\\.', tmp[,3])) != 0){
tmp <- tmp[-c(grep('\\.', tmp[,3])),]
}
# remove NAs
if(sum(is.na(tmp[,3])) != 0){
tmp <- tmp[-c(which(is.na(tmp[,3]) == T)),]	
}

# put into multiple dataframes
# with gene name and snp (drops tissue)
unique(tmp[,c(1,3)]) -> etab2
# with gene name, tissue, and snp
unique(tmp) -> etab3

# STEP 2: snps with maximum number of genes & summary table for top
# count is a vector with all snps and the number of genes it is an eQTL before - this is independent of tissue and is built from the etab2 object
count <- list()
# top is a subset of count for those snps that have the maximum number of snps
top <- list()
# summs is a subset of etab3 for those snps in the top object 
summs <- list()
# filtered is a vector of those snp ids that appear for more than one gene - independent of tissue; this is used to reduce the dimension in plotting heatmaps
filtered <- list()
for(x in 1:length(tissue)){
## use etab2 if independent of tissue type (i.e. if tissue = 'All')	
	if(tissue[x] == 'All'){
		table(etab2[,2]) -> count[[x]]
		count[[x]] <- sort(count[[x]], decreasing = T)
		top[[x]] <- count[[x]][which(count[[x]] == max(count[[x]]))]
		etab3[which(etab3[,3] %in% names(top[[x]])),] -> summs[[x]]
		summs[[x]] <- summs[[x]][order(summs[[x]][,3]),]
		filtered[[x]] <- names(count[[x]])[which(count[[x]] >= 2)]
	}
## use etab3 if there is a specific tissue of interest
	else{
		count[[x]] <- table(etab3[which(etab3[,2] == tissue[x]) ,3])
		count[[x]] <- sort(count[[x]], decreasing = T)
		top[[x]] <- count[[x]][which(count[[x]] == max(count[[x]]))]
		etab3[which(etab3[,2] == tissue[x] & etab3[,3] %in% names(top[[x]])),] -> summs[[x]]
		summs[[x]] <- summs[[x]][order(summs[[x]][,3]),]
		filtered[[x]] <- names(count[[x]])[which(count[[x]] >= 2)]
	}
}

# STEP 3: plotting
if (makePlot == TRUE){
	
pdf('MultiEqtl.pdf')
for(x in 1:length(tissue)){

## if independent of tissue (i.e. tissue = 'All')
	if(tissue[x] == 'All'){
# restrict to only those SNPs that appear more than once (independent of tissue)
etab2.f <- etab2[which(etab2[,2] %in% filtered[[x]]),]

# stop if there are no eQTLs
if(dim(etab2.f)[1] == 0) stop('No multi eQTLs') 

# make matrix
emat <- as.matrix(table(etab2.f[,1], etab2.f[,2]))

# plot
gplots::heatmap.2(emat, dendrogram = 'row', density.info = 'none', trace = 'none', cexRow = 0.8, cexCol = 0.8,  labCol = 'SNPs', srtCol = 0, col = c('black', 'blue'), key.xlab = 'Eqtl Status', key.xtickfun = function(){side <- 1; line <- 0; col <- par("col.axis"); font <- par("font.axis"); mtext("Not an Eqtl", side=side, at=0, adj=0, line=line, cex=0.6, col=col, font=font); mtext("Eqtl for gene", side=side, at=1, adj=1, line=line, cex=0.6, col=col, font=font); return(list(labels=FALSE, tick=FALSE))}, main = paste('Shared Eqtls in ', tissue[x], ' Tissue', sep = ""))
	}

## if there is a specific tissue of interest
	else{
# restrict to only those SNPs that appear more than once (independent of tissue)
		etab3.f <- etab3[which(etab3[,3] %in% filtered[[x]]),]

# stop if there are no eqtls for this tissue
if(dim(etab3.f)[1] == 0) stop(paste('No multi eQTLs for', tissue[x]))

# make matrix and restrict to those eQTLs that apply to the tissue of interest
emat <- as.matrix(table(etab3.f[which(etab3.f[,2] == tissue[x]),1], etab3.f[which(etab3.f[,2] == tissue[x]),3]))

# make sure the color is blue if all the snps are shared
if( unique(emat) == 1){
gplots::heatmap.2(emat, dendrogram = 'row', density.info = 'none', trace = 'none', cexRow = 0.8, cexCol = 0.8,  labCol = 'SNPs', srtCol = 0, col = c('blue', 'black'), key.xlab = 'Eqtl Status', key.xtickfun = function(){side <- 1; line <- 0; col <- par("col.axis"); font <- par("font.axis"); mtext("Not an Eqtl", side=side, at=0, adj=0, line=line, cex=0.6, col=col, font=font); mtext("Eqtl for gene", side=side, at=1, adj=1, line=line, cex=0.6, col=col, font=font); return(list(labels=FALSE, tick=FALSE))}, main = paste('Shared Eqtls in ', tissue[x], ' Tissue', sep = ""))	
}

# if there are both 0's and 1's 
if(unique(emat) != 1){
gplots::heatmap.2(emat, dendrogram = 'row', density.info = 'none', trace = 'none', cexRow = 0.8, cexCol = 0.8,  labCol = 'SNPs', srtCol = 0, col = c('black', 'blue'), key.xlab = 'Eqtl Status', key.xtickfun = function(){side <- 1; line <- 0; col <- par("col.axis"); font <- par("font.axis"); mtext("Not an Eqtl", side=side, at=0, adj=0, line=line, cex=0.6, col=col, font=font); mtext("Eqtl for gene", side=side, at=1, adj=1, line=line, cex=0.6, col=col, font=font); return(list(labels=FALSE, tick=FALSE))}, main = paste('Shared Eqtls in ', tissue[x], ' Tissue', sep = ""))
}	
	}
}
dev.off()
}

# STEP 4: saving
if (saveFile == TRUE){
# empty header to establish csv
write.table("", "MultiEqtl.csv" , sep = ',', col.names = FALSE, row.names = FALSE)

# loop to fill in data for each tissue
for(x in 1:length(tissue)){
# paste tissue name
	write.table(tissue[x], "MultiEqtl.csv", sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
# paste gene, tissue, snp id for the top snps 
	write.table(summs[[x]], "MultiEqtl.csv", sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
# empty line
	write.table("", "MultiEqtl.csv" , sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
}
}

# compile results into list
res <- list(fullTable = etab3, ranking = count, topSnps = top, topSummary = summs)

#message('Please acknowledge: The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.')

return(res)

}


