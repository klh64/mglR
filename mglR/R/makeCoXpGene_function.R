#' Tests pairwise co-expression for normalized expression data
#'
#' \code{makeCoXpGene} returns a list of cor.test output for each user-defined gene pair in a given tissue
#'
#' Uses the cor.test function in R to test for co-expression between genes.  Option to generate scatter plots and a summary table of the correlations to be saved as a pdf and csv, respectively.  Any number of genes and tissues can be included.  Normalized (data = 7), counts (data = 8), or RPKM (data = 10) can be used. 
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param genes Character vector indicating genes to be considered (note gene names must match names(mgl))
#'
#' @param data Number corresponding to element of mgl containing data of interest: 8 for count data and 10 for RPKM
#'
#' @param makePlot Logical flag - TRUE indicates pdf file named 'CoXpGene.pdf'  
#'    will be saved in the current directory
#'
#' @param saveFile Logical flag - TRUE indicates csv file named 'CoXpGene.csv'  
#'    will be saved in the current directory. 'CoXpGene.csv' contains cor.test results
#'
#' @param tissue Character vector indicating tissues to be considered (note tissue names must match those provided here: c('All', Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood'))
#'

#' @examples
#' makeCoXpGene(myMgl, gene = c('RP11-109G10.2', 'NCOA4'), makePlot = FALSE, 
#'     saveFile = FALSE, data = 7, tissue = c('Prostate')) -> results
#'
#'@export
#'@importFrom stats cor.test
#'@importFrom grDevices pdf dev.off
#'@importFrom graphics layout par plot plot.new legend
#'@importFrom utils combn write.table
 
makeCoXpGene <- function(mgl, genes = c(), makePlot = FALSE, saveFile = FALSE, data = c(7, 8, 10), tissue = c('All', 'Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood')){

if (data %in% c(8, 10)){
	if ('All' %in% tissue){
		tissue <- c('Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood')
}
}

if (data == 7){
	if ('All' %in% tissue){
		tissue <- c('Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Breast_Mammary_Tissue', 'Cells_EBV-transformed_lymphocytes', 'Cells_Transformed_fibroblasts', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Liver', 'Lung', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood')
	}
	if('Adipose - Subcutaneous' %in% tissue){
		tissue[which(tissue == 'Adipose - Subcutaneous')] <- 'Adipose_Subcutaneous'
	}
	if('Adipose - Visceral (Omentum)' %in% tissue){
		tissue[which(tissue == 'Adipose - Visceral (Omentum)')] <- 'Adipose_Visceral_Omentum'
	}
	if('Adrenal Gland' %in% tissue){
		tissue[which(tissue == 'Adrenal Gland')] <- 'Adrenal_Gland'
	}
	if('Artery - Aorta' %in% tissue){
		tissue[which(tissue == 'Artery - Aorta')] <- 'Artery_Aorta'
	}
	if('Artery - Coronary' %in% tissue){
		tissue[which(tissue == 'Artery - Coronary')] <- 'Artery_Coronary'
	}
	if('Artery - Tibial' %in% tissue){
		tissue[which(tissue == 'Artery - Tibial')] <- 'Artery_Tibial'
	}
	if('Bladder' %in% tissue){
		stop('Bladder not available for normalized expression')
	}
	if ('Brain - Amygdala' %in% tissue){
		stop('Brain - Amygdala not available for normalized expression')
	}
	if('Brain - Anterior cingulate cortex (BA24)' %in% tissue){
		tissue[which(tissue == 'Brain - Anterior cingulate cortex (BA24)')] <- 'Brain_Anterior_cingulate_cortex_BA24'
	}
	if('Brain - Caudate (basal ganglia)' %in% tissue){
		tissue[which(tissue == 'Brain - Caudate (basal ganglia)')] <- 'Brain_Caudate_basal_ganglia'
	}
	if('Brain - Cerebellar Hemisphere' %in% tissue){
		tissue[which(tissue == 'Brain - Cerebellar Hemisphere')] <- 'Brain_Cerebellar_Hemisphere'
	}
	if('Brain - Cerebellum' %in% tissue){
		tissue[which(tissue == 'Brain - Cerebellum')] <- 'Brain_Cerebellum'
	}
	if('Brain - Cortex' %in% tissue){
		tissue[which(tissue == 'Brain - Cortex')] <- 'Brain_Cortex'
	}
	if('Brain - Frontal Cortex (BA9)' %in% tissue){
		tissue[which(tissue == 'Brain - Frontal Cortex (BA9)')] <- 'Brain_Frontal_Cortex_BA9'
	}
	if('Brain - Hippocampus' %in% tissue){
		tissue[which(tissue == 'Brain - Hippocampus')] <- 'Brain_Hippocampus'
	}
	if('Brain - Hypothalamus' %in% tissue){
		tissue[which(tissue == 'Brain - Hypothalamus')] <- 'Brain_Hypothalamus'
	}
	if('Brain - Nucleus accumbens (basal ganglia)' %in% tissue){
		tissue[which(tissue == 'Brain - Nucleus accumbens (basal ganglia)')] <- 'Brain_Nucleus_accumbens_basal_ganglia'
	}
	if('Brain - Putamen (basal ganglia)' %in% tissue){
		tissue[which(tissue == 'Brain - Putamen (basal ganglia)')] <- 'Brain_Putamen_basal_ganglia'
	}
	if('Brain - Spinal cord (cervical c-1)' %in% tissue){
		stop('Brain - Spinal cord (cervical c-1) not available for normalized expression')
	}
	if('Brain - Substantia nigra' %in% tissue){
		stop('Brain - Substantia nigra not available for normalized expression')
	}
	if('Breast - Mammary Tissue' %in% tissue){
		tissue[which(tissue == 'Breast - Mammary Tissue')] <- 'Breast_Mammary_Tissue'
	}
	if('Cells - EBV-transformed lymphocytes' %in% tissue){
		tissue[which(tissue == 'Cells - EBV-transformed lymphocytes')] <- 'Cells_EBV-transformed_lymphocytes'
	}
	if('Cells - Transformed fibroblasts' %in% tissue){
		tissue[which(tissue == 'Cells - Transformed fibroblasts')] <- 'Cells_Transformed_fibroblasts'
	}
	if('Cervix - Ectocervix' %in% tissue){
		stop('Cervix - Ectocervix not available for normalized expression')
	}
	if('Cervix - Endocervix' %in% tissue){
		stop ('Cervix - Endocervix not available for normalized expression')
	}
	if('Colon - Sigmoid' %in% tissue){
		tissue[which(tissue == 'Colon - Sigmoid')] <- 'Colon_Sigmoid'
	}
	if('Colon - Transverse' %in% tissue){
		tissue[which(tissue == 'Colon - Transverse')] <- 'Colon_Transverse'
	}
	if('Esophagus - Gastroesophageal Junction' %in% tissue){
		tissue[which(tissue == 'Esophagus - Gastroesophageal Junction')] <- 'Esophagus_Gastroesophageal_Junction'
	}
	if('Esophagus - Mucosa' %in% tissue){
		tissue[which(tissue == 'Esophagus - Mucosa')] <- 'Esophagus_Mucosa'
	}
	if('Esophagus - Muscularis' %in% tissue){
		tissue[which(tissue == 'Esophagus - Muscularis')] <- 'Esophagus_Muscularis'
	}
	if('Fallopian Tube' %in% tissue){
		stop('Fallopian Tube not available for normalized expression')
	}
	if('Heart - Atrial Appendage' %in% tissue){
		tissue[which(tissue == 'AHeart - Atrial Appendage')] <- 'Heart_Atrial_Appendage'
	}
	if('Heart - Left Ventricle' %in% tissue){
		tissue[which(tissue == 'Heart - Left Ventricle')] <- 'Heart_Left_Ventricle'
	}
	if('Kidney - Cortex' %in% tissue){
		stop('Kidney - Cortex not available for normalized expression')
	}
	if('Minor Salivary Gland' %in% tissue){
		stop("Minor Salivary Gland not available for normalized expression")
	}
	if('Muscle - Skeletal' %in% tissue){
		tissue[which(tissue == 'Muscle - Skeletal')] <- 'Muscle_Skeletal'
	}
	if('Nerve - Tibial' %in% tissue){
		tissue[which(tissue == 'Nerve - Tibial')] <- 'Nerve_Tibial'
	}
	if('Skin - Not Sun Exposed (Suprapubic)' %in% tissue){
		tissue[which(tissue == 'Skin - Not Sun Exposed (Suprapubic)')] <- 'Skin_Not_Sun_Exposed_Suprapubic'
	}
	if('Skin - Sun Exposed (Lower leg)' %in% tissue){
		tissue[which(tissue == 'Skin - Sun Exposed (Lower leg)')] <- 'Skin_Sun_Exposed_Lower_leg'
	}
	if('Small Intestine - Terminal Ileum' %in% tissue){
		tissue[which(tissue == 'Small Intestine - Terminal Ileum')] <- 'Small_Intestine_Terminal_Ileum'
	}
	if('Whole Blood' %in% tissue){
		tissue[which(tissue == 'Whole Blood')] <- 'Whole_Blood'
	}
}

# table with all possible pairwise combinations of genes
combn(1:length(genes), 2) -> ns

res <- list()
# STEP 1: Correlation test
for(x in 1:dim(ns)[2]){
	res[[x]] <- list()
	for(t in 1:length(tissue)){
# define which pair
n1 <- ns[1,x]
n2 <- ns[2,x]

# print message if normalized expression is missing
if(data == 7){
	if(stringr::str_sub(rownames(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][1,]), 1,4) !=  'ENSG'){
		stop(paste('Normalized expression not available for', genes[n1], 'in', tissue[t]))
	}
	if(stringr::str_sub(rownames(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][1,]), 1,4) !=  'ENSG'){
		stop(paste('Normalized expression not available for', genes[n2], 'in', tissue[t]))
	}	
}

# pull data as numeric
g1 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][1,]))
g2 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][1,]))

# run correlation test 
res[[x]][[t]] <-cor.test(g1, g2)
}}

if (makePlot == TRUE){
# STEP 2: Plot
pdf('CoXpGene.pdf')	
for(x in 1:dim(ns)[2]){
	for(t in 1:length(tissue)){
# define which pair
n1 <- ns[1,x]
n2 <- ns[2,x]
# pull data as numeric
g1 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][1,]))
g2 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][1,]))

# alter layout so legend with correlation statistics will fit below
layout(rbind(1,2), heights = c(7,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(g1, g2, main = paste('CoExpression of ', genes[n1],' and ', genes[n2], ' in ', tissue[t], sep = ""), xlab = paste(genes[n1], " (", names(mgl[[which(names(mgl) == genes[n1])]])[data], ")", sep = ""), ylab = paste(genes[n2], " (", names(mgl[[which(names(mgl) == genes[n1])]])[data], ")", sep = ""))

# add correlation statistics to plot
par(mar = c(0,0,0,0))
plot.new()
legend('center', legend = paste(res[[x]][[t]]$method, ' estimate: ', round(res[[x]][[t]]$estimate, digits = 2), ' , p-value: ', round(res[[x]][[t]]$p.value, digits = 4), sep = ""))
		}
	}
dev.off()
}


# STEP 3: saving
if (saveFile == TRUE){

tmp <- list()
for(x in 1:dim(ns)[2]){
	tmp[[x]] <- list()
	for(t in 1:length(tissue)){
		tmp[[x]][[t]] <- unlist(res[[x]][[t]])
	}
	tmp[[x]] <- do.call(rbind, tmp[[x]])
	tmp[[x]][,8] <- tissue
}
do.call(rbind, tmp) -> tmp
cbind(rep(paste(genes[ns[1,]], genes[ns[2,]], sep = "|"), each = length(tissue)), tmp) -> tmp
colnames(tmp)[1] <- 'Gene Pair'
write.table(tmp, 'CoXpGene.csv', sep = ',', row.names = FALSE)
}

# adding names
names(res) <- paste(genes[ns[1,]], genes[ns[2,]], sep = "|")

for(x in 1:length(res)){
	names(res[[x]]) <- tissue
}

#message('Please acknowledge: The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n')


return(res)

}



