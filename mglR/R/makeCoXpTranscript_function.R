#' Tests pairwise co-expression between user defined transcripts in specified tissues
#'
#' \code{makeCoXpTranscript} returns a list of cor.test output for each transcript pair in a given tissue
#'
#' Uses the cor.test function in R to test for co-expression between transcripts.  Option to generate scatter plots and a summary table of the correlations to be saved as a pdf and csv, respectively.  Any number of transcripts and tissues can be included. Counts (data = 9) or RPKM (data = 11) can be used. 
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param genes Character vector indicating genes corresponding to transcripts to be considered (note gene names must match names(mgl))
#'
#' @param transcripts Character vector indicating transcripts to be considered (see mgl[[x]][[2]] where x is number 
#'    corresponding to gene of interest)
#'
#' @param data Number corresponding to element of mgl containing data of interest: 8 for count data and 10 for RPKM
#'
#' @param makePlot A logical flag indicating whether a pdf file ('CoXpGene.pdf') should be saved in the current directory
#'
#' @param saveFile A logical flag indicating whether a csv file ('CoXpGene.csv') should be saved in the current directory. 'CoXpGene.csv' contains cor.test results
#'
#' @param tissue Character vector indicating tissues to be considered (note tissue names must match those provided here: c('All', 'Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood'))
#'
#' @examples
#' makeCoXpTranscript(myMgl, transcripts = 
#'     c('ENST00000344348', 'ENST00000414907', 'ENST00000478719'), 
#'     genes = c('NCOA4', 'NCOA4', 'MSMB'), tissue = c('Prostate'), 
#'     data = 9, makePlot = FALSE, saveFile = FALSE)
#'
#'@export
#'@importFrom stats cor.test
#'@importFrom grDevices pdf dev.off
#'@importFrom graphics layout par plot plot.new legend
#'@importFrom utils combn write.table
 
makeCoXpTranscript <- function(mgl, transcripts = c(), genes = c(), data = c(9, 11), makePlot = FALSE, saveFile = FALSE, tissue = c('All', 'Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood')){

if ('All' %in% tissue){
	tissue <- c('Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland', 'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala', 'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)', 'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex', 'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra', 'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes', 'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction', 'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube', 'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood')
}

# table with all possible pairwise combinations of genes
combn(1:length(transcripts), 2) -> ns

res <- list()
# STEP 1: Correlation test
for(x in 1:dim(ns)[2]){
	res[[x]] <- list()
	for(t in 1:length(tissue)){
# define which pair
n1 <- ns[1,x]
n2 <- ns[2,x]

# pull data as numeric
g1 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][which(rownames(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]]) == transcripts[n1]),]))
g2 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][which(rownames(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]]) == transcripts[n2]),]))

# run correlation test 
res[[x]][[t]] <- suppressWarnings(cor.test(g1, g2))
}}

if (makePlot == TRUE){
# STEP 2: Plot
pdf('CoXpTranscript.pdf')	
for(x in 1:dim(ns)[2]){
	for(t in 1:length(tissue)){
# define which pair
n1 <- ns[1,x]
n2 <- ns[2,x]
# pull data as numeric
g1 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][which(rownames(mgl[[which(names(mgl) == genes[n1])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]]) == transcripts[n1]),]))
g2 <- as.numeric(as.character(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]][which(rownames(mgl[[which(names(mgl) == genes[n2])]][[data]][[which(names(mgl[[1]][[data]]) == tissue[t])]]) == transcripts[n2]),]))

# alter layout so legend with correlation statistics will fit below
layout(rbind(1,2), heights = c(7,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(g1, g2, main = paste('CoExpression of ', transcripts[n1], '(', genes[n1], ')', ' and ', transcripts[n2], '(', genes[n2], ')', ' in ', tissue[t], sep = ""), cex.main = 0.9, xlab = paste(genes[n1], " (", names(mgl[[which(names(mgl) == genes[n1])]])[data], ")", sep = ""), ylab = paste(genes[n2], " (", names(mgl[[which(names(mgl) == genes[n1])]])[data], ")", sep = ""))

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
cbind(rep(paste(paste(transcripts[ns[1,]], '(', genes[ns[1,]], ')', sep = ""), paste(transcripts[ns[2,]], '(', genes[ns[2,]], ')', sep = ""), sep = "|"), each = length(tissue)), tmp) -> tmp
colnames(tmp)[1] <- 'Transcript Pair'
write.table(tmp, 'CoXpTranscript.csv', sep = ',', row.names = FALSE)
}

# adding names
names(res) <- paste(paste(transcripts[ns[1,]], '(', genes[ns[1,]], ')', sep = ""), paste(transcripts[ns[2,]], '(', genes[ns[2,]], ')', sep = ""), sep = "|")

for(x in 1:length(res)){
	names(res[[x]]) <- tissue
}


#message('Please acknowledge: The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n')

return(res)


}



