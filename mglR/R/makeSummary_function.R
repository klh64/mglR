#' Returns a summary table of the number of eQTL and GWAS SNPs for each gene
#'
#' \code{makeSummary} returns a dataframe with the number of rows corresponding to the number of genes in the list and columns: gene, Number of eQTL SNPs, Number of GWAS Catalog SNPs, Number of GRASP SNPs. 
#'
#' Provides a brief summary of variants that have been associated with expression and clinical traits. Of interest maybe genes that have evidence of regulatory variants but have not yet been tied to a clinical phenotype or those genes that have evidence of a clinical phenotype but have not yet been shown to have regulatory variants
#'
#' @family output
#'
#' @param mgl Object of class 'mgl'; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param saveFile Logical flag indicating whether a csv ('Summary.csv') should be saved in the current directory
#'
#' @examples
#' makeSummary(myMgl, saveFile = TRUE)
#'
#'@export
#'@importFrom utils write.csv
 
makeSummary <- function(mgl, saveFile = FALSE){

# stop process if none of the SNP based elements have been filled in
tmp <- c(unique(unlist(lapply(mgl, function(x) class(x[[12]])))), unique(unlist(lapply(mgl, function(x) class(x[[13]])))), unique(unlist(lapply(mgl, function(x) class(x[[14]])))), unique(unlist(lapply(mgl, function(x) class(x[[15]])))), unique(unlist(lapply(mgl, function(x) class(x[[16]])))), unique(unlist(lapply(mgl, function(x) class(x[[17]])))), unique(unlist(lapply(mgl, function(x) class(x[[18]])))), unique(unlist(lapply(mgl, function(x) class(x[[19]])))))

if(length(unique(tmp)) == 1){
	if(unique(tmp) == 'integer')
	stop('SNP data not available. See function listElements to check which elements are present.  See functions addDnase, addTransEqtl, addCisEqtl, addGrasp, addgwasCatalog, addPtv, addSqtlAltrans, addSqtlSeek to add pertinent information.') 
}


### get only those elements that are filled in
dane <-as.list(1:length(mgl))
for(i in 1:length(mgl)){
	dane[[i]] <- as.list(1:8)
	names(dane[[i]]) <- c('DNAse', 'TranseQTL', 'CiseQTL', 'sqtlSeek', 'sqtlAltrans', 'Ptv', 'GwasCatalog', 'GRASP')
}
names(dane) <- names(mgl)

for(i in 1:length(mgl)){
if(tmp[2] != "integer"){
	dane[[i]][[2]] <- as.character(unique(do.call(rbind,mgl[[i]][[13]])[,2]))}
if(tmp[3] != "integer"){
	a1 <- listSeparate(mgl, 14)[[i]]
	a1 <- a1[which(unlist(lapply(a1, class)) == 'data.frame')]
	if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) != 1 ){
		a1 <- a1[which(unlist(lapply(a1, function(x) dim(x)[1])) > 0)]
		a1 <- unique(do.call(rbind,a1))
		if(dim(a1)[2] == 12){
			dane[[i]][[3]] <- a1[,2]
			message("Warning: It maybe preferable to convert SNP Ids to rs identifiers where possible.  See fixSnpIds.")}
		if(dim(a1)[2] == 13){
			dane[[i]][[3]] <- a1[,13]}}
	if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) == 1 ) {	
	dane[[i]][[3]] <- 'None'}}
if(tmp[4] != "integer"){
	dane[[i]][[4]] <- unique(do.call(rbind,mgl[[i]][[15]])[,1])}		
if(tmp[5] != "integer"){
	dane[[i]][[5]] <- unique(do.call(rbind,mgl[[i]][[16]])[,1])}
if(tmp[6] != "integer"){
	dane[[i]][[6]] <- unique(do.call(rbind,mgl[[i]][[17]])[,1])}
if(tmp[7] != "integer"){
	dane[[i]][[7]] <- unique(mgl[[i]][[18]][,22])}
if(tmp[8] != "integer"){
	dane[[i]][[8]] <- unique(mgl[[i]][[19]][,7])}
if(tmp[1] != "integer"){
	dane[[i]][[1]] <- unique(mgl[[i]][[12]][which(mgl[[i]][[12]][,13] == 'imbalanced_(5%_FDR)' | mgl[[i]][[12]][,13] == 'imbalanced_(0.1%_FDR)'),4])}
	}

# Number of SNPs
nums <-list()
for (i in 1:length(dane)){
	nums[[i]] <- lapply(dane[[i]], length)
}
names(nums) <- names(mgl)


# fix length for places with 'None' as CiseQTL
if(length(unlist(lapply(dane, function(x) grep('None',x)))) != 0){
	pos <- unlist(lapply(dane, function(x) grep('None',x)))
	for(z in 1:length(pos)){
		nums[[pos[z]]][grep('None', dane[[pos[z]]])] <- 0
	}
}

# Making a table
do.call(rbind, nums) -> summary
rownames(summary) <- names(mgl)

# Add message about data missing from mgl for those elements that are not filled in	
for(x in 1:8){
	if(tmp[x] == 'integer'){
		for(i in 1:length(mgl)){
			summary[i,x] <- NA
		}
				message (paste(names(dane[[1]])[x], 'missing from mgl list'))
	}
}

# Saving if flag
if (saveFile == TRUE){
write.csv(summary, 'summary.csv')}

#message('\nPlease cite datasources as appropriate. \n \n Element 12 - DNAse: \n Maurano, Matthew T, Eric Haugen, Richard Sandstrom, Jeff Vierstra, Anthony Shafer, Rajinder Kaul, and John A Stamatoyannopoulos. 2015. “Large-Scale Identification of Sequence Variants Influencing Human Transcription Factor Occupancy in Vivo.” Nature Genetics 47 (12): 1393–1401. doi:10.1038/ng.3432.\n \n Elements 13-17 GTEx: \n The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n \n Element 15 - sqtlSeek: \n (1)Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 5:4698 doi: 10.1038/ncomms5698 (2014) \n (2) "Multi-tissue analysis of gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014. \n \n  Element 16 - sqtlAltrans: \n (1) "Multi-tissue analysisof gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014.\n (2) Ongen, H., & Dermitzakis, E. T. (2015). Alternative Splicing QTLs in European and African Populations. American Journal of Human Genetics, 97(4), 567–575. http://doi.org/10.1016/j.ajhg.2015.09.004 \n \n Element 17 - pqtl: \n Rivas MA, Pirinen M, Conrad DF, et al. Impact of predicted protein-truncating genetic variants on the human transcriptome. Science (New York, NY). 2015;348(6235):666-669. doi:10.1126/science.1261877. \n \n Element 18 - NHGRI-EBI GWAS Catalog: \n Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H.The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006 \n \n Element 19 - GRASP: \n (1) Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390 genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP Build 2.0.0.0\n (2) Carey V (2016). grasp2db: grasp2db, sqlite wrap of GRASP 2.0. R package version 0.1.14.\n')

return(summary)

}

