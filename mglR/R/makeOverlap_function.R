#' Returns SNPs that appear in two user-defined groups for a given gene
#'
#' \code{makeSummary} returns a list with the first element (named overlap) being a list corresponding to the number of genes in mgl with each element as a vector of SNPs that overlap between the two user-defined groups.  The second two elements (named snpsA and snpsB)are the user-defined groups 
#'
#' Provides a character vector of SNPs for each gene that appear in two different groups as defined by the user.  These groups can be built from any element in mgl containing SNP information: eqtl, sqtlSeek, sqtlAltrans, pqtl, dnase, gwasCatalog, grasp.  Can be used to determine those SNPs that appear as both regulatory variants and trait associated variants, for example.
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param snpsA A character vector indicating which elements of mgl should be used to define group A - any combination of 'eqtl', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'dnase', 'gwasCatalog', 'grasp'
#'
#' @param snpsB A character vector indicating which elements of mgl should be used to define group B - any combination of 'eqtl', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'dnase', 'gwasCatalog', 'grasp'
#'
#' @param saveFile A logical flag indicating whether a csv file ('Overlap[groups].csv') should be saved in the current directory
#'
#' @examples
#' makeOverlap(myMgl, snpsA = c('cisEqtls'), snpsB = c('gwasCatalog'), saveFile = TRUE) -> myOverlap
#'
#'@export
#'@importFrom utils write.csv
 
makeOverlap <- function(mgl, snpsA = c('transEqtls','cisEqtls', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'dnase', 'gwasCatalog', 'grasp'), snpsB = c('transEqtls','cisEqtls', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'dnase', 'gwasCatalog', 'grasp'), saveFile = FALSE){

# Empty list for pulling data
elements <- as.list(1:8)
names(elements) <- c('transEqtls', 'cisEqtls', 'sqtlSeek', 'sqtlAltrans', 'pqtl', 'gwasCatalog', 'grasp', 'dnase')

# Numbers of list elements for each group
which(names(elements) %in% snpsA) -> posA
which(names(elements) %in% snpsB) -> posB

# stop process if the necessary elements have not been filled in

tmp <- c(unique(unlist(lapply(mgl, function(x) class(x[[13]])))), unique(unlist(lapply(mgl, function(x) class(x[[14]])))), unique(unlist(lapply(mgl, function(x) class(x[[15]])))), unique(unlist(lapply(mgl, function(x) class(x[[16]])))), unique(unlist(lapply(mgl, function(x) class(x[[17]])))), unique(unlist(lapply(mgl, function(x) class(x[[18]])))), unique(unlist(lapply(mgl, function(x) class(x[[19]])))), unique(unlist(lapply(mgl, function(x) class(x[[12]])))))

lapply(posA, function(x) if(tmp[x] == 'integer') stop(paste(names(elements)[x], ' data not available. see functions addDnase, addTransEqtl, addCisEqtl, addGrasp, addgwasCatalog, addPtv, addSqtlAltrans, addSqtlSeek', sep = "")))

lapply(posB, function(x) if(tmp[x] == 'integer') stop(paste(names(elements)[x], ' data not available. see functions addDnase, addTransEqtl, addCisEqtl, addGrasp, addgwasCatalog, addPtv, addSqtlAltrans, addSqtlSeek', sep = ""))) 


# Pull snps for each relevant element
if(1 %in% c(posA, posB)){lapply(mgl, function(x) unique(as.character(do.call(rbind, x[[13]])[,2]))) -> elements[[1]]}

if(2 %in% c(posA,posB)){
	elements[[2]] <- list()
	for(i in 1:length(mgl)){
		a1 <- listSeparate(mgl, 14)[[i]]
		a1 <- a1[which(unlist(lapply(a1, class)) == 'data.frame')]
		if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) != 1 ){
			a1 <- a1[which(unlist(lapply(a1, function(x) dim(x)[1])) > 0)]
			a1 <- unique(do.call(rbind,a1))
			if(dim(a1)[2] == 12){
				elements[[2]][[i]] <- a1[,2]
				message("Warning: It maybe preferable to convert SNP Ids to rs identifiers where possible.  See fixSnpIds.")}
			if(dim(a1)[2] == 13){
				elements[[2]][[i]] <- a1[,13]}}
}}

if(3 %in% c(posA, posB)){
lapply(mgl, function(x) unique(do.call(rbind, x[[15]])[,1])) -> elements[[3]]}

if(4 %in% c(posA, posB)){
lapply(mgl, function(x) unique(do.call(rbind, x[[16]])[,1])) -> elements[[4]]}

if(5 %in% c(posA, posB)){
lapply(mgl, function(x) unique(do.call(rbind, x[[17]])[,1])) -> elements[[5]]}

if(6 %in% c(posA, posB)){
lapply(mgl, function(x) unique(x[[18]][,22])) -> elements[[6]]
lapply(elements[[6]], function(x) unlist(stringr::str_split(x, pattern = ';'))) -> elements[[6]]}

if(7 %in% c(posA, posB)){
lapply(mgl, function(x) unique(x[[19]][,7])) -> elements[[7]]}

if(8 %in% c(posA, posB)){
lapply(mgl, function(x) unique(x[[12]][which(x[[12]][,13] == 'imbalanced_(5%_FDR)' | x[[12]][,13] == 'imbalanced_(0.1%_FDR)'),4])) -> elements[[8]]}

# Pulling elements of interest
# group A
daneA <- list()
for (i in 1:length(mgl)){
	daneA[[i]] <- as.list(1:length(posA))
	for(j in 1:length(posA)){
		if (is.null(elements[[posA[j]]][[i]]) == F){
		daneA[[i]][[j]] <- elements[[posA[j]]][[i]]}
		if (is.null(elements[[posA[j]]][[i]]) == T){
		daneA[[i]][[j]] <- NA
		}
	}
	daneA[which(unlist(lapply(daneA, length)) == 0)] <- NA
	names(daneA[[i]]) <- names(elements)[posA]
}
names(daneA) <- names(mgl)

# group B
daneB <- list()
for (i in 1:length(mgl)){
	daneB[[i]] <- list()
	for(j in 1:length(posB)){
		if (is.null(elements[[posB[j]]][[i]]) == F){
		daneB[[i]][[j]] <- elements[[posB[j]]][[i]]}
		if (is.null(elements[[posB[j]]][[i]]) == T){
		daneB[[i]][[j]] <- NA
		}
	}
	daneB[which(unlist(lapply(daneB, length)) == 0)] <- NA
	names(daneB[[i]]) <- names(elements)[posB]
}
names(daneB) <- names(mgl)

# Determing SNPs in both group A and group B
overlap <- list()
for(i in 1:length(mgl)){
	x = unique(unlist(daneA[[i]]))
	y = unique(unlist(daneB[[i]]))
	overlap[[i]] <- as.character(x[which(x %in% y)])
}
names(overlap) <- names(mgl)

# Remove NA and 'character(0)'
overlap[which(unlist(lapply(overlap, length)) == 0)] <- 'None'
overlap[which(unlist(lapply(overlap, is.na)) == T)] <- 'None'

# Restructuring so overlapping SNPs can be written as a csv
doc <- lapply(1:length(mgl), function(x) cbind(names(overlap)[[x]], overlap[[x]]))
doc <- doc[which(unlist(lapply(doc, function(x) dim(x)[2])) == 2)]
do.call(rbind, doc) -> doc

# Saving if flag
if (saveFile == TRUE)
write.csv(doc, paste('Overlap_', paste(snpsA, collapse = '_'), '_With_', paste(snpsB, collapse = '_'), '.csv', sep = ''))

all <- list(overlap = overlap, snpsA = snpsA, snpsB = snpsB)

#message('\nPlease cite datasources as appropriate. \n \n Element 12 - DNAse: \n Maurano, Matthew T, Eric Haugen, Richard Sandstrom, Jeff Vierstra, Anthony Shafer, Rajinder Kaul, and John A Stamatoyannopoulos. 2015. “Large-Scale Identification of Sequence Variants Influencing Human Transcription Factor Occupancy in Vivo.” Nature Genetics 47 (12): 1393–1401. doi:10.1038/ng.3432.\n \n Elements 13-17 GTEx: \n The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n \n Element 15 - sqtlSeek: \n (1)Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 5:4698 doi: 10.1038/ncomms5698 (2014) \n (2) "Multi-tissue analysis of gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014. \n \n  Element 16 - sqtlAltrans: \n (1) "Multi-tissue analysis of gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014.\n (2) Ongen, H., & Dermitzakis, E. T. (2015). Alternative Splicing QTLs in European and African Populations. American Journal of Human Genetics, 97(4), 567–575. http://doi.org/10.1016/j.ajhg.2015.09.004 \n \n Element 17 - pqtl: \n Rivas MA, Pirinen M, Conrad DF, et al. Impact of predicted protein-truncating genetic variants on the human transcriptome. Science (New York, NY). 2015;348(6235):666-669. doi:10.1126/science.1261877. \n \n Element 18 - NHGRI-EBI GWAS Catalog: \n Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H.The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006 \n \n Element 19 - GRASP: \n (1) Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390 genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP Build 2.0.0.0\n (2) Carey V (2016). grasp2db: grasp2db, sqlite wrap of GRASP 2.0. R package version 0.1.14.\n')

return(all)

}

