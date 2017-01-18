#' Returns information for SNPs that appear in two user-defined groups for a given gene
#'
#' \code{makeSummary} returns a nested list with a dataframe for each gene, snp, source combination (i.e. myOverlapTable[[gene]][[snp]][[source]]).    
#'
#' For every gene with SNPs that overlap between two different user defined groups (see \code{\link{makeOverlap}}), provides a list of dataframes for each overlapping SNP that details the association whether it be eQTL, splicing, or trait-based.  Provides the details for overlapping SNPs. 
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param overlap A list of overlapping snps made using \code{\link{makeOverlap}}
#'
#' @param saveFile A logical flag indicating whether a csv ('Overlap[groups].csv') should be saved in the current directory
#'
#' @examples
#' makeOverlap(myMgl, snpsA = c('cisEqtls'), snpsB = c('gwasCatalog'), saveFile = TRUE) -> myOverlap
#' makeOverlapTable(myMgl, myOverlap, saveFile = TRUE) -> myOverlapTable
#'
#'@export
#'@importFrom utils write.table
 
makeOverlapTable <- function(mgl, overlap, saveFile = TRUE){
	
# which elements of mgl to pull from
unlist(lapply(overlap$snpsA, function(x) which(names(mgl[[1]]) == x))) -> posA
unlist(lapply(overlap$snpsB, function(x) which(names(mgl[[1]]) == x))) -> posB
c(posA, posB) -> pos

# pulling data
res <- list()
nums <- 1:length(overlap[[1]])
if (length(grep('None', overlap[[1]])) != 0){
	nums <- nums[-c(grep('None', overlap[[1]]))]}
# for each gene where there is an overlapping snp
for (i in 1:length(nums)){
	res[[i]] <- list()
	# for each overlapping snp
	for(j in 1:length(overlap[[1]][[nums[i]]])){
		res[[i]][[j]] <- list()
		# for each element defining group A
		for (k in 1:length(pos)){
			if(pos[k] == 13){
				res[[i]][[j]][[k]] <- do.call(rbind, lapply(mgl[[nums[i]]][[13]], function(x) x[which(x[,2] == overlap[[1]][[nums[i]]][[j]]),]))
			}
			if (pos[k] == 14){
				a1 <- listSeparate(mgl, 14)[[nums[i]]]
				a1 <- a1[which(unlist(lapply(a1, class)) == 'data.frame')]
				if(length(unique(unlist(lapply(a1, function(x) dim(x)[1])))) != 1 ){
					a1 <- a1[which(unlist(lapply(a1, function(x) dim(x)[1])) > 0)]
					a1 <- unique(do.call(rbind,a1))
					if(dim(a1)[2] == 12){
						res[[i]][[j]][[k]] <- a1[which(a1[,2] == overlap[[1]][[nums[i]]][[j]]),]}
					if(dim(a1)[2] == 13){
						res[[i]][[j]][[k]] <- a1[which(a1[,13] == overlap[[1]][[nums[i]]][[j]]),]}		
			}}
			if (pos[k] == 15){
				res[[i]][[j]][[k]] <- do.call(rbind,lapply(mgl[[nums[i]]][[15]], function(x) x[which(x[,1] == overlap[[1]][[nums[i]]][[j]]),]))			
			}
			if (pos[k] == 16){
				res[[i]][[j]][[k]] <- do.call(rbind,lapply(mgl[[nums[i]]][[16]], function(x) x[which(x[,1] == overlap[[1]][[nums[i]]][[j]]),]))			
			}
			if (pos[k] == 17){
				res[[i]][[j]][[k]] <- do.call(rbind,lapply(mgl[[nums[i]]][[17]], function(x) x[which(x[,1] == overlap[[1]][[nums[i]]][[j]]),]))	
			}
			if (pos[k] == 18){
				res[[i]][[j]][[k]] <- mgl[[nums[i]]][[18]][which(mgl[[nums[i]]][[18]][,22] == overlap[[1]][[nums[i]]][[j]]),]
			}
			if (pos[k] == 19){
				res[[i]][[j]][[k]] <- mgl[[nums[i]]][[19]][which(mgl[[nums[i]]][[19]][,7] == overlap[[1]][[nums[i]]][[j]]),]
			}
			if (pos[k] == 12){
				res[[i]][[j]][[k]] <- mgl[[nums[i]]][[12]][which(mgl[[nums[i]]][[12]][,4] == overlap[[1]][[nums[i]]][[j]]),]				
			}
		}
		names(res[[i]][[j]]) <- c(overlap$snpsA, overlap$snpsB)
	}
	names(res[[i]]) <- overlap[[1]][[nums[i]]]	
}
names(res) <- names(mgl)[nums]

# Removing empty slots - a SNP may appear in one but not all categories given to group A or B
# first a list mirroring results with the dimensions
dimensions <- as.list(1:length(res))
for(i in 1:length(res)){
	dimensions[[i]] <- as.list(1:length(res[[i]]))
	for(j in 1:length(res[[i]])){
		dimensions[[i]][[j]] <- lapply(res[[i]][[j]], dim)
		res[[i]][[j]] <- res[[i]][[j]][which(do.call(rbind, dimensions[[i]][[j]])[,1] != 0)]
		for(k in 1:length(res[[i]][[j]])){
			rownames(res[[i]][[j]][[k]]) <- NULL
		}
	}
	names(dimensions[[i]]) <- names(res[[i]])
}
names(dimensions) <- names(res)

# Saving if flag
if (saveFile == TRUE){
fname = paste('OverlapTable_', paste(overlap$snpsA, collapse = '_'), '_With_', paste(overlap$snpsB, collapse = '_'), '.csv', sep = '')
write.table("", fname , sep = ',', col.names = FALSE, row.names = FALSE)
for (i in 1:length(res)){
	for(j in 1:length(res[[i]])){
		for(k in 1:length(res[[i]][[j]])){
			write.table(paste(names(res)[i], names(res[[i]])[j], names(res[[i]][[j]])[k], sep = " | "), fname, sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
			suppressWarnings(write.table(res[[i]][[j]][[k]], fname, row.names = FALSE, col.names = TRUE, sep = ',', append = TRUE))
			write.table("", fname , sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
		}
		write.table("", fname , sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
	}
	write.table("", fname , sep = ',', append = TRUE, col.names = FALSE, row.names = FALSE)
}
}

#message('\nPlease cite datasources as appropriate. \n \n Element 12 - DNAse: \n Maurano, Matthew T, Eric Haugen, Richard Sandstrom, Jeff Vierstra, Anthony Shafer, Rajinder Kaul, and John A Stamatoyannopoulos. 2015. “Large-Scale Identification of Sequence Variants Influencing Human Transcription Factor Occupancy in Vivo.” Nature Genetics 47 (12): 1393–1401. doi:10.1038/ng.3432.\n \n Elements 13-17 GTEx: \n The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n \n Element 15 - sqtlSeek: \n (1)Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 5:4698 doi: 10.1038/ncomms5698 (2014) \n (2) "Multi-tissue analysis of gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014. \n \n  Element 16 - sqtlAltrans: \n (1) "Multi-tissue analysisof gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014.\n (2) Ongen, H., & Dermitzakis, E. T. (2015). Alternative Splicing QTLs in European and African Populations. American Journal of Human Genetics, 97(4), 567–575. http://doi.org/10.1016/j.ajhg.2015.09.004 \n \n Element 17 - pqtl: \n Rivas MA, Pirinen M, Conrad DF, et al. Impact of predicted protein-truncating genetic variants on the human transcriptome. Science (New York, NY). 2015;348(6235):666-669. doi:10.1126/science.1261877. \n \n Element 18 - NHGRI-EBI GWAS Catalog: \n Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H.The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006 \n \n Element 19 - GRASP: \n (1) Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390 genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP Build 2.0.0.0\n (2) Carey V (2016). grasp2db: grasp2db, sqlite wrap of GRASP 2.0. R package version 0.1.14.\n')

return(res)

}

