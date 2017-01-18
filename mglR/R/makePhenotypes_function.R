#' Make list of GWAS-based traits
#'
#' \code{makePhenotypes} returns a list with four elements summarizing the traits or phenotypes associated with each gene via GWAS results reported in both the GWAS catalog and GRASP.  The first element (phenRes) is a list with length corresponding to the number of candidate genes and includes all GWAS-based phenotypes split by source: NHGRI-EBI GWAS Catalog and GRASP.  The second element (phenList) collpases the first so only phenotypes for each gene are listed, regardless of the source.  The third element (phenTable) is a dataframe with two columns: Phenotypes and gene names.  The fourth element (phenCount) is a table with the number of elements corresponding to the number of unique phenotypes in the genelist - reported for each phenotype is the number of times it appears.  It is sorted in descending order.  The structure is similar to the \code{\link{makeGo}} and \code{\link{makeSnps}} functions.
#'
#' This gives a summary of all traits associated with the gene in the NHGRI-EBI GWAS catalog and GRASP.
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param saveFile A logical flag indicating whether a csv file ('Phenotypes.csv') should be saved in the current directory. 
#'
#' @examples
#' myMgl <- makePhenotypes(myMgl, saveFile = FALSE)
#'
#'@export
#'@importFrom utils write.table write.csv
 
makePhenotypes <- function(mgl, saveFile = FALSE){

# stop process if none of the elements of interest have been filled in

tmp <- c(unique(unlist(lapply(mgl, function(x) class(x[[18]])))), unique(unlist(lapply(mgl, function(x) class(x[[19]])))))

if(length(unique(tmp)) == 1){
	if(unique(tmp) == 'integer')
	stop('Trait data not available. See function listElements to check which elements are present.  See functionss addGrasp and addgwasCatalog to add pertinent information.') 
}

### get only those elements that are filled in
dane <-as.list(1:length(mgl))
for(i in 1:length(mgl)){
	dane[[i]] <- as.list(1:2)
	names(dane[[i]]) <- c('GwasCatalog', 'GRASP')
}
names(dane) <- names(mgl)

for(i in 1:length(mgl)){
if(tmp[1] != "integer"){
	dane[[i]][[1]] <- unique(mgl[[i]][[18]][,8])}
if(tmp[2] != "integer"){
	dane[[i]][[2]] <- unique(mgl[[i]][[19]][,11])}
}

# Add message about no phenotypes found for those elements that are empty
for (i in 1:length(mgl)){
	for(x in which(unlist(lapply(dane[[i]], length)) == 0)){
		dane[[i]][[x]] <- 'None'
	}
}

# Add message about data missing from mgl for those elements that are not filled in	
for(x in 1:2){
	if(tmp[x] == 'integer'){
		for(i in 1:length(mgl)){
			dane[[i]][[x]] <- 'Data missing from mgl list'
		}
	}
}

### Collapse so just get SNPs and format
daneCollapsed <- as.list(1:length(mgl))
names(daneCollapsed) <- names(mgl)
for(i in 1:length(mgl)){
unique(unlist(dane[[i]])) -> x
if (sum(is.na(x)) != 0) {
	x <- x[-c(which(is.na(x) == T))]
}
x <- x[x != ""]
if (length(grep('None', x)) != 0){
	x <- x[-c(which(x == 'None'))]
}
if (length(grep('Data missing from mgl list', x)) != 0){
	x <- x[-c(which(x == 'Data missing from mgl list'))]
}
x <- sort(x)
daneCollapsed[[i]] <- x
}


### making table with first column as phenotype and second as gene
# Removing those genes that have no go information i.e. an empty elementy five
tab <- daneCollapsed
for (x in 1:length(tab)){
	if (length(tab[[x]]) > 0) {
		tab[[x]] <- cbind(tab[[x]], names(tab)[x])
		}
		}
# Making it into a datable instead of a list 
do.call(rbind, tab) -> phenTable
# Adding names
colnames(phenTable) <- c("Phenotype", "Gene")
rownames(phenTable) <- NULL

### Tabling the number of times a phenotype appears
table(phenTable[,1]) -> phenCount
phenCount <- sort(phenCount, decreasing = T)

### Saving		
if(saveFile == TRUE){
write.table("GWAS-based Traits", 'Phenotypes.csv', sep = ',', col.names = FALSE, row.names = FALSE)
for(i in 1:length(dane)){
	write.table(names(mgl)[i], "Phenotypes.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table(matrix(dane[[i]], ncol = 5), "Phenotypes.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table("", "Phenotypes.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
write.csv(phenTable, 'phenTable.csv')
write.csv(phenCount, 'phenCount.csv')
}
}

phen <- list(phenRes = dane, phenList = daneCollapsed, phenTable = phenTable, phenCount = phenCount)

#message('\nPlease cite the NHGRI-EBI GWAS Catalog or GRASP as appropriate. \n \n NHGRI-EBI GWAS Catalog: \n Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H.The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006 \n \n GRASP: \n (1) Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390 genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP Build 2.0.0.0\n (2) Carey V (2016). grasp2db: grasp2db, sqlite wrap of GRASP 2.0. R package version 0.1.14.\n')

return(phen)

}

