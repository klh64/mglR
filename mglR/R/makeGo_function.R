#' Returns a list with Gene Ontology data
#'
#' \code{makeGo} returns a list with three elements summarizing Gene Ontology (GO) terms.  The first element (goRes) is the GO results subsetted from the mgl.  The second element (goTable) is a dataframe with two columns: GO terms and gene names.  The third element (goCount) is a table with the number of elements corresponding to the number of unique GO terms in the genelist - reported for each GO term is the number of times it appears.  It is sorted in descending order.The structure is similar to the \code{\link{makePhenotypes}} and \code{\link{makeSnps}} functions.
#'
#' Provides a brief summary of GO terms that have been associate with candidate genes. Of interest maybe groups of genes that have the same GO term see \code{makeGoSearch}.
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param saveFile A logical flag indicating whether two csv files: 'GoTable.csv' and 'GoCount.csv' should be saved in the current directory
#'
#' @examples
#' makeGo(myMgl, saveFile = TRUE) -> myGo
#'
#'@export
#'@importFrom utils write.table write.csv
 
makeGo <- function(mgl, saveFile = TRUE){

# Pulling data from mgl
goRes <- NULL
for (x in 1:length(mgl)){
	goRes[[x]] <- mgl[[x]][[5]]}
names(goRes) <- names(mgl)

# Removing those genes that have no go information i.e. an empty elementy five
goTable <- goRes
for (x in 1:length(goTable)){
	if (dim(mgl[[x]][[5]])[1] > 0) {
		goTable[[x]] <- cbind(goTable[[x]], names(goTable)[x])
		}
		}
		
# Making it into a datable instead of a list 
do.call(rbind, goTable) -> goTable

# Adding names
rownames(goTable) <- NULL
colnames(goTable)[1] <- 'GoTerm'
colnames(goTable)[2] <- 'Gene'

# Tabling the number of times a GO term appears
table(goTable[,1]) -> goCount
goCount <- sort(goCount, decreasing = T)

# Saving if flag
if (saveFile == TRUE){
write.table("GO Results", 'goRes.csv', sep = ',', col.names = FALSE, row.names = FALSE)
for(i in 1:length(goRes)){
	write.table(names(mgl)[i], "goRes.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table(matrix(goRes[[i]], ncol = 5), "goRes.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
	write.table("", "goRes.csv", sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
write.csv(goTable, 'goTable.csv')
write.csv(goCount, 'goCount.csv')
}
}

# Compiling into one object
go <- list(goRes = goRes, goTable = goTable, goCount = goCount)

#message('Please cite biomaRt: \n (1) Durinck S, Spellman P, Birney E and Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, pp. 1184–1191. \n(2) Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A and Huber W (2005). “BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.”Bioinformatics, 21, pp. 3439–3440\n')

return(go)

}

