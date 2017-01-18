#' Add allele specific expression data to list
#'
#' \code{addAei} returns a list with the twentieth element as dataframe with allele specific data from GTEx. 
#'
#' This gives allele specific expression RNAsequencing data as reported by GTEx for the gene of interest.  See \url{http://www.gtexportal.org/}.  It pulls data based on the ENSG identifier.
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param fpsource A character string of with the filepath where the data has been downloaded
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addAei(myMgl, fpsource = "/Downloads")}
#'
#'@export
#'@importFrom utils write.table read.table
 
 
addAei <- function(mgl, fpsource){

# make vector of ENSGs
ensgs <- list(); for (x in 1:length(mgl)){ensgs[[x]] <- mgl[[x]][[1]][,2]}
ensgs <- unlist(ensgs)

# if NAs exist remove them
if (sum(is.na(ensgs)) != 0) {
	ensgs <- ensgs[-c(which(is.na(ensgs) == T))]
}

# write table to use in grep command
write.table(ensgs, 'ensgs.txt', quote = F, col.names = F, row.names = F)

# compile & run code
code <- paste('cd ', fpsource, ' ; grep CHR ', list.files(fpsource)[1], ' > ', getwd(), '/aei.txt ; grep -r -f ', getwd(), '/ensgs.txt ', ' >> ', getwd(), '/aei.txt', sep = "")
system(code)

# read into R
read.table('aei.txt', header = T, stringsAsFactors = F, sep = '\t') -> aei

# reformat a bit
aei[,1] <- stringr::word(aei[,1], 2, sep = stringr::fixed(':'))

# establish structure in mgl - one element of list for each tissue
tissues <- unique(aei$TISSUE_ID)
for (i in 1:length(mgl)) { 
	mgl[[i]][[20]] <- as.list(1:length(tissues)) 
	names(mgl[[i]][[20]]) <- tissues
}

# add tissue specific data to mgl
for(i in 1:length(mgl)) { 
	for(j in 1:length(tissues)) {
		mgl[[i]][[20]][[j]] <- aei[which(aei$TISSUE_ID == tissues[j] & aei$GENE_ID == mgl[[i]][[1]][1,2]),]
	}
}

unlink('aei.txt')
unlink('ensgs.txt')
	
return(mgl)

}

