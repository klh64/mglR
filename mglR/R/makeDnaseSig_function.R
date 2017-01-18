#' Filter Dnase (Maurano) Results
#'
#' \code{makeDnaseSig} returns a list with length corresponding to the number of candidate genes.  Each element is a dataframe of all the Dnase results filtered for those that are significant. 
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @examples
#' myMgl <- makeDnaseSig(myMgl)
#'
#'@export
 
makeDnaseSig <- function(mgl){

# stop process if none of the SNP based elements have been filled in


if(unique(unlist(lapply(mgl, function(x) class(x[[12]])))) == 'integer')
	stop('data not available. see function addDnase') 

res <- list()
for (i in 1:length(mgl)){
	res[[i]] <- unique(mgl[[i]][[12]][which(mgl[[i]][[12]][,13] == 'imbalanced_(5%_FDR)' | mgl[[i]][[12]][,13] == 'imbalanced_(0.1%_FDR)'),4])
	}

names(res) <- names(mgl)

#message('Please cite: Maurano, Matthew T, Eric Haugen, Richard Sandstrom, Jeff Vierstra, Anthony Shafer, Rajinder Kaul, and John A Stamatoyannopoulos. 2015. “Large-Scale Identification of Sequence Variants Influencing Human Transcription Factor Occupancy in Vivo.” Nature Genetics 47 (12): 1393–1401. doi:10.1038/ng.3432.\n')
	
return(res)

}

