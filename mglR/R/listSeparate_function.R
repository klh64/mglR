#' Returns an individual element for each gene
#'
#' \code{listSeparate} returns a list of any given element 
#'
#' @family list
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param element A numeric indicating which element to be isolated
#'
#' @examples
#' listSeparate(myMgl, element = 3) -> sep
#'
#'@export
 
listSeparate <- function(mgl, element = c('')){

tmp <- list(); 
for (x in 1:length(mgl)){
	tmp[[x]] <- mgl[[x]][[element]]
	}

names(tmp) <- names(mgl)

return(tmp)

}

