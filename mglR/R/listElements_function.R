#' Returns names of elements that are filled in
#'
#' \code{listElements} returns a list of any given element 
#'
#' @family list
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param added A logical flag indicating whether elements that have already been added should be displayed.
#'
#' @examples
#' listElements(myMgl) -> added
#'
#'@export
 
listElements <- function(mgl, added = TRUE){

if (added == TRUE){
tmp <- which(lapply(mgl[[1]], class) != 'integer')}

if (added == FALSE){
tmp <- which(lapply(mgl[[1]], class) == 'integer')}

return(tmp)

}

