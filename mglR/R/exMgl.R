#' exMgl
#'
#' \code{exMgl} returns example mgl list for vignette myMgl
#'
#' Example of fully built list for five candidate genes(note: element 20 - AEI is still missing as this data requires IRB approval and access via dbGaP).
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'myMgl.rda'
#'
#' @examples
#' exMgl() -> myMgl
#'
#'@export

exMgl <- function(saveDownload = F){
	# Download data
	download.file('https://github.com/klh64/mglR/raw/master/mglR/data/myMgl.rda', 'myMgl.rda')	
	# Read into R
	load('myMgl.rda')
	if(saveDownload == FALSE){
		unlink('myMgl.rda')}
	return(myMgl)
}