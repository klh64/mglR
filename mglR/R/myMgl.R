#' ex
#'
#' \code{ex} returns example mgl list
#'
#' Example of fully built list for five candidate genes(note: element 20 - AEI is still missing as this data requires IRB approval and access via dbGaP).
#'
#' @param saveDownload A logical vector indicating if the data should be saved as 'myMgl.RData'
#'
#' @examples
#' ex() -> myMgl
#'
#'@export

ex <- function(saveDownload = F){
	# Download data
	download.file('http://www.nature.com/ng/journal/v47/n12/extref/ng.3432-S5.txt', temp)	
	# Read into R
	load('temp')
	if(saveDownload == TRUE){
		save(myMgl, file = 'myMgl.RData')}
	unlink(temp)		
}