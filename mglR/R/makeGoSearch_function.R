#' Returns gene names with a given GO term
#'
#' \code{makeGoSearch} returns a character vector of gene names. 
#'
#' Of interest maybe groups of genes that have the same GO term.
#'
#' @family output
#'
#' @seealso \code{\link{makeGo}}
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param term Character vector representing GO term of interest 
#'
#' @param go List produced by \code{\link{makeGo}}
#'
#' @param saveFile A logical flag indicating whether a csv file ('GoSearch_[term].csv') should be saved in the current directory
#'
#' @examples
#' makeGo(myMgl, saveFile = TRUE) -> myGo
#' makeGoSearch(myMgl, term = 'small molecule metabolic process', 
#'     go = myGo, saveFile = TRUE) -> myGoSearch
#'
#'@export
#'@importFrom utils write.csv
 
makeGoSearch <- function(mgl, term = c(''), go, saveFile = TRUE){

# use to display genes with certain GOID
as.character(go[[2]][which(go[[2]][,1] == term),2]) -> gnames

# Saving if flag
if (saveFile == TRUE){
write.csv(gnames, paste('GoSearch_', stringr::str_replace_all(term, " ", "_"), '.csv', sep = ""))}

#message('Please cite biomaRt: \n (1) Durinck S, Spellman P, Birney E and Huber W (2009). “Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt.” Nature Protocols, 4, pp. 1184–1191. \n(2) Durinck S, Moreau Y, Kasprzyk A, Davis S, De Moor B, Brazma A and Huber W (2005). “BioMart and Bioconductor: a powerful link between biological databases and microarray data analysis.”Bioinformatics, 21, pp. 3439–3440\n')

return(gnames)

}

