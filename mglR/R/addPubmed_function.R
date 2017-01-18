#' Adds number of titles from Pubmed to list
#'
#' \code{addGo} returns an 'mgl' list with the sixth element as the number of titles from scientific publications indexed in Pubmed that reference the gene of interest.  Please note that this function is a modification of GetPubMed from the NCBI2R package written by Scott Melville and archived on CRAN (see \url{https://ncbi2r.wordpress.com/})
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param max Maximum number of papers to return (i.e. this value will be returned if the number of papers exceeds it)
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addPubmed(myMgl)}
#'
#'@export

 
addPubmed <- function(mgl, max = 30000){
	
for(i in 1:length(mgl)){

#Set-up URL
	searchterm <- names(mgl)[i]
    searchterm <- gsub(" ", "+", searchterm)
    getURL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", "esearch.fcgi?db=pubmed&term=", 
        searchterm, "&retmax=", max, "&rettype=FASTA", "&tool=NCBI2R&email=ncbi2r@gmail.com", 
        sep = "")

# Pull details
    sep = "\n"
    quiet = TRUE
    retry = 20
    error = TRUE
    thispage <- try(scan(getURL, what = "character", sep = sep, quiet = quiet), silent = TRUE)
    retrycount <- 0
    while(class(thispage) == "try-error" & retrycount < retry){
    	TimeStampA <- Sys.time()
    	while(Sys.time()<TimeStampA+0.5)
    		OnlyForDelay <- 1
    	retrycount <- retrycount+1
    	thispage <- try(scan(getURL, what = 'character', sep = sep, quiet = quiet), silent = TRUE)
    	}
    	if(class(thispage) == 'try-error' & error == TRUE)
    		stop("NCBI2R error: The website was not found.")
    	if(class(thispage) == 'try-error' & error == FALSE)
    		thispage <- "Page not found"

# Calculate number
	mgl[[i]][[6]] <- tryCatch({as.numeric(gsub("^[[:print:]]*<Count>([[:digit:]]+)[[:print:]]*$", "\\1", thispage[grep("<Count>", thispage)][[1]]))}, error = function(e){return(NA)}) 
}
	return(mgl)
}


