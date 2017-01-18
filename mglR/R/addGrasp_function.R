#hack to avoid "no visible global function definition"
utils::globalVariables(c("chr_hg19", "pos_hg19"))

#' Add GWAS data from GRASP to list
#'
#' \code{addGrasp} returns a list with the nineteenth element as dataframe with GWAS data from the GRASP database \url{grasp.nhlbi.nih.gov/}. 
#'
#' This gives basic information on trait associated variants as reported by the GRASP database for the gene of interest.  It pulls data based on the position (i.e. any SNP that falls between the start and stop position).  Note a wider range can be taken using the \emph{range} flag.   
#'
#' @family elements
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param range Number indicating distance upstream of start and downstream of stop that should be used for filtering SNPs.
#'
#' @examples
#' \dontrun{buildFromRegion(chr = 2, start = 102314000, stop = 103435000) -> myMgl}
#' \dontrun{myMgl <- addGrasp(myMgl, range = 0)}
#'
#'@export
#'@importFrom grasp2db GRASP2
#'@importFrom dplyr tbl filter
#'@importFrom magrittr %>%
 
addGrasp <- function(mgl, range){

grasp2 <- GRASP2()
	
# add to mgl
for (i in 1:length(mgl)) {mgl[[i]][[19]] <- 	
as.data.frame(tbl(grasp2, 'variant') %>% filter(chr_hg19 == mgl[[i]][[3]][[1]] & pos_hg19 <= (mgl[[i]][[3]][[3]] + range) & pos_hg19 >= (mgl[[i]][[3]][[2]] - range)))
	}
		
return(mgl)

}

#


