#hack to avoid "no visible global function definition"
utils::globalVariables(c("TISSUE", "REF_RATIO", "VARIANT_ID", "VARIANT_ANNOTATION"))

#' Plot GTEx AEI ratios 
#'
#' \code{makeCorGene} plots AEI ratios person by person and returns a dataframe with allelic ratios for a given gene
#'
#' Uses ggplot function to plot allelic ratios published by GTEx.  Each individual is plotted a separate graph. Tissues are displayed on the x-axis and corresponding ratios on the y.  SNPs are displayed in different colors and annotations in different shapes.Option to save a summarytable (csv) of allelic ratios.  Only one gene at a time.
#'
#' @family output
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param gene String indicating gene to be considered (note gene name must match names(mgl)) 
#'
#' @param saveFile A logical flag indicating whether a csv file ('AeiPlot.csv') should be saved in the current directory. 
#'
#' @examples
#' \dontrun{makeAeiPlot(myMgl, gene = c('NCOA4'), saveFile = FALSE)}
#'
#'@export
#'@importFrom grDevices topo.colors pdf dev.off
#'@importFrom utils write.table
#'@import ggplot2
 
makeAeiPlot <- function(mgl, gene = c(), saveFile = FALSE){

x <- which(names(mgl) == gene)

# Pull data from mgl
dane <- do.call(rbind,mgl[[x]][[20]])

# Add tissue as a column name
cbind(stringr::word(rownames(dane), 1, sep = stringr::fixed('.')), dane) -> dane
colnames(dane)[1] <- 'TISSUE'

# Make ratio numeric
dane$REF_RATIO <- as.numeric(as.character(dane$REF_RATIO))

# Add color 
snps <- unique(dane$VARIANT_ID)
cbind(snps, topo.colors(length(snps))) -> colors
cbind(colors[match(dane$VARIANT_ID, colors[,1]),2], dane) -> dane
dane[,1] <- as.character(dane[,1])
colnames(dane)[1] <- 'COLOR'

# Add shape
locs <- unique(dane$VARIANT_ANNOTATION)
shps <- 15:25
cbind(locs, shps[1:length(locs)])-> shapes
cbind(shapes[match(dane$VARIANT_ANNOTATION, shapes[,1]),2], dane) -> dane
dane[,1] <- as.numeric(as.character(dane[,1]))
colnames(dane)[1] <- 'SHAPE'

# Split into list by person
peeps <- unique(dane$SUBJECT_ID)
lapply(peeps, function(y) dane[which(dane$SUBJECT_ID == y),]) -> dane
names(dane) <- peeps

# Plot by person

#requireNamespace("ggplot2", quiet = TRUE)

pdf(paste('AeiPlot_', gene, '.pdf', sep = ""))
lapply(1:length(peeps), function(i) 
	print(
	ggplot(aes(x = TISSUE, y = REF_RATIO), data = dane[[i]]) + 
	ylim(0,1) + 
	geom_point() + 
	geom_point(aes(color = VARIANT_ID, shape = VARIANT_ANNOTATION), size = 5) + 
	scale_color_manual(values = c(dane[[i]][,2])) +  
	scale_shape_manual(values = c(dane[[i]][,1])) + 
	labs(title = peeps[i]) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	))
dev.off()

# Save
do.call(rbind, dane) -> dane
if (saveFile == TRUE){
write.table(dane, paste('AeiPlot', gene, '.csv', sep = ""), sep = ',', row.names = FALSE)}

#message('Please acknowledge: The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n')

return(dane)

}



