#' Returns gene names with a given SNP
#'
#' \code{makePhenotypeSearch} returns a character vector of gene names. 
#'
#' Of interest maybe groups of genes that have the same associated SNP 
#'
#' @family output
#'
#' @seealso \code{\link{makePhenotypes}}
#'
#' @param mgl List; see \code{\link{buildFromNames}}, \code{\link{buildFromRegion}}, or \code{\link{buildFromEnsgs}}
#'
#' @param term Character vector representing phenotype of interest 
#'
#' @param snp List produced by \code{\link{makePhenotypes}}
#'
#' @param saveFile A logical flag indicating whether a csv file ('SnpSearch_[term].csv') should be saved in the current directory
#'
#' @examples
#' makeSnps(myMgl, saveFile = FALSE) -> mySnps
#' makeSnpSearch(myMgl, term = '10_51518910_G_A_b37', snp = mySnps, saveFile = FALSE) -> mySnpSearch
#'
#'@export
#'@importFrom utils write.csv
 
makeSnpSearch <- function(mgl, term = c(''), snp, saveFile = FALSE){

# use to display genes with certain GOID
as.character(snp[[3]][which(snp[[3]][,1] == term),2]) -> snames

# Saving if flag
if (saveFile == TRUE)
write.csv(snames, paste('SnpSearch_', stringr::str_replace_all(term, " ", "_"), '.csv', sep = ""))

#message('\nPlease cite datasources as appropriate. \n \n Element 12 - DNAse: \n Maurano, Matthew T, Eric Haugen, Richard Sandstrom, Jeff Vierstra, Anthony Shafer, Rajinder Kaul, and John A Stamatoyannopoulos. 2015. “Large-Scale Identification of Sequence Variants Influencing Human Transcription Factor Occupancy in Vivo.” Nature Genetics 47 (12): 1393–1401. doi:10.1038/ng.3432.\n \n Elements 13-17 GTEx: \n The Genotype-Tissue Expression (GTEx) Project was supported by the Common Fund  of the Office of the Director of the National Institutes of Health. Additional funds were provided by the NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. Donors were enrolled at Biospecimen Source Sites funded by NCISAIC-Frederick, Inc. (SAIC-F) subcontracts to the National Disease Research Interchange (10XS170), Roswell Park Cancer Institute (10XS171), and Science Care, Inc. (X10S172). The Laboratory, Data Analysis, and Coordinating Center (LDACC) was funded through a contract (HHSN268201000029C) to The Broad Institute, Inc. Biorepository operations were funded through an SAIC-F subcontract to Van Andel Institute (10ST1035). Additional data repository and project management were provided by SAIC-F (HHSN261200800001E). The Brain Bank was supported by a supplements to University of Miami grants DA006227 & DA033684 and to contract N01MH000028. Statistical Methods development grants were made to the University of Geneva (MH090941 & MH101814), the University of Chicago (MH090951, MH090937, MH101820, MH101825), the University of North Carolina - Chapel Hill (MH090936 & MH101819), Harvard University (MH090948), Stanford University (MH101782), Washington University St Louis (MH101810), and the University of Pennsylvania (MH101822). The data used for the analyses described in this manuscript were obtained from: [insert, where appropriate] the GTEx Portal on MM/DD/YY and/or dbGaP  accession number phs000424.vN.pN  on MM/DD/YYYY.\n \n Element 15 - sqtlSeek: \n (1)Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun. 5:4698 doi: 10.1038/ncomms5698 (2014) \n (2) "Multi-tissue analysis of gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014. \n \n  Element 16 - sqtlAltrans: \n (1) "Multi-tissue analysisof gene regulation in a human population sample: the Genotype-Tissue Expression (GTEx) pilot study", The GTEx Consortium, under review at Science, Dec. 2014.\n (2) Ongen, H., & Dermitzakis, E. T. (2015). Alternative Splicing QTLs in European and African Populations. American Journal of Human Genetics, 97(4), 567–575. http://doi.org/10.1016/j.ajhg.2015.09.004 \n \n Element 17 - pqtl: \n Rivas MA, Pirinen M, Conrad DF, et al. Impact of predicted protein-truncating genetic variants on the human transcriptome. Science (New York, NY). 2015;348(6235):666-669. doi:10.1126/science.1261877. \n \n Element 18 - NHGRI-EBI GWAS Catalog: \n Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, Flicek P, Manolio T, Hindorff L, and Parkinson H.The NHGRI GWAS Catalog, a curated resource of SNP-trait associations. Nucleic Acids Research, 2014, Vol. 42 (Database issue): D1001-D1006 \n \n Element 19 - GRASP: \n (1) Leslie R, O’Donnell CJ, Johnson AD (2014) GRASP: analysis of genotype-phenotype results from 1,390 genome-wide association studies and corresponding open access database. Bioinformatics 30(12), i185-94. GRASP Build 2.0.0.0\n (2) Carey V (2016). grasp2db: grasp2db, sqlite wrap of GRASP 2.0. R package version 0.1.14.\n')

return(snames)

}

