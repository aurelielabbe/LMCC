#' data_methyl
#'
#' Extracted from the methylation data set of Ribaud et al. from site 119521761 to 119527907
#'
#'
#' \itemize{
#' \item ID.source are the ID of source hosts
#' \item ID.recep are the ID of receiver hosts
#'   \item y are the vector of transmission probabilities source -> receiver
#'   \item other columns are the factors
#' }
#'
#' @docType data
#' @keywords equineDiffFactors
#' @name data_methyl
#' @usage data(data_methyl)
#' @format A list with three elements: Y the matrix with the methylation level for each sample and each site,
#' sites the genomic position of sites,
#' X the covariate with the cell type
#' @author Melina Ribaud
#' @references Melina Ribaud, Aurélie Labbe and Karim Oualkacha.
#' Imputation in genetic methylation studies: A linear model of coregionalization (LMC) with informative covariates.
#' 2022. hal-00000000
NULL
