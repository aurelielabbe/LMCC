#' lmcc: a Linear Model of Coregionalization with informative Covariates
#' @import data.table
#' @import FastGP
#' @import FastGaSP
#' @importFrom stats optim
#' @importFrom data.table :=
#' @importFrom ZIprop scale_01
#' @description
#' We propose a method to predict the missing methylation levels.
#' This method catches correlation structures among the methylation levels across genome sites and across samples.
#'The regression function linking the methylation level to the covariates
#'is modeled through a linear combination of the covariates together with latent factors.
#'We assume the covariates' and latent factors' effects to be Gaussian random processes (GP).
#' @author Melina Ribaud
#' @references Melina Ribaud, Aurélie Labbe and Karim Oualkacha.
#' Imputation in genetic methylation studies: A linear model of coregionalization (LMC) with informative covariates.
#' 2022. hal-00000000
#'
#' @docType package
#' @name lmcc
NULL
#> NULL


