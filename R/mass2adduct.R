"_PACKAGE"
#> [1] "_PACKAGE"
#'
#' @title
#' Mass Spectrometry Imaging Molecular Adduct Finder
#'
#' @description
#' Tools for calculating mass differences in mass spectrometry data to explore
#' molecular adducts. To accompany Janda et al. (MS in prep.)
#'
#' @details
#' Import mass spectrometry imaging intensity data with \code{\link{msimat}}.
#' 
#' Tabulate mass differences from simple list of mass peaks with
#' \code{\link{massdiff}}.
#' 
#' Bin mass differences into histogram with \code{\link{hist.massdiff}}.
#' 
#' Get pairs of mass peaks corresponding to a specific mass difference with
#' \code{\link{diffGetPeaks}}.
#' 
#' Find the closest-matching mass differences for a set of known molecular
#' adducts with \code{\link{adductMatch}}; conversely, find the closest-
#' matching known molecular adducts for the most-abundant mass differences
#' with \code{\link{topAdducts}}.
#'
#' Find significant correlations between specified pairs of mass peaks in
#' mass spectrometry imaging intensity data using \code{\link{corrPairsMSI}}.
#'
#' Accompanying data: Lists of common molecular \code{\link{adducts}} with their
#' respective formulae and masses
#' 
#' @docType package
#' @name mass2adduct
#' @useDynLib mass2adduct
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline hist plot points image
#' @importFrom stats cor.test ecdf
#' @importFrom utils read.csv
NULL
