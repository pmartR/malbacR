# pmart_amide --------------------------------------------------------

#' pmartR-friendly "Amide" Data
#'
#' These data contains information for 500 features and 642 samples. These data
#' originated from the R package "WaveICA2.0" (see \code{\link[WaveICA2.0]{Amide_data}} for 
#' details) which initially contained 6461 features and 642 samples with three batches.
#'
#' @format A metabData object (see
#'   \code{\link[pmartR]{as.metabData}} for details)
#' \describe{
#'   \item{e_data}{a \eqn{500 \times 643} data frame of expression data.
#'   Each row corresponds to data for each metabolite}
#'   \item{f_data}{a data frame with \eqn{642 \times 4} data frame of feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names, the injection order, the group information, and the batch information.}
#' }
#' @rdname pmart_amide
#' @name pmart_amide
NULL

# pmart_amideFilt --------------------------------------------------------

#' Filtered pmartR-friendly "Amide" Data
#'
#' These data originates from the pmart_amide dataset. However, all necessary
#' data transformations and filters have been applied to the data such that
#' all batch correction methods that in theory can be applied to this dataset can
#' without any additional alterations.
#'
#' @format A metabData object (see
#'   \code{\link[pmartR]{as.metabData}} for details)
#' \describe{
#'   \item{e_data}{a \eqn{500 \times 643} data frame of expression data.
#'   Each row corresponds to data for each metabolite}
#'   \item{f_data}{a data frame with \eqn{642 \times 4} data frame of feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names, the injection order, the group information, and the batch information.}
#' }
#' @rdname pmart_amideFilt
#' @name pmart_amideFilt
NULL

# pmart_mix --------------------------------------------------------

#' pmartR-friendly "mix" Data
#'
#' These data contains information for 46 features and 42 samples. These data
#' originated from the R package "crmn" (see \code{\link[crmn]{mix}} for 
#' details).
#'
#' @format A metabData object (see
#'   \code{\link[pmartR]{as.metabData}} for details)
#' \describe{
#'   \item{e_data}{a \eqn{46 \times 42} data frame of expression data.
#'   Each row corresponds to data for each metabolite}
#'   \item{f_data}{a data frame with \eqn{42 \times 2} data frame of feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names and for the batch information}
#'   \item{e_meta}{a data frame with \eqn{46 \times 7} data frame of meta information.
#'   Each row corresponds to meta information for a metabolite with a column giving the metabolite, marking
#'   information, internal standards information, synonyms, RI value, query, and known entries}
#' }
#' @rdname pmart_mix
#' @name pmart_mix
NULL

# pmart_mixFilt --------------------------------------------------------

#' Filtered pmartR-friendly "mix" Data
#'
#' These data originates from the pmart_mix dataset. However, all necessary
#' data transformations and filters have been applied to the data such that
#' all batch correction methods that in theory can be applied to this dataset can
#' without any additional alterations.
#'
#' @format A metabData object (see
#'   \code{\link[pmartR]{as.metabData}} for details)
#' \describe{
#'   \item{e_data}{a \eqn{46 \times 42} data frame of expression data.
#'   Each row corresponds to data for each metabolite}
#'   \item{f_data}{a data frame with \eqn{42 \times 2} data frame of feature data.
#'   Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
#'   column names and for the batch information}
#'   \item{e_meta}{a data frame with \eqn{46 \times 7} data frame of meta information.
#'   Each row corresponds to meta information for a metabolite with a column giving the metabolite, marking
#'   information, internal standards information, synonyms, RI value, query, and known entries}
#' }
#' @rdname pmart_mixFilt
#' @name pmart_mixFilt
NULL