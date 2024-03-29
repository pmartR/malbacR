% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{pmart_amide}
\alias{pmart_amide}
\title{pmartR-friendly "Amide" Data}
\format{
A metabData object (see
  \code{\link[pmartR]{as.metabData}} for details)
\describe{
  \item{e_data}{a \eqn{500 \times 643} data frame of expression data.
  Each row corresponds to data for each metabolite}
  \item{f_data}{a data frame with \eqn{642 \times 4} data frame of feature data.
  Each row corresponds to a sample with a column giving the unique sample identifiers found in e_data
  column names, the injection order, the group information, and the batch information.}
}
}
\description{
These data contains information for 500 features and 642 samples. These data
originated from the R package "WaveICA2.0" (see \code{\link[WaveICA2.0]{Amide_data}} for 
details) which initially contained 6461 features and 642 samples with three batches.
}
