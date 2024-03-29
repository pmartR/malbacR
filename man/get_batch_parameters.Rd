% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_fns.R
\name{get_batch_parameters}
\alias{get_batch_parameters}
\alias{get_batch_method}
\title{Fetch the Batch Effect Parameters}
\usage{
get_batch_parameters(omicsData)

get_batch_method(omicsData)
}
\arguments{
\item{omicsData}{An object of class pepData, proData, metabData, lipidData,
or nmrData.}
}
\value{
A list for each of the different parameters used to run batch correction

A character vector with the function name for batch correction
}
\description{
Retrieves the values from the batch effect parameters attributes

Retrieves the values from the batch effect method name
}
