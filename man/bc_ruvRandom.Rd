% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bc_ruvRandom.R
\name{bc_ruvRandom}
\alias{bc_ruvRandom}
\title{RUV-random normalized adjusted object}
\usage{
bc_ruvRandom(omicsData, nc_cname, nc_val, k = 3)
}
\arguments{
\item{omicsData}{an object of the class 'pepData', 'proData', 'metabData',
'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
\code{\link{as.proData}}, \code{\link{as.metabData}},
\code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.}

\item{nc_cname}{character string giving the name of the column in omicsData$e_meta
that contains the factor variable indicating whether the sample is a negative control or not}

\item{nc_val}{value giving the name/value in the nc_cname column that indicates
that the molecule is a negative control or not}

\item{k}{number of factors of unwanted variation, the default is 3, if the only
factor of unwanted variation is the batch effect, we recommend setting k to be the number
of batches in the data if visual inspection is not an option}
}
\value{
Object of same class as omicsData that has been undergone
  RUV-random normalization
}
\description{
This function returns a pmart object that has been undergone for RUV-random
normalization regarding batch effect correction
}
\examples{
library(malbacR)
library(pmartR)
data("pmart_mix")
pmart_mix <- edata_transform(pmart_mix,"log2")
pmart_mix <- group_designation(pmart_mix,main_effects = "BatchNum",batch_id = "BatchNum")
mix_ruv <- bc_ruvRandom(omicsData = pmart_mix, nc_cname = "tag",nc_val = "IS", k = 3)

}
\author{
Damon Leach
}
