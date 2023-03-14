# malbacR

This R package provides functionality for batch correction methods for mass spectrometry (MS) omics data. This includes Pareto Scaling, Power Scaling, Range Scaling, ComBat, EigenMS, RUV-random, NOMIS, QC-RLSC, WaveICA2.0, SERRF, and TIGER). Additional functions include an imputation function (using expectation maximization algorithms) and a function to filter data for TIGER batch correction. Example data is also provided in the package.

## Installation:

``` r
devtools::install_github("pmartR/malbacR")
```

pmart required


It’s been noted, here https://stackoverflow.com/questions/51257009/is-rtools-incompatible-with-r-version-3-5-1 , that there are issues with Rtools for R version 3.5.1.
 
If you run into the same problem when trying to install using `devtools::install_github()` and the suggested fix in the above link does not work, you can clone or download the malbacR library to your computer and install it locally.

It's also been noted, here https://github.com/rcppsmc/rcppsmc/issues/27 , that there are issues with xcode. This manifests itself as an error stating "math.h" not found when trying to install and restart a package. The following fix, at the terminal, has worked for us:

xcode-select --install


## Tutorial:

To get started, see the package documentation and function reference located [here](https://pmartr.github.io/pmartR/index.html).

## Data:

Example peptide (both filtered and unfiltered versions) of two metabolite datasets. The first data set (pmart_amide/pmart_amideFilt) is a subset of a dataset originally came from the the R package __WaveICA2.0__ which can be found [here](https://github.com/dengkuistat/WaveICA_2.0). The second dataset (pmart_mix/pmart_mixFilt) originally came from the R package __crmn__ which can be found [here](https://cran.rstudio.com/web/packages/crmn/index.html).
 

