# malbacR
<!-- badges: start -->
[![DOI](https://zenodo.org/badge/69275428.svg)](https://zenodo.org/badge/latestdoi/69275428)
[![R-CMD-check](https://github.com/pmartR/pmartR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pmartR/pmartR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This R package provides functionality for batch correction methods for mass spectrometry (MS) omics data. This includes Pareto Scaling, Power Scaling, Range Scaling, ComBat, EigenMS, RUV-random, NOMIS, QC-RLSC, WaveICA2.0, SERRF, and TIGER). Additional functions include an imputation function (using expectation maximization algorithms) and a function to filter data for TIGER batch correction. Example data is also provided in the package.

## Installation:

``` r
devtools::install_github("pmartR/malbacR")
```


### Problems with rcppArmadillo and gfortran on mac

The R package pmartR is a dependency for malbacR and will install when installing malbacR if not already installed. There is a problem that causes pmartR to fail compiling cpp code, which has something to do with rcppArmadillo and certain installations of gfortran.  See these posts that try to explain the issue:  [1](https://stackoverflow.com/questions/64992467/mac-clang-installation-seems-to-override-gcc-install) [2](https://stackoverflow.com/questions/29992066/rcpp-warning-directory-not-found-for-option-l-usr-local-cellar-gfortran-4-8/29993906#29993906) [3](https://community.rstudio.com/t/setting-up-travis-ci-on-linux-with-an-r-package-that-uses-rcpparmadillo/53910/3).  Two solutions we have found:

1.  Install gfortran from a recommended source (not homebrew): 
    - [This CRAN-approved resource for build tools on mac](https://mac.r-project.org/tools/) lists two versions of gfortran and how to install them.
    - On Catalina 10.15.7 I downloaded and installed gfortran 8.2 from the link provided in [this blog post](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/#google_vignette)  
2.  When using the homebrew gfortran installation, add the line **FLIBS = -L\`gfortran -print-file-name=libgfortran.dylib | xargs dirname\`** to ~/.R/Makevars (a plain text file with no extention)


## Tutorial:

To get started, see the package documentation and function reference located [here](https://pmartr.github.io/malbacR/).

## Data:

Example peptide (both filtered and unfiltered versions) of two metabolite datasets. The first data set (pmart_amide/pmart_amideFilt) is a subset of a dataset originally came from the the R package __WaveICA2.0__ which can be found [here](https://github.com/dengkuistat/WaveICA_2.0). The second dataset (pmart_mix/pmart_mixFilt) originally came from the R package __crmn__ which can be found [here](https://cran.rstudio.com/web/packages/crmn/index.html).

## Citation:

To cite this package, please the following:

Degnan, D. J.; Stratton, K. G.; Richardson, R.; Claborne, D.; Martin, E. A.; Johnson, N. A.; Leach, D.; Webb-Robertson, B.-J. M.; Bramer, L. M. PmartR 2.0: A Quality Control, Visualization, and Statistics Pipeline for Multiple Omics Datatypes. J. Proteome Res. 2023, 22 (2), 570–576. https://doi.org/10.1021/acs.jproteome.2c00610.

**BibTex:**

```
@article{degnan2023pmartr,
  title={pmartR 2.0: A Quality Control, Visualization, and Statistics Pipeline for Multiple Omics Datatypes},
  author={Degnan, David J and Stratton, Kelly G and Richardson, Rachel and Claborne, Daniel and Martin, Evan A and Johnson, Nathan A and Leach, Damon and Webb-Robertson, Bobbie-Jo M and Bramer, Lisa M},
  doi = {10.1021/acs.jproteome.2c00610},
  journal={Journal of Proteome Research},
  year={2023},
  publisher={ACS Publications}
}) 
```
