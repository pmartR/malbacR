context('Run imputation')

test_that('imputation function to impute the data as expected',{
  
  # Load the reduced peptide data frames ---------------------------------------
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Run as.pepData with agreeable data frames ----------------------------------
  
  # Construct a pepData object with the edata, fdata, and emeta data frames.
  mdata <- pmartR::as.metabData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Metabolite',
                                fdata_cname = 'SampleID',
                                emeta_cname = 'Metabolite')
  attributes(mdata)$data_info$data_scale_orig <- "log2"
  attributes(mdata)$data_info$data_scale <- "log2"
  
  # Run through the potential error messages -----------------------------------
  
  # need to remove molecules that have all NA values
  expect_error(imputation(mdata),
               "At least one row has all NA values")
  # so we remove those values
  molfilt <- pmartR::molecule_filter(mdata)
  mdataFilt <- pmartR::applyFilt(molfilt,mdata)
  
  # object must be an omicsData
  expect_error(imputation(mdataFilt$e_data),
               "omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData'")
  
  # Get the imputation object --------------------------------------------------
  # run the function
  impObj <- imputation(mdataFilt)
  
  # dimensions should be the same (but edata will have one extra column)
  expect_equal(nrow(impObj),
                   nrow(mdataFilt$e_data))
  expect_equal(ncol(impObj) + 1,
                   ncol(mdataFilt$e_data))
  
  # ensure that there are no NA values
  expect_equal(sum(is.na(impObj)),0)
  
  # make sure that non imputed values remain the same between both objects
  impObjVals <- impObj[!is.na(mdataFilt$e_data[,-1])]
  nonimpObjVals <- mdataFilt$e_data[,-1][!is.na(mdataFilt$e_data[,-1])]
  
  expect_equal(impObjVals,nonimpObjVals)
  
  # apply the imp object back into a pmart object ------------------------------
  # apply the imputation to omics object
  mdataImp <- apply_imputation(impObj,mdataFilt)
  
  # expect error if we run imputation on a dataset with no missing values
  expect_error(imputation(mdataImp),
               "There are no missing values in this dataset. Imputation is not necessary.")
  
  # make sure we have the right class of the data
  expect_equal(class(impObj),c('imputeData','data.frame'))
})
