context('Run apply_imputation')

test_that('apply_imputation impute the data as expected',{
  
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
  
  # Get the imputation object --------------------------------------------------
  # so we remove those values
  molfilt <- pmartR::molecule_filter(mdata)
  mdataFilt <- pmartR::applyFilt(molfilt,mdata)
  impObj <- imputation(mdataFilt)
  
  # apply the imp object back into a pmart object ------------------------------
  # potential warnings ---------------------------------------------------------
  expect_error(apply_imputation("dog",mdataFilt),
               "The input to the imputeData argument must be an imputeData object")
  expect_error(apply_imputation(impObj,mdataFilt$e_data),
               "omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData'")
  
  # here is the right one
  mdataImp <- apply_imputation(impObj,mdataFilt)
  
  # results checking -----------------------------------------------------------
  
  # the dimensions should be the same
  expect_equal(dim(mdataImp$e_data),dim(mdataFilt$e_data))
  expect_equal(dim(mdataImp$f_data),dim(mdataFilt$f_data))
  expect_equal(dim(mdataImp$e_meta),dim(mdataFilt$e_meta))
  
  # in fact the pmart object's fdata and emeta should be identical
  expect_identical(mdataImp$f_data,mdataFilt$f_data)
  expect_identical(mdataImp$e_meta,mdataFilt$e_meta)
  
  # edata should be the same except with regards to missing values
  # find which are the missing values
  valNA <- which(is.na(as.numeric(unlist(mdataFilt$e_data[,-1]))))
  # vector of observed values in both
  vectorFilt <- as.numeric(unlist(mdataFilt$e_data[,-1]))
  vectorFilt <- vectorFilt[-valNA]
  vectorImp <- as.numeric(unlist(mdataImp$e_data[,-1]))
  vectorImp <- vectorImp[-valNA]
  expect_identical(vectorFilt,vectorImp)
  
  # we should have pmart object return
  expect_equal(class(mdataImp),'metabData')
})
