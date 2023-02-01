context('Run filtering for TIGER batch correction')

test_that('tiger_filter returns the correct custom filter and attributes',{
  
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
  
  # must be a pmart object
  expect_error(tiger_filter(omicsData = mdata$e_data,sampletype_cname = "QC",test_val = "QC.NIST"),
               "omicsData must be of class")
  
  # what if we enter parameters wrong
  # sampletype_cname
  # cannot have more than 1 element in string
  expect_error(tiger_filter(omicsData = mdata,sampletype_cname = c("Sex","Age"),test_val = "QC.NIST"),
               "Input parameter sampletype_cname must be of length 1")
  # must be numeric
  expect_error(tiger_filter(omicsData = mdata,sampletype_cname = 10,test_val = "QC.NIST"),
               "Input parameter sampletype_cname must be of class 'character'")
  # must be a column in fdata
  expect_error(tiger_filter(omicsData = mdata, sampletype_cname = "Cat",test_val = "QC.NIST"),
               "Input parameter sampletype_cname must be a column found in f_data of omicsData")
  
  # test_val
  # must have one element
  expect_error(tiger_filter(omicsData = mdata, sampletype_cname = "Sex", test_val = c("QC.NIST","QC.Pool")),
               "Input parameter test_val must be of length 1")
  # must be numeric
  expect_error(tiger_filter(omicsData = mdata,sampletype_cname = "Sex",test_val = 10),
               "Input parameter test_val must be of class 'character'")
  # must be a column in fdata
  expect_error(tiger_filter(omicsData = mdata, sampletype_cname = "Sex",test_val = "Dog"),
               "Input parameter test_val must be a value in sampletype_cname")
  
  # run actual version of the function -----------------------------------------
  
  # run tiger filter on the data
  tigerFilt <- tiger_filter(omicsData = mdata, sampletype_cname = "QC", test_val = "QC.NIST")
  
  # compare to manual version
  qc_order <- mdata$f_data$QC
  edata <- mdata$e_data[,-1]
  qc_edata <- edata[,qc_order == "QC.NIST"]
  nonqc_edata <- edata[,qc_order != "QC.NIST"]
  
  # apply the filter so we can actually run bc_tiger()
  mdataFilt <- apply_tigerFilt(tigerFilt,mdata)
  # need to run batch id group designation for bc_tiger
  mdataFilt <- pmartR::group_designation(omicsData = mdataFilt, main_effects = "Sex", batch_id = "BatchNum")
  
  # all things considered - this function should run smoothly then
  mdata_tiger <- bc_tiger(omicsData = mdataFilt,sampletype_cname = "QC", test_val = "QC.NIST")
})
