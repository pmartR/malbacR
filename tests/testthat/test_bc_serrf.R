context('Run SERRF Normalization')

test_that('bc_serrf returns the correct data frame and attributes',{
  
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
  
  # omicsData needs to have undergone batch correction
  expect_error(bc_serrf(omicsData = mdata,sampletype_cname = "Sex",test_val = "QC.NIST"),
               "omicsData must have batch_id attribute for batch correction")
  # update the batch correction
  mdata <- pmartR::group_designation(mdata,main_effects = "Age",batch_id = "BatchName")
  
  # what if we enter parameters wrong
  # sampletype_cname
  # cannot have more than 1 element in string
  expect_error(bc_serrf(omicsData = mdata,sampletype_cname = c("Sex","Age"),test_val = "QC.NIST"),
               "Input parameter sampletype_cname must be of length 1")
  # must be numeric
  expect_error(bc_serrf(omicsData = mdata,sampletype_cname = 10,test_val = "QC.NIST"),
               "Input parameter sampletype_cname must be of class 'character'")
  # must be a column in fdata
  expect_error(bc_serrf(omicsData = mdata, sampletype_cname = "Cat",test_val = "QC.NIST"),
               "Input parameter sampletype_cname must be a column found in f_data of omicsData")
  
  # test_val
  # must have one element
  expect_error(bc_serrf(omicsData = mdata, sampletype_cname = "Sex", test_val = c("QC.NIST","QC.Pool")),
               "Input parameter test_val must be of length 1")
  # must not be numeric
  expect_error(bc_serrf(omicsData = mdata,sampletype_cname = "Sex",test_val = 10),
               "Input parameter test_val must be of class 'character'")
  # must be a column in fdata
  expect_error(bc_serrf(omicsData = mdata, sampletype_cname = "Sex",test_val = "Dog"),
               "Input parameter test_val must be a value in sampletype_cname")
  
  # cannot have missing values in at least one sample of train and testing
  expect_error(bc_serrf(omicsData = mdata, sampletype_cname = "Sex",test_val = "QC.NIST"),
               "SERRF requires no missing observations")
  
  # Check the dimensions of results --------------------------------------------
  
  # remove the missing values
  
  molfilt <- pmartR::molecule_filter(mdata)
  mdataFilt <- pmartR::applyFilt(molfilt,mdata,min_num = 2)
  impObj <- imputation(mdataFilt)
  mdataImp <- apply_imputation(impObj,mdataFilt)

  # run serrf
  udn_serrf <- bc_serrf(omicsData = mdataImp, sampletype_cname = "QC", test_val = "QC.NIST")
  
  # how many QC samples are there
  numQC = sum(mdataImp$f_data$QC == "QC.NIST")
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(nrow(mdataImp$e_data),
               nrow(udn_serrf$e_data))
  expect_equal(ncol(mdataImp$e_data),
               ncol(udn_serrf$e_data)+numQC)
  expect_equal(nrow(mdata$f_data),
               nrow(udn_serrf$f_data)+numQC)
  expect_equal(ncol(mdata$f_data),
               ncol(udn_serrf$f_data))
  expect_equal(nrow(mdataImp$e_meta),
               nrow(udn_serrf$e_meta))
  expect_equal(ncol(mdataImp$e_meta),
               ncol(udn_serrf$e_meta))
  
  # Inspecticate the attributes of the bc_serrf data frame --------------------
  
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(c("e_data","f_data","e_meta"),attr(udn_serrf,"names"))
  # check cnames
  expect_equal(attr(mdataImp,"cnames"),attr(udn_serrf,"cnames"))
  # check data info except for batch info
  # check data info
  expect_equal(attributes(mdata)$data_info[1:3],
               attributes(udn_serrf)$data_info[1:3])
  expect_equal(attr(udn_serrf,"data_info")$num_edata, nrow(udn_serrf$e_data))
  expect_equal(attr(udn_serrf,"data_info")$num_miss_obs, sum(is.na(udn_serrf$e_data)))
  expect_equal(attr(udn_serrf,"data_info")$prop_missing, sum(is.na(udn_serrf$e_data))/(nrow(udn_serrf$e_data)*(ncol(udn_serrf$e_data)-1)))
  expect_equal(attr(udn_serrf,"data_info")$num_samps, nrow(udn_serrf$f_data))
  
  # check check.names
  # expect_equal(attr(udnFilt,"check.names"),attr(udn_serrf,"check.names"))
  # check meta_info
  expect_equal(attr(mdataImp,"meta_info"),attr(udn_serrf,"meta_info"))
  # check filters
  expect_equal(length(attr(mdataImp,"filters")),length(attr(udn_serrf,"filters")))
  # check group_DF
  og_group_without_QC = attr(mdataImp,"group_DF") %>% dplyr::filter(Group != "QC.NIST")
  expect_identical(data.frame(og_group_without_QC),
                   data.frame(attr(udn_serrf,"group_DF")))
  
  # batch info should be updated
  expect_identical(attributes(udn_serrf)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "serrf", params = list()))
})
