context('Run TIGER Normalization')

test_that('bc_tiger returns the correct data frame and attributes',{
  
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
  
  # omicsData needs to have undergone batch id group designation
  expect_error(bc_tiger(omicsData = mdata,sampletype_cname = "Sex",test_val = "QC.NIST",
                        group_cname = "Age"),
               "omicsData must have batch_id attribute for batch correction")
  # update the batch correction
  mdata <- pmartR::group_designation(mdata,main_effects = "Age",batch_id = "BatchName")
  
  # what if we enter parameters wrong
  # sampletype_cname
  # cannot have more than 1 element in string
  expect_error(bc_tiger(omicsData = mdata,sampletype_cname = c("Sex","Age"),test_val = "QC.NIST",
                        group_cname = "Age"),
               "Input parameter sampletype_cname must be of length 1")
  # must be numeric
  expect_error(bc_tiger(omicsData = mdata,sampletype_cname = 10,test_val = "QC.NIST",
                        group_cname = "Age"),
               "Input parameter sampletype_cname must be of class 'character'")
  # must be a column in fdata
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Cat",test_val = "QC.NIST",
                        group_cname = "Age"),
               "Input parameter sampletype_cname must be a column found in f_data of omicsData")
  
  # test_val
  # must have one element
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex", test_val = c("QC.NIST","QC.Pool"),
                        group_cname = "Age"),
               "Input parameter test_val must be of length 1")
  # must be numeric
  expect_error(bc_tiger(omicsData = mdata,sampletype_cname = "Sex",test_val = 10,
                        group_cname = "Age"),
               "Input parameter test_val must be of class 'character'")
  # must be a column in fdata
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex",test_val = "Dog",
                        group_cname = "Age"),
               "Input parameter test_val must be a value in sampletype_cname")
  
  # position_cname
  # must have one element
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex", test_val = "QC.NIST",
                        position_cname = c("Sex","Age"),
                        group_cname = "Age"),
               "Input parameter position_cname must be of length 1")
  # must be numeric
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex", test_val = "QC.NIST",
                        position_cname = 10,
                        group_cname = "Age"),
               "Input parameter position_cname must be of class 'character'")
  # must be a column in sampletype_cname column
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex",test_val = "QC.NIST",
                        position_cname = "Well_Pos",
                        group_cname = "Age"),
               "Input parameter position_cname must be a column found in f_data of omicsData")
  
  # injection_cname
  # must have one element
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex", test_val = "QC.NIST",
                        injection_cname = c("RunOrderOverall","RunOrderBatch"),
                        group_cname = "Age"),
               "Input parameter injection_cname must be of length 1")
  # must be numeric
  expect_error(bc_tiger(omicsData = mdata,sampletype_cname = "Sex",test_val = "QC.NIST",
                        injection_cname = 10,
                        group_cname = "Age"),
               "Input parameter injection_cname must be of class 'character'")
  # must be a column in fdata
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex",test_val = "QC.NIST",
                        injection_cname = "RunningOrder",
                        group_cname = "Age"),
               "Input parameter injection_cname must be a column found in f_data of omicsData")
  
  # group cname must be 1 or 2 length
  # must have one element
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex", test_val = "QC.NIST",
                        injection_cname = "RunOrderOverall",
                        group_cname = c("Age","Sex","Health")),
               "Input parameter group_cname must be of length 1 or 2")
  # must be numeric
  expect_error(bc_tiger(omicsData = mdata,sampletype_cname = "Sex",test_val = "QC.NIST",
                        injection_cname = "RunOrderOverall",
                        group_cname = 10),
               "Input parameter group_cname must be of class 'character'")
  # must be a column in fdata
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex",test_val = "QC.NIST",
                        injection_cname = "RunningOrder",
                        group_cname = "Health"),
               "Input parameter group_cname must be a column found in f_data of omicsData")
  
  # cannot have missing values in at least one sample of train and testing
  expect_error(bc_tiger(omicsData = mdata, sampletype_cname = "Sex",test_val = "QC.NIST",
                        injection_cname = "RunOrderOverall",
                        group_cname = "Age"),
               "At least one sample associated with the test_val and at least one sample not associated")
  
  # Check the dimensions of results --------------------------------------------
  
  # remove the missing values
  molfilt <- pmartR::molecule_filter(mdata)
  mdata <- pmartR::applyFilt(molfilt,mdata)
  mdataImpObj <- imputation(mdata)
  mdataImp <- apply_imputation(mdataImpObj,mdata)

  # run tiger
  udn_tiger <- bc_tiger(omicsData = mdataImp, sampletype_cname = "QC", test_val = "QC.NIST",
                        injection_cname = "RunOrderOverall",group_cname = "Age")
  
  # how many QC samples are there
  numQC = sum(mdataImp$f_data$QC == "QC.NIST")
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(nrow(mdataImp$e_data),
               nrow(udn_tiger$e_data))
  expect_equal(ncol(mdataImp$e_data),
               ncol(udn_tiger$e_data)+numQC)
  expect_equal(nrow(mdata$f_data),
               nrow(udn_tiger$f_data)+numQC)
  expect_equal(ncol(mdata$f_data),
               ncol(udn_tiger$f_data))
  expect_equal(nrow(mdataImp$e_meta),
               nrow(udn_tiger$e_meta))
  expect_equal(ncol(mdataImp$e_meta),
               ncol(udn_tiger$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame --------------------
  
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(c("e_data","f_data","e_meta"),attr(udn_tiger,"names"))
  # check cnames
  expect_equal(attr(mdataImp,"cnames"),attr(udn_tiger,"cnames"))
  # check data info except for batch info
  # check data info
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(udn_tiger)$data_info[1:2])
  expect_equal(attributes(udn_tiger)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_tiger,"data_info")$num_edata, nrow(udn_tiger$e_data))
  expect_equal(attr(udn_tiger,"data_info")$num_miss_obs, sum(is.na(udn_tiger$e_data)))
  expect_equal(attr(udn_tiger,"data_info")$prop_missing, sum(is.na(udn_tiger$e_data))/(nrow(udn_tiger$e_data)*(ncol(udn_tiger$e_data)-1)))
  expect_equal(attr(udn_tiger,"data_info")$num_samps, nrow(udn_tiger$f_data))
  
  # check check.names
  # expect_equal(attr(mdataImp,"check.names"),attr(udn_tiger,"check.names"))
  # check meta_info
  expect_equal(attr(mdataImp,"meta_info"),attr(udn_tiger,"meta_info"))
  # check filters
  expect_equal(length(attr(mdataImp,"filters")),length(attr(udn_tiger,"filters")))
  # check group_DF
  mdataImp_groupDF_filt = attr(mdataImp,"group_DF") %>% dplyr::filter(!stringr::str_detect(SampleID,"QC.NIST"))
  expect_equal(attributes(mdataImp_groupDF_filt)[1:4],attributes(attr(udn_tiger,"group_DF"))[1:4])
  expect_equal(attributes(mdataImp_groupDF_filt)$nonsingleton_groups[1:5],attributes(attr(udn_tiger,"group_DF"))$nonsingleton_groups)
  batchImp_df <- attributes(attr(mdataImp,"group_DF"))$batch_id %>% dplyr::filter(!stringr::str_detect(SampleID,"QC.NIST"))
  expect_equal(batchImp_df,attributes(attr(udn_tiger,"group_DF"))$batch_id)
  
  # batch info should be updated
  expect_identical(attributes(udn_tiger)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_tiger", params = list(sampletype_cname = "QC",
                                                                           test_val = "QC.NIST",
                                                                           group_cname = "Age",
                                                                           position_cname = NULL,
                                                                           injection_cname = "RunOrderOverall")))
})
