context('Run WaveICA 2.0 Normalization')

test_that('bc_waveica returns the correct data frame and attributes',{
  
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
  # INJECTION CNAME
  # give injection_cname two columns
  expect_error(bc_waveica(omicsData = mdata,injection_cname = c("RunOrderOverall","RunOrderBatch")),
               "Input parameter injection_cname must be of length 1")
  # injection_cname must be a character string
  expect_error(bc_waveica(omicsData = mdata,injection_cname = 1),
               "Input parameter injection_cname must be of class 'character'")
  # injection_cname must be a column in fdata
  expect_error(bc_waveica(omicsData = mdata, injection_cname = "RunningOrder"),
               "Input parameter injection_col must be a column found in f_data of omicsData")
  
  # ALPHA
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",alpha = c(0,1)),
               "Input parameter alpha must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",alpha = "Batch"),
               "Input parameter alpha must be of class 'numeric'")
  # alpha has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",alpha = -1),
               "Input parameter alpha must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",alpha = 5),
               "Input parameter alpha must be between 0 and 1 inclusive")
  
  # CUTOFF
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",cutoff = c(0,1)),
               "Input parameter cutoff must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",cutoff = "Batch"),
               "Input parameter cutoff must be of class 'numeric'")
  # cutoff has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",cutoff = -1),
               "Input parameter cutoff must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",cutoff = 5),
               "Input parameter cutoff must be between 0 and 1 inclusive")
  
  # K
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",alpha = c(0,1)),
               "Input parameter alpha must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",K = "Batch"),
               "Input parameter K must be of class 'numeric'")
  # cutoff has to be greater than 0
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",K = 0),
               "Input parameter K must be greater than 0")
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall",K = -1),
               "Input parameter K must be greater than 0")
  
  # cannot have missing values
  expect_error(bc_waveica(omicsData = mdata,injection_cname = "RunOrderOverall"),
               "WaveICA requires no missing observations. Remove molecules with missing samples")
  
  # Check the dimensions of results --------------------------------------------
  # remove the missing values
  keep <- mdata$e_data[which(complete.cases(mdata$e_data[,-1])),1]
  cfilt <- pmartR::custom_filter(mdata,e_data_keep = keep)
  udn_complete <- pmartR::applyFilt(cfilt,mdata)
  # run wave ica
  udn_wave <- bc_waveica(udn_complete,injection_cname = "RunOrderOverall")

  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(udn_complete$e_data),
               dim(udn_wave$e_data))
  expect_equal(dim(mdata$f_data),
               dim(udn_wave$f_data))
  expect_equal(dim(udn_complete$e_meta),
               dim(udn_wave$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame --------------------
  
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(c("e_data","f_data","e_meta"),attr(udn_wave,"names"))
  # check cnames
  expect_equal(attr(udn_complete,"cnames"),attr(udn_wave,"cnames"))
  # check data info except for batch info
  # check data info
  expect_equal(attributes(mdata)$data_info[1:3],
               attributes(udn_wave)$data_info[1:3])
  expect_equal(attr(udn_wave,"data_info")$num_edata, nrow(udn_wave$e_data))
  expect_equal(attr(udn_wave,"data_info")$num_miss_obs, sum(is.na(udn_wave$e_data)))
  expect_equal(attr(udn_wave,"data_info")$prop_missing, sum(is.na(udn_wave$e_data))/(nrow(udn_wave$e_data)*(ncol(udn_wave$e_data)-1)))
  expect_equal(attr(udn_wave,"data_info")$num_samps, nrow(udn_wave$f_data))
  
  # check check.names
  # expect_equal(attr(udn_complete,"check.names"),attr(udn_wave,"check.names"))
  # check meta_info
  expect_equal(attr(udn_complete,"meta_info"),attr(udn_wave,"meta_info"))
  # check filters
  expect_equal(length(attr(udn_complete,"filters")),length(attr(udn_wave,"filters")))
  # check group_DF
  expect_equal(attr(udn_complete,"group_DF"),attr(udn_wave,"group_DF"))
  
  # batch info should be updated
  expect_identical(attributes(udn_wave)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "waveica", params = list()))
})
