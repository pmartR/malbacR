context('Run WaveICA Normalization')

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
  # WAVICA OG
  # WAVEICA 2.0
  # BATCH CNAME
  # give batch_cname two columns
  expect_error(bc_waveica(omicsData = mdata,batch_cname = c("batch","Batch")),
               "Input parameter batch_cname must be of length 1")
  # batch_cname must be a character string
  expect_error(bc_waveica(omicsData = mdata,batch_cname = 1),
               "Input parameter batch_cname must be of class 'character'")
  # batch_cname must be a column in fdata
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "Batch"),
               "Input parameter batch_cname must be a column found in f_data of omicsData")
  
  # ALPHA
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",alpha = c(0,1)),
               "Input parameter alpha must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",alpha = "Batch"),
               "Input parameter alpha must be of class 'numeric'")
  # alpha has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",alpha = -1),
               "Input parameter alpha must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",alpha = 5),
               "Input parameter alpha must be between 0 and 1 inclusive")
  
  # CUTOFF_BATCH
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_batch = c(0,1)),
               "Input parameter cutoff_batch must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_batch = "Batch"),
               "Input parameter cutoff_batch must be of class 'numeric'")
  # cutoff_batch has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_batch = -1),
               "Input parameter cutoff_batch must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_batch = 5),
               "Input parameter cutoff_batch must be between 0 and 1 inclusive")
  
  # CUTOFF_GROUP
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_group = c(0,1)),
               "Input parameter cutoff_group must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_group = "Batch"),
               "Input parameter cutoff_group must be of class 'numeric'")
  # cutoff_group has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_group = -1),
               "Input parameter cutoff_group must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",cutoff_group = 5),
               "Input parameter cutoff_group must be between 0 and 1 inclusive")
  
  # K
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",alpha = c(0,1)),
               "Input parameter alpha must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",K = "Batch"),
               "Input parameter K must be of class 'numeric'")
  # K has to be greater than 0
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",K = 0),
               "Input parameter K must be greater than 0")
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum",K = -1),
               "Input parameter K must be greater than 0")
  
  # cannot have missing values
  expect_error(bc_waveica(omicsData = mdata,batch_cname = "BatchNum"),
               "WaveICA requires no missing observations. Remove molecules with missing samples")
  
  # Check the dimensions of results --------------------------------------------
  # remove the missing values
  keep <- mdata$e_data[which(complete.cases(mdata$e_data[,-1])),1]
  cfilt <- pmartR::custom_filter(mdata,e_data_keep = keep)
  udn_complete <- pmartR::applyFilt(cfilt,mdata)
  # run wave ica
  udn_wave <- bc_waveica(udn_complete,batch_cname = "BatchNum")
  
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
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(udn_wave)$data_info[1:2])
  expect_equal(attributes(udn_wave)$data_info$norm_info$is_norm,TRUE)
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
                   list(is_bc = TRUE,bc_method = "bc_waveica", params = list(version = "WaveICA",
                                                                             injection_cname = NULL,
                                                                             batch_cname = "BatchNum",
                                                                             cutoff_injection = 0.1,
                                                                             cutoff_batch = 0.05,
                                                                             cutoff_group = 0.05,
                                                                             alpha = 0,
                                                                             K = 20)))
  
  # WAVEICA 2.0
  # INJECTION CNAME
  # give injection_cname two columns
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = c("RunOrderOverall","RunOrderBatch")),
               "Input parameter injection_cname must be of length 1")
  # injection_cname must be a character string
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = 1),
               "Input parameter injection_cname must be of class 'character'")
  # injection_cname must be a column in fdata
  expect_error(bc_waveica(omicsData = mdata, version = "WaveICA2.0", injection_cname = "RunningOrder"),
               "Input parameter injection_cname must be a column found in f_data of omicsData")
  
  # ALPHA
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",alpha = c(0,1)),
               "Input parameter alpha must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",alpha = "Batch"),
               "Input parameter alpha must be of class 'numeric'")
  # alpha has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",alpha = -1),
               "Input parameter alpha must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",alpha = 5),
               "Input parameter alpha must be between 0 and 1 inclusive")
  
  # CUTOFF_INJECTION
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",cutoff_injection = c(0,1)),
               "Input parameter cutoff_injection must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",cutoff_injection = "Batch"),
               "Input parameter cutoff_injection must be of class 'numeric'")
  # cutoff_injection has to be between 0 and 1
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",cutoff_injection = -1),
               "Input parameter cutoff_injection must be between 0 and 1 inclusive")
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",cutoff_injection = 5),
               "Input parameter cutoff_injection must be between 0 and 1 inclusive")
  
  # K
  # cannot have two elements in vector
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",alpha = c(0,1)),
               "Input parameter alpha must be of length 1")
  # must be numeric
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",K = "Batch"),
               "Input parameter K must be of class 'numeric'")
  # K has to be greater than 0
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",K = 0),
               "Input parameter K must be greater than 0")
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall",K = -1),
               "Input parameter K must be greater than 0")
  
  # cannot have missing values
  expect_error(bc_waveica(omicsData = mdata,version = "WaveICA2.0",injection_cname = "RunOrderOverall"),
               "WaveICA requires no missing observations. Remove molecules with missing samples")
  
  # Check the dimensions of results --------------------------------------------
  # remove the missing values
  keep <- mdata$e_data[which(complete.cases(mdata$e_data[,-1])),1]
  cfilt <- pmartR::custom_filter(mdata,e_data_keep = keep)
  udn_complete <- pmartR::applyFilt(cfilt,mdata)
  # run wave ica
  udn_wave <- bc_waveica(udn_complete,version = "WaveICA2.0",injection_cname = "RunOrderOverall")

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
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(udn_wave)$data_info[1:2])
  expect_equal(attributes(udn_wave)$data_info$norm_info$is_norm,TRUE)
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
                   list(is_bc = TRUE,bc_method = "bc_waveica", params = list(version = "WaveICA2.0",
                                                                             injection_cname = "RunOrderOverall",
                                                                             batch_cname = NULL,
                                                                             cutoff_injection = 0.1,
                                                                             cutoff_batch = 0.05,
                                                                             cutoff_group = 0.05,
                                                                             alpha = 0,
                                                                             K = 20)))
})
