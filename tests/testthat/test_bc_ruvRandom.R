context('Run RUV-Random Normalization')

test_that('bc_ruvRandom returns the correct data frame and attributes',{
  
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
  # what if we enter parameters wrong
  # give nc_cname two columns
  expect_error(bc_ruvRandom(omicsData = mdata,nc_cname = c("IS","QC"),nc_val = "IS"),
               "Input parameter nc_cname must be of length 1")
  # give nc_val two values
  expect_error(bc_ruvRandom(omicsData = mdata,nc_cname = "IS",nc_val = c("IS","NotIS")),
               "Input parameter nc_val must be of length 1")
  # nc_cname is not a column in fdata
  expect_error(bc_ruvRandom(omicsData = mdata,nc_cname = "IntStand",nc_val = "IS"),
               "Input parameter nc_cname must be a column found in e_meta of omicsData.")
  # nc_val not found in the nc_cname column
  expect_error(bc_ruvRandom(omicsData = mdata,nc_cname = "IS",nc_val = "IntStand"),
               "Input parameter nc_val must be a value found in the nc_cname column")
  # pick nc_cname to not be a character
  expect_error(bc_ruvRandom(omicsData = mdata,nc_cname = 3,nc_val = "IntStand"),
               "Input parameter nc_cname must be of class 'character'.")
  # if we have missing data ruv-random won't run
  expect_error(bc_ruvRandom(omicsData = mdata,nc_cname = "IS",nc_val = "IS"),
               "RUV-random requires no missing observations.")
  
  # Check the dimensions of results --------------------------------------------
  # remove the missing values
  keep <- mdata$e_data[which(complete.cases(mdata$e_data[,-1])),1]
  cfilt <- pmartR::custom_filter(mdata,e_data_keep = keep)
  udn_complete <- pmartR::applyFilt(cfilt,mdata)
  # run ruv-random
  udn_ruv <- bc_ruvRandom(udn_complete,"IS","IS")
  numNC = sum(udn_complete$e_meta$IS == "IS")
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(ncol(udn_complete$e_data),
               ncol(udn_ruv$e_data))
  expect_equal(nrow(udn_complete$e_data),
               nrow(udn_ruv$e_data)+numNC)
  expect_equal(dim(mdata$f_data),
               dim(udn_ruv$f_data))
  expect_equal(ncol(udn_complete$e_meta),
               ncol(udn_ruv$e_meta))
  expect_equal(nrow(udn_complete$e_meta),
               nrow(udn_ruv$e_meta)+numNC)
  
  # Inspecticate the attributes of the bc_combat data frame --------------------
  
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(c("e_data","f_data","e_meta"),attr(udn_ruv,"names"))
  # check cnames
  expect_equal(attr(udn_complete,"cnames"),attr(udn_ruv,"cnames"))
  # check data info except for batch info
  # check data info
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(udn_ruv)$data_info[1:2])
  expect_equal(attributes(udn_ruv)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_ruv,"data_info")$num_edata, nrow(udn_ruv$e_data))
  expect_equal(attr(udn_ruv,"data_info")$num_miss_obs, sum(is.na(udn_ruv$e_data)))
  expect_equal(attr(udn_ruv,"data_info")$prop_missing, sum(is.na(udn_ruv$e_data))/(nrow(udn_ruv$e_data)*(ncol(udn_ruv$e_data)-1)))
  expect_equal(attr(udn_ruv,"data_info")$num_samps, nrow(udn_ruv$f_data))
  
  # check check.names
  # expect_equal(attr(udn_complete,"check.names"),attr(udn_ruv,"check.names"))
  # check meta_info
  expect_equal(attr(udn_complete,"meta_info"),attr(udn_ruv,"meta_info"))
  # check filters
  expect_equal(length(attr(udn_complete,"filters")),length(attr(udn_ruv,"filters")))
  # check group_DF
  expect_equal(attr(udn_complete,"group_DF"),attr(udn_ruv,"group_DF"))
  
  # batch info should be updated
  expect_identical(attributes(udn_ruv)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "ruv_random", params = list()))
})
