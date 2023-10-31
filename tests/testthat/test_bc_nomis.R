context('Run NOMIS Normalization')

test_that('bc_nomis returns the correct data frame and attributes',{
  
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
  
  # We have not added batch_id group designation in the data
  expect_error(bc_nomis(omicsData = mdata,is_cname = "IS",is_val = "IS"),
               "omicsData must have batch_id")
  # so we add batch information
  udn_batched <- pmartR::group_designation(mdata,main_effects = "Age",batch_id = "BatchName")
  
  # what if we enter parameters wrong
  # give is_cname two columns
  expect_error(bc_nomis(omicsData = udn_batched,is_cname = c("IS","QC"),is_val = "IS"),
               "Input parameter is_cname must be of length 1")
  # give is_val two values
  expect_error(bc_nomis(omicsData = udn_batched,is_cname = "IS",is_val = c("IS","NotIS")),
               "Input parameter is_val must be of length 1")
  # is_cname is not a column in fdata
  expect_error(bc_nomis(omicsData = udn_batched,is_cname = "IntStand",is_val = "IS"),
               "Input parameter is_cname must be a column found in e_meta of omicsData.")
  # is_val not found in the is_cname column
  expect_error(bc_nomis(omicsData = udn_batched,is_cname = "IS",is_val = "IntStand"),
               "Input parameter is_val must be a value found in the is_cname column")
  # pick is_cname to not be a character
  expect_error(bc_nomis(omicsData = udn_batched,is_cname = 3,is_val = "IntStand"),
               "Input parameter is_cname must be of class 'character'.")
  # if we have missing data NOMIS won't run
  expect_error(bc_nomis(omicsData = udn_batched,is_cname = "IS",is_val = "IS"),
               "NOMIS requires no missing observations.")
  
  # Check the dimensions of results --------------------------------------------
  # remove the missing values
  keep <- udn_batched$e_data[which(complete.cases(udn_batched$e_data[,-1])),1]
  cfilt <- pmartR::custom_filter(udn_batched,e_data_keep = keep)
  udn_complete <- pmartR::applyFilt(cfilt,udn_batched)
  # run NOMIS
  udn_nomis <- bc_nomis(udn_complete,"IS","IS")
  
  # find the number of samples that are internal standards
  numIS = sum(udn_complete$e_meta$IS == "IS")

  # find the number of molecules with missing data
  keep <- udn_complete$e_data[which(complete.cases(udn_batched$e_data[,-1])),1]
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(ncol(udn_complete$e_data),
               ncol(udn_nomis$e_data))
  expect_equal(nrow(udn_complete$e_data),
               nrow(udn_nomis$e_data) + numIS)
  expect_equal(dim(udn_batched$f_data),
               dim(udn_nomis$f_data))
  expect_equal(ncol(udn_complete$e_meta),
               ncol(udn_nomis$e_meta))
  expect_equal(nrow(udn_complete$e_meta),
               nrow(udn_nomis$e_meta) + numIS)
  
  # Inspecticate the attributes of the bc_combat data frame --------------------
  
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(c("e_data","f_data","e_meta"),attr(udn_nomis,"names"))
  # check cnames
  expect_equal(attr(udn_complete,"cnames"),attr(udn_nomis,"cnames"))
  # check data info except for batch info
  # check data info
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(udn_nomis)$data_info[1:2])
  expect_equal(attributes(udn_nomis)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_nomis,"data_info")$num_edata, nrow(udn_nomis$e_data))
  expect_equal(attr(udn_nomis,"data_info")$num_miss_obs, sum(is.na(udn_nomis$e_data)))
  expect_equal(attr(udn_nomis,"data_info")$prop_missing, sum(is.na(udn_nomis$e_data))/(nrow(udn_nomis$e_data)*(ncol(udn_nomis$e_data)-1)))
  expect_equal(attr(udn_nomis,"data_info")$num_samps, nrow(udn_nomis$f_data))
  
  # check check.names
  # expect_equal(attr(udn_complete,"check.names"),attr(udn_nomis,"check.names"))
  # check meta_info
  expect_equal(attr(udn_complete,"meta_info"),attr(udn_nomis,"meta_info"))
  # check filters
  expect_equal(length(attr(udn_complete,"filters")),length(attr(udn_nomis,"filters")))
  # check group_DF
  expect_equal(attr(udn_complete,"group_DF"),attr(udn_nomis,"group_DF"))
  
  # batch info should be updated
  expect_identical(attributes(udn_nomis)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "nomis", params = list()))
})
