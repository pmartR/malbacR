context('Run range scaling batch correction')

test_that('bc_pareto returns the correct data frame and attributes',{
  
  # Load the reduced peptide data frames ---------------------------------------
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Run as.metabData with agreeable data frames ----------------------------------
  
  # Construct a metabData object with the edata, fdata, and emeta data frames.
  fdata$SampleID <- as.character(fdata$SampleID)
  mdata <- pmartR::as.metabData(e_data = edata,
                              f_data = fdata,
                              e_meta = emeta,
                              edata_cname = 'Metabolite',
                              fdata_cname = 'SampleID',
                              emeta_cname = 'Metabolite',
                              data_scale = 'log2')
  
  # Warnings -------------------------------------------------------------------
  # we need at each biomolecule present for at least 2 samples
  expect_error(bc_range(mdata),
               "Range Scaling requires that each biomolecule be present in at least 2 samples")

  # remove the missing values
  keep <- mdata$e_data[which(complete.cases(mdata$e_data[,-1])),1]
  cfilt <- pmartR::custom_filter(mdata,e_data_keep = keep)
  mdata_complete <- pmartR::applyFilt(cfilt,mdata)
  
  # what if data was in abundance
  mdata_abundance <- pmartR::edata_transform(mdata_complete,"abundance")
  expect_error(bc_range(mdata_abundance),
               "Range Scaling must be ran with log2 abundance values")
  
  # run with log2 data
  mdata_scaled <- bc_range(mdata_complete)
  
  # Check the dimensions of results --------------------------------------------
  mdata_scaled <- bc_range(mdata_complete)
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(mdata_complete$e_data),
               dim(mdata_scaled$e_data))
  expect_equal(dim(mdata_complete$f_data),
               dim(mdata_scaled$f_data))
  expect_equal(dim(mdata_complete$e_meta),
               dim(mdata_scaled$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame.
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(attr(mdata_complete,"names"),attr(mdata_scaled,"names"))
  # check cnames
  expect_equal(attr(mdata_complete,"cnames"),attr(mdata_scaled,"cnames"))
  # check data info
  # check data info
  expect_equal(attributes(mdata_complete)$data_info[1:2],
               attributes(mdata_scaled)$data_info[1:2])
  expect_equal(attributes(mdata_scaled)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(mdata_scaled,"data_info")$num_edata, nrow(mdata_scaled$e_data))
  expect_equal(attr(mdata_scaled,"data_info")$num_miss_obs, sum(is.na(mdata_scaled$e_data)))
  expect_equal(attr(mdata_scaled,"data_info")$prop_missing, sum(is.na(mdata_scaled$e_data))/(nrow(mdata_scaled$e_data)*(ncol(mdata_scaled$e_data)-1)))
  expect_equal(attr(mdata_scaled,"data_info")$num_samps, nrow(mdata_scaled$f_data))
  
  # check check.names
  expect_equal(attr(mdata_complete,"check.names"),attr(mdata_scaled,"check.names"))
  # check meta_info
  expect_equal(attr(mdata_complete,"meta_info"),attr(mdata_scaled,"meta_info"))
  # check filters
  expect_equal(attr(mdata_complete,"filters"),attr(mdata_scaled,"filters"))
  # check group_DF
  expect_equal(attr(mdata_complete,"group_DF"),attr(mdata_scaled,"group_DF"))
  
  # batch info should be updated
  expect_identical(attributes(mdata_scaled)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_range", params = list()))
})

