context('Run pareto scaling batch correction')

test_that('bc_pareto returns the correct data frame and attributes',{
  
  # Load the reduced peptide data frames ---------------------------------------
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Run as.metabData with agreeable data frames ----------------------------------

  # Construct a pemdata object with the edata, fdata, and emeta data frames.
  fdata$SampleID <- as.character(fdata$SampleID)
  mdata <- pmartR::as.metabData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Metabolite',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Metabolite')
  
  # Check the dimensions of results --------------------------------------------
  mdata_scaled <- bc_pareto(mdata)
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(mdata$e_data),
               dim(mdata_scaled$e_data))
  expect_equal(dim(mdata$f_data),
               dim(mdata_scaled$f_data))
  expect_equal(dim(mdata$e_meta),
               dim(mdata_scaled$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame.
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(attr(mdata,"names"),attr(mdata_scaled,"names"))
  # check cnames
  expect_equal(attr(mdata,"cnames"),attr(mdata_scaled,"cnames"))
  # check data info
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(mdata_scaled)$data_info[1:2])
  expect_equal(attributes(mdata_scaled)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(mdata_scaled,"data_info")$num_edata, nrow(mdata_scaled$e_data))
  expect_equal(attr(mdata_scaled,"data_info")$num_miss_obs, sum(is.na(mdata_scaled$e_data)))
  expect_equal(attr(mdata_scaled,"data_info")$prop_missing, sum(is.na(mdata_scaled$e_data))/(nrow(mdata_scaled$e_data)*(ncol(mdata_scaled$e_data)-1)))
  expect_equal(attr(mdata_scaled,"data_info")$num_samps, nrow(mdata_scaled$f_data))
  
  # check check.names
  expect_equal(attr(mdata,"check.names"),attr(mdata_scaled,"check.names"))
  # check meta_info
  expect_equal(attr(mdata,"meta_info"),attr(mdata_scaled,"meta_info"))
  # check filters
  expect_equal(attr(mdata,"filters"),attr(mdata_scaled,"filters"))
  # check group_DF
  expect_equal(attr(mdata,"group_DF"),attr(mdata_scaled,"group_DF"))
  
  # batch info should be updated
  expect_identical(attributes(mdata_scaled)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_pareto", params = list()))
})

