context('Run combat batch correction')

test_that('bc_combat returns the correct data frame and attributes',{

  # Load the reduced peptide data frames ---------------------------------------
  load(system.file('testdata',
                   'little_pdata.RData',
                   package = 'malbacR'))
  
  # Run as.pepData with agreeable data frames ----------------------------------
  
  # Construct a pepData object with the edata, fdata, and emeta data frames.
  fdata$SampleID <- as.character(fdata$SampleID)
  pdata <- pmartR::as.pepData(e_data = edata,
                      f_data = fdata,
                      e_meta = emeta,
                      edata_cname = 'Mass_Tag_ID',
                      fdata_cname = 'SampleID',
                      emeta_cname = 'Protein')
  # add a grouping variable and a batch variable
  pdata$f_data$Grouping <- c(rep("Infection",6),rep("Mock",6))
  pdata$f_data$Batch <- rep(c(1,2),6)
  
  # Error messages if loaded improperly ----------------------------------------
  
  # pdata does not have group designation
  expect_error(bc_combat(pdata),
               "omicsData must have batch_id attribute for batch correction")
  
  # pdata has group_designation but not specifically for batch
  pdata_g <- pmartR::group_designation(pdata,main_effects = "Grouping")
  expect_error(bc_combat(pdata_g),
               "omicsData must have batch_id attribute for batch correction")
  
  # pdata has group_designation and batch but does not have molecule filter with batch applied
  pdata_bg <- pmartR::group_designation(pdata,main_effects = "Grouping",batch_id = "Batch")
  
  # now error out that the data is not normalized
  expect_error(bc_combat(pdata_bg),
               "omicsData must be normalized prior to running ComBat")
  
  # normalize the data
  pdata_bgn <- pmartR::edata_transform(pdata_bg,"log2")
  pdata_bgn <- pmartR::normalize_global(pdata_bgn, subset_fn = "all", norm_fn = "median",
                                        apply_norm = T, backtransform = T)
  # now we error out because we don't have the proper filters
  expect_error(bc_combat(pdata_bgn),
               "omicsData must have molecule filter applied with use_batch = TRUE and min_num = 2")
  # there is the exception that we already have the data we need (like lipid data has fewer missing values

  # need to get code for if molecule_filter not run but the data is as intended
  complete_pep <- pmartR:::complete_mols(e_data = edata,edata_id = "Mass_Tag_ID")
  custfilt <- pmartR::custom_filter(pdata_bgn,e_data_keep = complete_pep)
  pdataComplete <- pmartR::applyFilt(custfilt,pdata_bgn)
  
  # check to make sure molecule filter has not been applied to the data (even molecule filter)
  expect_equal("moleculeFilt" %in% unlist(lapply(attr(pdataComplete,"filters"),function(x){x$type})),
               FALSE)
  
  # run batch correction
  pdataCompleteBC <- bc_combat(pdataComplete)
  expect_equal(dim(pdataCompleteBC$e_data),dim(pdataComplete$e_data))
  expect_equal(dim(pdataCompleteBC$f_data),dim(pdataComplete$f_data))
  expect_equal(dim(pdataCompleteBC$e_meta),dim(pdataComplete$e_meta))
  
  # run batch correction and show the results are as intended ------------------
  # now we run the more likely scenario that we need a molecule filter to filter out the bad eggs
  # apply molecule filter with use_batch = TRUE
  mol_filt <- pmartR::molecule_filter(pdata_bgn, use_batch = TRUE)
  pdata_bgnFilt <- pmartR::applyFilt(mol_filt,pdata_bgn)

  # run bc_combat on the filtered data
  pdata_bc <- bc_combat(pdata_bgnFilt)
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_bc$e_data),
               dim(pdata_bgnFilt$e_data))
  expect_equal(dim(pdata_bc$f_data),
               dim(pdata$f_data))
  expect_equal(dim(pdata_bc$e_meta),
               dim(pdata_bgnFilt$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame.
  expect_equal(attributes(attr(pdata_bc, 'group_DF'))$main_effects,
               'Grouping')
  expect_equal(attributes(attr(pdata_bc, 'group_DF'))$nonsingleton_groups,
               c('Infection', 'Mock'))
  expect_null(attributes(attr(pdata_bc, 'group_DF'))$covariates)
  expect_equal(attributes(attr(pdata_bc, 'group_DF'))$batch_id,
               pdata$f_data[,c(1,4)])
  expect_null(attributes(attr(pdata_bc, 'group_DF'))$time_course)
  
  # Ensurate the remaining attributes have not changed.
  expect_identical(attr(pdata_bgnFilt, 'cnames'),
                   attr(pdata_bc, 'cnames'))
  # all information other than batch info should stay the same
  expect_identical(attributes(pdata_bc)$data_info[1:8],
                   attributes(pdata_bgnFilt)$data_info[1:8])
  # update batch info
  expect_identical(attributes(pdata_bc)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_combat", params = list(use_groups = FALSE)))
  expect_identical(attr(pdata_bgnFilt, 'meta_info'),
                   attr(pdata_bc, 'meta_info'))
  expect_identical(attr(pdata_bgnFilt, 'filters'),
                   attr(pdata_bc, 'filters'))
  
  
  # run batch correction with use_groups ---------------------------------------
  # run bc_combat on the filtered data
  mol_filt <- pmartR::molecule_filter(pdata_bg, use_batch = TRUE)
  pdata_bgnFilt <- pmartR::applyFilt(mol_filt,pdata_bgn)
  pdata_gbc <- bc_combat(pdata_bgnFilt,use_groups = TRUE)

  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(pdata_gbc$e_data),
               dim(pdata_bgnFilt$e_data))
  expect_equal(dim(pdata_gbc$f_data),
               dim(pdata$f_data))
  expect_equal(dim(pdata_gbc$e_meta),
               dim(pdata_bgnFilt$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame.
  expect_equal(attributes(attr(pdata_gbc, 'group_DF'))$main_effects,
               'Grouping')
  expect_equal(attributes(attr(pdata_gbc, 'group_DF'))$nonsingleton_groups,
               c('Infection', 'Mock'))
  expect_null(attributes(attr(pdata_gbc, 'group_DF'))$covariates)
  expect_equal(attributes(attr(pdata_gbc, 'group_DF'))$batch_id,
               pdata$f_data[,c(1,4)])
  expect_null(attributes(attr(pdata_gbc, 'group_DF'))$time_course)
  
  # Ensurate the remaining attributes have not changed.
  expect_identical(attr(pdata_bgnFilt, 'cnames'),
                   attr(pdata_gbc, 'cnames'))
  # all information other than batch info should stay the same
  expect_identical(attributes(pdata_gbc)$data_info[1:8],
                   attributes(pdata_bgnFilt)$data_info[1:8])
  # update batch info
  expect_identical(attributes(pdata_gbc)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "bc_combat", params = list(use_groups = TRUE)))
  expect_identical(attr(pdata_bgnFilt, 'meta_info'),
                   attr(pdata_gbc, 'meta_info'))
  expect_identical(attr(pdata_bgnFilt, 'filters'),
                   attr(pdata_gbc, 'filters'))
  
})
  
