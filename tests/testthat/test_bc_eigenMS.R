context('Run EigenMS Normalization')

test_that('bc_eigenMS returns the correct data frame and attributes',{
  
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
  
  # Check for any warning errors -----------------------------------------------
  
  # batch id error message
  expect_error(bc_eigenMS(mdata),
               "omicsData must have group_designation ran")
  
  # Check the dimensions of results --------------------------------------------
  
  # run the real version
  mdata <- pmartR::group_designation(mdata,main_effects = "Age",
                                     batch_id = "BatchName")
  
  # we should now get a warning because we have not applied molecule filter with
  # use_groups = TRUE
  expect_error(bc_eigenMS(mdata),
               "omicsData must have molecule filter applied with use_groups = TRUE and min_num = 1")
  
  # run molecule filter
  molfilt <- pmartR::molecule_filter(mdata,use_groups = TRUE)
  mdataFilt <- pmartR::applyFilt(molfilt,mdata,min_num=1)
  udn_scaled <- bc_eigenMS(omicsData = mdataFilt)
  
  # get information about the molecules that we drop along the way
  edat = mdata$e_data[,-1]
  # batch
  replicateCol = attr(mdata,"group_DF")$Group
  replicate = as.factor(replicateCol)
  meta = mdata$e_meta
  # set seed
  set.seed(12345)
  # run bias trends
  bias_trends = suppressWarnings(ProteoMM::eig_norm1(m = edat,
                                    treatment = replicate,
                                    prot.info = meta))
  numMissing <- nrow(mdata$e_data) - nrow(bias_trends$present)
  
  # Investigate the e_data, f_data, and e_meta data frames.
  expect_equal(dim(mdataFilt$e_data),
               dim(udn_scaled$e_data))
  expect_equal(dim(mdata$f_data),
               dim(udn_scaled$f_data))
  expect_equal(dim(mdataFilt$e_meta),
               dim(udn_scaled$e_meta))
  
  # Inspecticate the attributes of the bc_combat data frame.
  # all of these should be the same between old and new data objects
  # check names
  expect_equal(attr(mdata,"names"),attr(udn_scaled,"names"))
  # check cnames
  expect_equal(attr(mdata,"cnames"),attr(udn_scaled,"cnames"))
  # check data info
  expect_equal(attributes(mdata)$data_info[1:2],
               attributes(udn_scaled)$data_info[1:2])
  expect_equal(attributes(udn_scaled)$data_info$norm_info$is_norm,TRUE)
  expect_equal(attr(udn_scaled,"data_info")$num_edata, nrow(udn_scaled$e_data))
  expect_equal(attr(udn_scaled,"data_info")$num_miss_obs, sum(is.na(udn_scaled$e_data)))
  expect_equal(attr(udn_scaled,"data_info")$prop_missing, sum(is.na(udn_scaled$e_data))/(nrow(udn_scaled$e_data)*(ncol(udn_scaled$e_data)-1)))
  expect_equal(attr(udn_scaled,"data_info")$num_samps, nrow(udn_scaled$f_data))
  
  # check check.names
  expect_equal(attr(mdataFilt,"check.names"),attr(udn_scaled,"check.names"))
  # check meta_info
  expect_equal(attr(mdataFilt,"meta_info"),attr(udn_scaled,"meta_info"))
  # check filters
  expect_equal(attr(mdataFilt,"filters"),attr(udn_scaled,"filters"))
  # check group_DF
  expect_equal(attr(mdataFilt,"group_DF"),attr(udn_scaled,"group_DF"))
  
  # batch info should be updated
  expect_identical(attributes(udn_scaled)$data_info$batch_info,
                   list(is_bc = TRUE,bc_method = "eigenMS", params = list()))
})

