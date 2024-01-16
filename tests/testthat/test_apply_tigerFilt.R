context('Run applying the filter for TIGER batch correction')

test_that('tiger_filter returns the correct custom filter and attributes',{
  
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
                                emeta_cname = 'Metabolite',
                                data_scale = 'log2')
  
  # run tiger filter
  tigerFilt <- tiger_filter(omicsData = mdata, sampletype_cname = "QC", test_val = "QC.NIST")
  
  # Run through the potential error messages -----------------------------------

  # must be a pmart object
  expect_error(apply_tigerFilt(filter_object = tigerFilt,omicsData = mdata$e_data),
               "omicsData must be of class 'pepData', 'proData', 'metabData', 'lipidData', 'nmrData', or 'seqData'")
  
  # filter object must be of class tigerFilt
  expect_error(apply_tigerFilt(filter_object = tigerFilt$e_data_remove,omicsData = mdata),
               "filter_object must be of class tigerFilt")
  # we could even try a different filter and it would still be wrong
  molfilt <- pmartR::molecule_filter(mdata)
  expect_error(apply_tigerFilt(filter_object = molfilt,omicsData = mdata),
               "filter_object must be of class tigerFilt")
  
  # Run the data and obtain results --------------------------------------------
  mdataFilt <- apply_tigerFilt(filter_object = tigerFilt, omicsData = mdata)
  
  # make sure the dimensions are updated accurately
  numRemove = length(tigerFilt$e_data_remove)
  # edata
  expect_equal(ncol(mdataFilt$e_data),
               ncol(mdata$e_data))
  expect_equal(nrow(mdataFilt$e_data),
               nrow(mdata$e_data) - numRemove)
  # fdata
  expect_equal(mdataFilt$f_data,
               mdata$f_data)
  # emeta
  expect_equal(ncol(mdataFilt$e_meta),
               ncol(mdata$e_meta))
  expect_equal(nrow(mdataFilt$e_meta),
               nrow(mdata$e_meta) - numRemove)
  
  # make sure the names of metabolites that should and should not be in the 
  # new dataset are correct
  expect_equal(mdata$e_data$Metabolite[!mdata$e_data$Metabolite %in% mdataFilt$e_data$Metabolite],
               tigerFilt$e_data_remove)
  
  # all things considered - this function should run smoothly then
  mdataFilt <- pmartR::group_designation(omicsData = mdataFilt, main_effects = "Sex", batch_id = "BatchNum")
  # must convert to abundance values
  mdata_abundance <- pmartR::edata_transform(mdataFilt,"abundance")
  mdata_tiger <- bc_tiger(omicsData = mdata_abundance,sampletype_cname = "QC", test_val = "QC.NIST",group_cname = "Age")
})
