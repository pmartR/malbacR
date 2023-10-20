#' TIGER batch correction
#' 
#' This function returns a pmart object that has been undergone TIGER
#'  batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param sampletype_cname character string giving the name of the column in omicsData$f_data
#'  that contains the sample type information (such as quality control samples)
#' @param test_val character string giving the name of the value within the column sampletype_cname to be used
#'  as the testing value for TIGER
#' @param position_cname character string giving the name of the column in omicsData$f_data that
#'  contains the well position, default is NULL
#' @param injection_cname character string giving the name of the column in omicsData$f_data that
#'  contains the injection order, default is NULL
#'  
#' @return Object of same class as omicsData that has been undergone
#'   TIGER normalization
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
#'                                apply_norm = TRUE,backtransform = TRUE)
#' tigerFilt <- tiger_filter(pmart_amide,sampletype_cname = "group",test_val = "QC")
#' pmart_amideFilt <- apply_tigerFilt(tigerFilt,pmart_amide)
#' amide_tiger <- bc_tiger(omicsData = pmart_amideFilt,sampletype_cname = "group",
#'                         test_val = "QC",injection_cname = "Injection_order")
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_tiger <- function(omicsData, sampletype_cname,test_val,group_cname,position_cname = NULL,injection_cname = NULL){
  
  # run through checks ---------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that omicsData has batch id information
  if (is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    stop (paste("omicsData must have batch_id attribute for batch correction",
                sep = ' '))
  }
  
  # sampletype_cname - type of each sample
  if (class(sampletype_cname) != "character") {
    stop("Input parameter sampletype_cname must be of class 'character'.")
  }
  
  if (length(sampletype_cname) != 1) {
    stop("Input parameter sampletype_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$f_data) == sampletype_cname)) {
    stop("Input parameter sampletype_cname must be a column found in f_data of omicsData.")
  }
  
  # test_val - the value that will be used to determine what is the testing and what is the training
  if (class(test_val) != "character") {
    stop("Input parameter test_val must be of class 'character'.")
  }
  if (length(test_val) != 1) {
    stop("Input parameter test_val must be of length 1 (e.g. vector containing a single element)")
  }
  if (!test_val %in% omicsData$f_data[,sampletype_cname]){
    stop("Input parameter test_val must be a value in sampletype_cname column in omicsData$f_data")
  }
  
  # group_cname - type of each sample
  if (class(group_cname) != "character") {
    stop("Input parameter group_cname must be of class 'character'.")
  }
  
  if (length(group_cname) > 2) {
    stop("Input parameter group_cname must be of length 1 or 2 (e.g. vector containing a one or two elements")
  }
  
  if (!any(names(omicsData$f_data) == group_cname)) {
    stop("Input parameter group_cname must be a column found in f_data of omicsData.")
  }
  
  # injection_cname - injection order or temporal information of each sample
  if(!is.null(injection_cname)){
    if (class(injection_cname) != "character") {
      stop("Input parameter injection_cname must be of class 'character'.")
    }
    
    if (length(injection_cname) != 1) {
      stop("Input parameter injection_cname must be of length 1 (e.g. vector containing a single element")
    }
    
    if (!any(names(omicsData$f_data) == injection_cname)) {
      stop("Input parameter injection_cname must be a column found in f_data of omicsData.")
    }
  }
  
  # position cname - well position of each sample
  if(!is.null(position_cname)){
    if (class(position_cname) != "character") {
      stop("Input parameter position_cname must be of class 'character'.")
    }
    
    if (length(position_cname) != 1) {
      stop("Input parameter position_cname must be of length 1 (e.g. vector containing a single element")
    }
    
    if (!any(names(omicsData$f_data) == position_cname)) {
      stop("Input parameter position_cname must be a column found in f_data of omicsData.")
    }
    
    # tiger only allows position information to be numeric
    if (!is.numeric(omicsData$f_data[,position_cname])){
      stop("The input parameter position_cname must be a column in f_data with numeric values")
    }
  }
  
  # useful information
  edata_cname <- pmartR::get_edata_cname(omicsData)
  fdata_cname <- pmartR::get_fdata_cname(omicsData)
  edata_cnameCol <- which(colnames(omicsData$e_data) == edata_cname)
  fdata_cnameCol <- which(colnames(omicsData$f_data) == fdata_cname)
  
  # TIGER does like missing values
  #which column has the edata cname
  # if (sum(is.na(omicsData$e_data[,-edata_cnameCol])) != 0) {
  #  stop ("TIGER requires no missing observations. Remove molecules with missing samples.")
  # }
  
  # run the TIGER calculations -------------------------------------------------
  
  # make sure we are aligning fdata and edata correctly
  fdata_sampleID_ordering <- omicsData$f_data[fdata_cnameCol]

  # pull out the edata information
  edata_og <- omicsData$e_data
  row.names(edata_og) <- edata_og[,edata_cnameCol]
  # create the transpose of the edata
  row.names(edata_og) <- edata_og[,1]
  edata = t(edata_og[,-1])
  edata <- data.frame(edata)
  edata <- tibble::rownames_to_column(edata,var = fdata_cname)
  proper_order <- match(edata[,1],omicsData$f_data[,fdata_cnameCol])
  edata <- edata[proper_order,]
  
  # simplify the fdata into only the columns that we need
  # important fdata columns
  fdata_info <- c(fdata_cname,sampletype_cname,injection_cname,position_cname)
  fdata <- omicsData$f_data %>%
    dplyr::select(dplyr::all_of(fdata_info))
  
  # add on batch information
  fdata$Batch <- attributes(attr(omicsData,"group_DF"))$batch_id[,2]
  fdata_info <- c(fdata_info,"Batch")

  # join the fdata with the edata for tiger operatios
  tigerSetup <- fdata %>%
    dplyr::left_join(edata, by = fdata_cname)
  
  # set up the train and test samples
  train_samples <- tigerSetup[tigerSetup[,sampletype_cname] == test_val,]
  test_samples  <- tigerSetup[tigerSetup[,sampletype_cname] != test_val,]
  
  if(sum(!unique(test_samples$Batch) %in% unique(train_samples$Batch)) != 0 |
     sum(!unique(train_samples$Batch) %in% unique(test_samples$Batch)) != 0 ){
    stop("All batches must have at least one sample in which the sample type is the value test_val and at least
         one sample in which the sample type is not test_val")
  }
  
  # find how many samples with NA values
  numNAtrain <- rowSums(is.na(dplyr::select_if(train_samples,is.numeric)))
  numNAtest <- rowSums(is.na(dplyr::select_if(test_samples,is.numeric)))
  # error out if we don't have enough observations in either the test or training data
  if(min(numNAtrain) != 0 | min(numNAtest) != 0){
    stop("At least one sample associated with the test_val and at least one sample not associated
         with the test_val must have no missing observations.")
  }
  
  # find out which columns are the ones that we will be keeping
  remaining_samps <- which(tigerSetup[,sampletype_cname] != test_val)
  # set seed
  set.seed(1)
  # include injection order and well position into feature set:
  udnTIGER <- run_TIGER_internal(test_samples = test_samples,
                        train_samples = train_samples,
                        col_sampleID  = fdata_cname,     # input column name
                        col_sampleType = sampletype_cname,  # input column name
                        col_batchID = "Batch",        # input column name
                        col_order = injection_cname,   # input column name
                        col_position = position_cname,  # input column name
                        parallel.cores = parallel::detectCores()-1)
  fdata_info2 <- c(sampletype_cname,injection_cname,position_cname,"Batch")

  rownames(udnTIGER) <- udnTIGER[,fdata_cname]
  edata_tiger <- udnTIGER %>%
    dplyr::select(-dplyr::all_of(fdata_info2))
  
  edata_tiger = t(edata_tiger[,-1])
  edata_tiger <- data.frame(edata_tiger)
  edata_tiger <- tibble::rownames_to_column(edata_tiger,var = edata_cname)
  edata_tiger[,1] <- omicsData$e_data[,edata_cnameCol]
  # find the proper order of the edata columns
  # first we find the qc samples
  qc_names <- train_samples[,fdata_cname]
  # remove them from original order
  proper_order <- colnames(omicsData$e_data)[!colnames(omicsData$e_data) %in% qc_names]
  # rearrange the columns
  edata_tiger <- edata_tiger[,match(proper_order,colnames(edata_tiger))]
  # sometime TIGER can cause a value to be infinite -> make it NA
  edata_tiger[sapply(edata_tiger, is.infinite)] <- NA
  
  # create fdata information
  samp_names <- udnTIGER[,fdata_cname]
  fdata_tiger <- omicsData$f_data[omicsData$f_data[,fdata_cname] %in% samp_names,]

  # move back into pmart object ------------------------------------------------
  
  # assume emeta is NULL unless otherwise stated
  emeta_cname = NULL
  emet = NULL
  if(!is.null(omicsData$e_meta)){
    emeta_cname = pmartR::get_emeta_cname(omicsData)
    emet = omicsData$e_meta
  }
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edata_tiger,
                                          edata_cname = edata_cname,
                                          f_data = fdata_tiger,
                                          fdata_cname = fdata_cname,
                                          e_meta = emet,
                                          emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edata_tiger,
                                  edata_cname = edata_cname,
                                  f_data = fdata_tiger,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edata_tiger,
                                  edata_cname = edata_cname,
                                  f_data = fdata_tiger,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edata_tiger,
                                    edata_cname = edata_cname,
                                    f_data = fdata_tiger,
                                    fdata_cname = fdata_cname,
                                    e_meta = emet,
                                    emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edata_tiger,
                                    edata_cname = edata_cname,
                                    f_data = fdata_tiger,
                                    fdata_cname = fdata_cname,
                                    e_meta = emet,
                                    emeta_cname = emeta_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edata_tiger,
                                  edata_cname = edata_cname,
                                  f_data = fdata_tiger,
                                  fdata_cname = fdata_cname,
                                  e_meta = emet,
                                  emeta_cname = emeta_cname)
  }
  
  # Update the data_info attribute.
  attr(pmartObj, 'data_info') <- pmartR:::set_data_info(
    e_data = pmartObj$e_data,
    edata_cname = pmartR::get_edata_cname(omicsData),
    data_scale_orig = pmartR::get_data_scale_orig(omicsData),
    data_scale = pmartR::get_data_scale(omicsData),
    data_types = pmartR::get_data_info(omicsData)$data_types,
    norm_info = pmartR::get_data_info(omicsData)$norm_info,
    is_normalized = pmartR::get_data_info(omicsData)$norm_info$is_normalized,
    batch_info = pmartR::get_data_info(omicsData)$batch_info,
    is_bc = pmartR::get_data_info(omicsData)$batch_info$is_bc
  )
  
  # check group designation
  # since we are removing samples if keep_qc != TRUE
  if(!is.null(attributes(attr(omicsData,"group_DF"))$batch_id)){
    batch_id_col = which(colnames(attributes(attr(omicsData,"group_DF"))$batch_id) != fdata_cname)
    batch_id_name = colnames(attributes(attr(omicsData,"group_DF"))$batch_id)[batch_id_col]
    pmartObj <- pmartR::group_designation(pmartObj,main_effects = group_cname,
                                          batch_id = batch_id_name)
  } else {
    pmartObj <- pmartR::group_designation(pmartObj,main_effects = group_cname)
  }
  
  # Update the data_info attribute.
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "tiger",
    params = list()
  )
  
  # Update the meta_info attribute.
  attr(pmartObj, 'meta_info') <- pmartR:::set_meta_info(
    e_meta = omicsData$e_meta,
    emeta_cname = pmartR::get_emeta_cname(omicsData)
  )
  
  # Update the filters attribute (but keep all the filtering)
  attr(pmartObj, 'filters') <- attr(omicsData,'filters')
  
  # return pmart object
  return(pmartObj)
}
