#' Run QC-RLSC Normalization
#'
#' Obtain abundance values in a pmart object for QC-RLSC Normalization
#'
#' @param omicsData an omicsData object (metabData, lipidData, pepData, or
#'  proData) created using the pmartR package, where any zero values have
#'  already been replaced with NAs and a log transformation has already been
#'  applied
#' @param block_cname character string giving name of column in omicsData$f_Data
#'  that contains block (or batch) information
#' @param qc_cname character string giving name of column in omicsData$f_data
#'  that contains the factor variable indicating whether sample is QC or not
#' @param qc_val character string giving the value from the qc_cname column that
#'  indicates a QC sample
#' @param order_cname character string giving name of column in omicsData$f_data
#'  that contains the run order
#' @param missing_thresh numeric threshold, between 0 and 1, used for filtering
#'  out biomolecules. See details for more information. A value of 0.5 is
#'  reasonable.
#' @param rsd_thresh numeric threshold used for filtering metabolites. See
#'  details for more information. A value of 0.3 is reasonable.
#' @param backtransform logical value indicated whether or not to backtransform
#'   the data to put the normalized data back on the same scale as the original
#'   data. Defaults to FALSE. If TRUE, the median of the y-values of the loess
#'   curve is added to the normalized value for each biomolecule for each batch
#'
#' @author Damon Leach
#' 
#' @return omicsData object of same class as omicsData_filt, where e_data contains the QC-RLSC normalized values
#' 
#' @examples
#' \dontrun{
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
#'                                apply_norm = TRUE,backtransform = TRUE)
#' amide_qcrlsc <- bc_qcrlsc(omicsData = pmart_amide,block_cname = "batch",
#'                           qc_cname = "group", qc_val = "QC", order_cname = "Injection_order",
#'                           missing_thresh = 0.5, rsd_thresh = 0.3, backtransform  = FALSE)
#'  
#' @export
bc_qcrlsc <- function(omicsData,block_cname,qc_cname,qc_val,
                   order_cname,missing_thresh = 0.5,rsd_thresh = 0.3,
                   backtransform = FALSE,keep_qc = FALSE) {
  
  # run a plethora of checks for this to run -----------------------------------
  if (!class(omicsData) %in% c("metabData", "lipidData", "pepData", "proData")) {
    stop(
      "omicsData must be an S3 object of class 'metabData', 'lipidData', 'pepData', 'proData', or 'nmrData'. See pmartR package for more details."
    )
  }
  
  # check that input parameters that should be character strings are indeed character strings, and each is a vector of length 1#
  if (class(block_cname) != "character") {
    stop("Input parameter block_cname must be of class 'character'.")
  }
  if (class(qc_cname) != "character") {
    stop("Input parameter qc_cname must be of class 'character'.")
  }
  if (class(qc_val) != "character") {
    stop("Input parameter qc_val must be of class 'character'.")
  }
  
  if (length(block_cname) != 1) {
    stop(
      "Input parameter block_cname must be of length 1 (e.g. vector containing a single element)."
    )
  }
  
  # values in block_cname need to be numeric for this to run I believe
  if (length(qc_cname) != 1) {
    stop(
      "Input parameter qc_cname must be of length 1 (e.g. vector containing a single element)."
    )
  }
  if (length(qc_val) != 1) {
    stop("Input parameter qc_val must be of length 1 (e.g. vector containing a single element).")
  }
  
  # check that column names specified exist in f_data #
  # batch information
  if (!any(names(omicsData$f_data) == block_cname)) {
    stop("The f_data component of omicsData must contain a column for 'Batch'.")
  }
  # quality control samples
  if (!any(names(omicsData$f_data) == qc_cname)) {
    stop("The f_data component of omicsData must contain a column for 'QCtype'.")
  }
  # run order 
  if (!any(names(omicsData$f_data) == order_cname)) {
    stop("The f_data component of omicsData must contain a column for 'RunOrder'")
  }
  
  # check that qc_val is present in qc_cname #
  if (!(qc_val %in% omicsData$f_data[, qc_cname])) {
    stop(
      "Input parameter qc_val must be a value found in the qc_cname column from the f_data of omicsData."
    )
  }
  
  # check that values in block_cname are numeric or a factor #
  if (!is.numeric(omicsData$f_data[, block_cname]) & !is.factor(omicsData$f_data[,block_cname])) {
    stop(
      "Values in block_cname column (in f_data component of omicsData) must be either numeric or factor."
    )
  }
  
  # check that the values in run order are numeric #
  if (!is.numeric(omicsData$f_data[, order_cname])) {
    stop(
      "Values in order_cname column (in f_data component of omicsData) must be numeric"
    )
  }
  
  # check that value in keep_qc is logical
  if (!is.logical(keep_qc)) {
    stop("Input parameter keep_qc must be logical (either TRUE or FALSE)")
  }
  
  # check that value in keep_qc is length 1
  if (length(keep_qc) != 1) {
    stop("Input parameter qc_val must be of length 1 (e.g. vector containing a single element).")
  }
  
  # check that backtransform is logical #
  if (!is.logical(backtransform)) {
    stop("Input parameter backtransform must be of class 'logical'.")
  }
  
  # check to see if we have enough non NA values per batch per sample
  
  # important values
  data <- omicsData$e_data
  f_data <- omicsData$f_data
  cnames <- attributes(omicsData)$cnames
  molecule_cname <- cnames$edata_cname
  samp_cname <- cnames$fdata_cname
  
  # maintain the original batch names
  original_batch_names <- omicsData$f_data[,block_cname]
  # if we have a factor Batch information, we need to do some data manipulation
  # save the original names of the Batches
  if(is.factor(omicsData$f_data[,block_cname])){
    # convert to numbers
    omicsData$f_data[,block_cname] <- as.numeric(original_batch_names)
  } else{
    # rank the numbers so they are in chronological order 1 to however many batches we have
    omicsData$f_data[,block_cname] <- match(f_data[,block_cname], 
                                  sort(unique(f_data[,block_cname])))
  }

  # separate out QC samples only from each data source #
  qc_cols = as.character(f_data[which(f_data[, qc_cname] == qc_val), samp_cname])
  qc_dat <-
    data[, which(names(data) %in% c(molecule_cname, qc_cols))]
  names(qc_dat)[1] = molecule_cname
  
  # remove features detected in less than 'missing_thresh'	proportion of samples #
  frac_missing = apply(is.na(qc_dat[, -1]), 1, function(x)
    sum(x) / length(x))
  bad_feats = qc_dat[which(frac_missing > missing_thresh), which(names(qc_dat) ==
                                                                   molecule_cname)]
  
  if (length(bad_feats) > 0) {
    # remove features with too many missing QC samples #
    filt_qc_data1 = qc_dat[-which(qc_dat[, 1] %in% bad_feats), ]
  } else{
    filt_qc_data1 = qc_dat
  }
  
  # split the f_data by batch
  batch_tib <- tibble::tibble(batchDat = split(omicsData$f_data,omicsData$f_data[,block_cname]))
  
  # find if we only have one QC total in a batch
  batch_tib <- batch_tib %>%
    dplyr::mutate(numQC = purrr::map(batchDat,function(bd){
      qc_colNum = which(colnames(bd) == qc_cname)
      qcLength = length(which(bd[,qc_colNum] == qc_val))
    }))
  numQC_perBatch = unlist(batch_tib$numQC)
  if(min(numQC_perBatch,na.rm=T) < 2) {
    stop ("Each batch must have at minimum 2 QC samples (the first and last sample run for each batch must be a QC sample)")
  }
  
  # now we calculate the number of non NA values from filt_qc_data1
  batch_tib <- batch_tib %>%
    dplyr::mutate(
      qcSample = purrr::map(batchDat,function(bd){
        qc_colNum = which(colnames(bd) == qc_cname)
        qcVec = which(bd[,qc_colNum] == qc_val)
        qcSamp = bd[,samp_cname][qcVec]
        return(qcSamp)
      }),
      qcEdata = purrr::map(qcSample,function(qcs){
        qc_edata = filt_qc_data1[,colnames(filt_qc_data1) %in% qcs]
        return(qc_edata)
      }),
      calc = purrr::map(qcEdata,function(edat){
        num_not_NA = rowSums(!is.na(edat))
        whichBad = which(num_not_NA < 6)
        badMolecules = omicsData$e_data[whichBad,][[pmartR::get_edata_cname(omicsData)]]
        return(badMolecules)
      }),
      begin_end_QC = purrr::map(batchDat,function(bd){
        begin_end_QC = purrr::map(batchDat,function(bd){
          minQC = bd[which(bd[,order_cname] == min(bd[,order_cname])),qc_cname]
          maxQC = bd[which(bd[,order_cname] == max(bd[,order_cname])),qc_cname]
          return(list(minQC,maxQC))
        })
      }))
  
  # see if we have any samples with too many samples still
  badMolecules = unique(unlist(batch_tib$calc,use.names = FALSE))
  if(length(badMolecules) > 0){
    stop(c(paste0("The following molecules have too few non-missing QC data points in at least one sample-batch \n combination. Please remove them prior to running bc_qcrlsc: ",'\n'),
           paste0(' - ', badMolecules,'\n')))
  }

  start_end_qc = sum(unlist(batch_tib$begin_end_QC) != qc_val)
  if(start_end_qc > 0){
    stop("The first and last sample run for each batch must be a QC sample")
  }
  
  # we cannot have values that equal 0
  if(any(!is.na(data) & data == 0)){
    stop("QC-RLSC cannot run on data that has been normalized without
         backtransforming the data.")
  }
  
  # begin the calculations -----------------------------------------------------
  # run get_param to get the bad features
  optParam <- get_params(omicsData,block_cname,qc_cname, qc_val,
                         order_cname,missing_thresh,rsd_thresh)
  
  # now we remove the bad feats from optParam
  omicsNoBad <- omicsData
  if(length(optParam$bad_feats) > 0){
    cfilter <- pmartR::custom_filter(omicsData,e_data_remove = optParam$bad_feats)
    omicsNoBad <- pmartR::applyFilt(cfilter,omicsData)
  }
  
  # run normalize_qcrlsc
  # supress the warning that a custom filter has already been applied
  pmartObj = suppressWarnings(normalize_qcrlsc(omicsData_filt = omicsNoBad,
                   optimal_params = optParam$final_ests,
                   block_cname = block_cname,
                   qc_cname = qc_cname,
                   qc_ind = qc_val,
                   backtransform = backtransform,
                   keep_qc = keep_qc))
  
  # reorder the data as it may have shifted around
  # find the original data order
  original_order = colnames(omicsData$e_data)
  # remove the QC samples that are no longer in the pmart object
  proper_order = c(pmartR::get_edata_cname(omicsData),original_order[original_order %in% pmartObj$f_data[,samp_cname]])
  # reorder by the proper order now
  pmartObj$e_data <- pmartObj$e_data[,proper_order]
  
  # update the batch column back
  if(keep_qc == FALSE){
    qc_samples <- which(omicsData$f_data[,qc_cname] == qc_val)
    original_batch_names <- original_batch_names[-qc_samples]
    pmartObj$f_data[,block_cname] <- original_batch_names
  }
  
  # pmart object creation time --------------------------------------------------

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
  
  # Add the group information to the group_DF attribute in the omicsData object.
  #attr(pmartObj, "group_DF") = attr(omicsData,"group_DF")
  
  # Update the data_info attribute.
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "qcrlsc_scaling",
    params = list()
  )
  
  # return pmart object
  return(pmartObj)
}




