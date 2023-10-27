#' WaveICA_2.0 batch correction
#' 
#' This function returns a pmart object that has been undergone WaveICA_2.0
#'  batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param injection_cname character string giving the name of the column in omicsData$f_data
#'  that contains the run order of the samples
#' @param alpha tradeoff value between the independence of samples and those of variables in ICA,
#'  and should be between 0 and 1, defaults to 0
#' @param cutoff threshold of the variation explained by the injection order for independent components,
#'  should be between 0 and 1, defaults to 0.1
#' @param K maximal component that ICA decomposes, defaults to 10
#'  
#' @return Object of same class as omicsData that has been undergone
#'   WaveICA normalization
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
#'                                apply_norm = TRUE,backtransform = TRUE)
#' impObj <- imputation(omicsData = pmart_amide)
#' amide_imp <- apply_imputation(imputeData = impObj, omicsData = pmart_amide)
#' amide_wave <- bc_waveica(omicsData = amide_imp, injection_cname = "Injection_order",
#'                          alpha = 0, cutoff = 0.1, K = 10)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_waveica <- function(omicsData,injection_cname, alpha = 0, cutoff = 0.1, K = 10){
  
  # run through checks ---------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # injection_cname
  if (class(injection_cname) != "character") {
    stop("Input parameter injection_cname must be of class 'character'.")
  }
  
  if (length(injection_cname) != 1) {
    stop("Input parameter injection_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$f_data) == injection_cname)) {
    stop("Input parameter injection_col must be a column found in f_data of omicsData.")
  }
  
  # alpha
  # must only be one number
  if (length(alpha) != 1) {
    stop("Input parameter alpha must be of length 1 (e.g. vector containing a single element")
  }
  # alpha must be numeric
  if (class(alpha) != "numeric") {
    stop("Input parameter alpha must be of class 'numeric'")
  }
  
  # alpha must also be between 0 and 1
  if (alpha < 0 | alpha > 1) {
    stop("Input parameter alpha must be between 0 and 1 inclusive, the default is 0")
  }
  
  # cutoff
  # must only be one number
  if (length(cutoff) != 1) {
    stop("Input parameter cutoff must be of length 1 (e.g. vector containing a single element")
  }
  # cutoff must be numeric
  if (class(cutoff) != "numeric") {
    stop("Input parameter cutoff must be of class 'numeric'")
  }
  
  # cutoff must also be between 0 and 1
  if (cutoff < 0 | cutoff > 1) {
    stop("Input parameter cutoff must be between 0 and 1 inclusive, the default is 0.1")
  }
  
  # K
  # must only be one number
  if (length(K) != 1) {
    stop("Input parameter K must be of length 1 (e.g. vector containing a single element")
  }
  # K must be numeric
  if (class(K) != "numeric") {
    stop("Input parameter K must be of class 'numeric'")
  }
  
  # k must also be be greater than 0
  if (K <= 0) {
    stop("Input parameter K must be greater than 0, the default is 10")
  }
  
  # WaveICA requires no missing values
  # which column has the edata cname
  cnameCol <- which(colnames(omicsData$e_data) == pmartR::get_edata_cname(omicsData))
  if (sum(is.na(omicsData$e_data[,-cnameCol])) != 0) {
    stop ("WaveICA requires no missing observations. Remove molecules with missing samples.")
  }
  
  # run the WaveICA2.0 calculations -------------------------------------------------
  
  # find the values needed for the normalize function
  # find the parameter Y (the e_data abundance values)
  edat <- as.matrix(omicsData$e_data[,-1]) %>%
    t() %>%
    data.frame()
  molecules <- omicsData$e_data[,1]
  colnames(edat) <- molecules
  
  # make sure fdata sampleid and edat are in the same order
  # find the injection order
  injection_orderCol <- which(colnames(omicsData$f_data) == injection_cname)
  injection_order <- omicsData$f_data[,injection_orderCol]
  # run WaveICA 2.0
  set.seed(1)
  waveica <- WaveICA2.0::WaveICA_2.0(data = edat, Injection_Order = injection_order,K = 10, alpha = 0, Cutoff = 0.1)
  
  edatWAVE <- waveica %>%
    data.frame() %>%
    t() %>%
    data.frame()
  
  edatWAVE <- cbind(molecules,edatWAVE)
  edat_colNames <- colnames(omicsData$e_data)
  colnames(edatWAVE) <- edat_colNames
  rownames(edatWAVE) <- NULL
  
  # move back into pmart object ------------------------------------------------
  
  # find the important values for pmart creation
  edat_cname = pmartR::get_edata_cname(omicsData)
  fdat_cname = pmartR::get_fdata_cname(omicsData)
  fdat = omicsData$f_data[omicsData$f_data[[fdat_cname]] %in% colnames(edatWAVE),]

  # assume emeta is NULL unless otherwise stated
  emet_cname = NULL
  emet = NULL
  if(!is.null(omicsData$e_meta)){
    # need to update the e_meta if it exists as NOMIS removes the IS values
    emet_cname = pmartR::get_emeta_cname(omicsData)
    emet_cnameNum = which(colnames(omicsData$e_meta) == emet_cname)
    emet = omicsData$e_meta[omicsData$e_data[,emet_cnameNum] %in% edatWAVE[,1],]
  }
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edatWAVE,
                                          edata_cname = edat_cname,
                                          f_data = fdat,
                                          fdata_cname = fdat_cname,
                                          e_meta = emet,
                                          emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edatWAVE,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edatWAVE,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edatWAVE,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edatWAVE,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edatWAVE,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  
  # Update the data_info attribute.
  attr(pmartObj, 'data_info') <- pmartR:::set_data_info(
    e_data = pmartObj$e_data,
    edata_cname = pmartR::get_edata_cname(omicsData),
    data_scale_orig = pmartR::get_data_scale_orig(omicsData),
    data_scale = pmartR::get_data_scale(omicsData),
    data_types = pmartR::get_data_info(omicsData)$data_types,
    norm_info = pmartR::get_data_info(omicsData)$norm_info,
    is_normalized = TRUE,
    batch_info = pmartR::get_data_info(omicsData)$batch_info,
    is_bc = pmartR::get_data_info(omicsData)$batch_info$is_bc
  )
  
  # Add the group information to the group_DF attribute in the omicsData object.
  attr(pmartObj, "group_DF") = attr(omicsData,"group_DF")
  
  # Update the data_info attribute.
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "waveica",
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
