#' WaveICA batch correction
#' 
#' This function returns a pmart object that has been undergone either WaveICA or WaveICA_2.0
#'  batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param batach_cname character string giving the name of the colum in omicsData$f_data
#'  that contains the batch information of the samples (to be used for WaveICA)
#' @param injection_cname character string giving the name of the column in omicsData$f_data
#'  that contains the run order of the samples (to be used for WaveICA2.0)
#' @param version character string that specifies whether the user wants to run "WaveICA"
#'  or "WaveICA2.0" (default is "WaveICA")
#' @param alpha tradeoff value between the independence of samples and those of variables in ICA,
#'  and should be between 0 and 1, defaults to 0
#' @param cutoff_injection threshold of the variation explained by the injection order for independent components,
#'  should be between 0 and 1, defaults to 0.1 (to be used for WaveICA2.0)
#' @param cutoff_batch threshold of the variation explained by the batch for independent components,
#'  should be between 0 and 1, defaults to 0.05 (to be used for WaveICA)
#' @param cutoff_group threshold of the variation explained by the group for independent components,
#'  should be between 0 and 1, defaults to 0.05 (to be used for WaveICA)
#' @param K maximal component that ICA decomposes, defaults to 20 (though recommendation for K = 10
#'  when using "WaveICA2.0")
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
#' impObj <- imputation(omicsData = pmart_amide)
#' amide_imp <- apply_imputation(imputeData = impObj, omicsData = pmart_amide)
#' amide_wave <- bc_waveica(omicsData = amide_imp, injection_cname = "Injection_order",
#'                          version = "WaveICA2.0",alpha = 0, cutoff_injection = 0.1, K = 10)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_waveica <- function(omicsData,batch_cname = NULL, injection_cname = NULL, version = "WaveICA", alpha = 0, cutoff_injection = 0.1, K = 20,
                       cutoff_batch = 0.05, cutoff_group = 0.05){
  # run through checks ---------------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that whether we are running WaveICA or WaveICA2.0
  if (class(version) != "character") {
    stop("Input parameter version must be of class 'character'")
  }
  if (!version %in% c("WaveICA","WaveICA2.0")) {
    stop("Input parameter version must be either 'WaveICA' or 'WaveICA2.0'")
  }
  
  # injection_cname (matters for WaveICA2.0)
  if (class(injection_cname) != "character" & version == "WaveICA2.0") {
    stop("Input parameter injection_cname must be of class 'character' for WaveICA2.0.")
  }
  
  if (length(injection_cname) != 1 & version == "WaveICA2.0") {
    stop("Input parameter injection_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$f_data) == injection_cname) & version == "WaveICA2.0") {
    stop("Input parameter injection_cname must be a column found in f_data of omicsData.")
  }
  
  # batch_cname (matters for WaveICA)
  if (class(batch_cname) != "character" & version == "WaveICA") {
    stop("Input parameter batch_cname must be of class 'character' for WaveICA")
  }
  
  if (length(batch_cname) != 1 & version == "WaveICA") {
    stop("Input parameter batch_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$f_data) == batch_cname) & version == "WaveICA") {
    stop("Input parameter batch_cname must be a column found in f_data of omicsData.")
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
  
  # cutoff_injection
  # must only be one number
  if (length(cutoff_injection) != 1 & version == "WaveICA2.0") {
    stop("Input parameter cutoff_injection must be of length 1 (e.g. vector containing a single element")
  }
  # cutoff_injection must be numeric
  if (class(cutoff_injection) != "numeric" & version == "WaveICA2.0") {
    stop("Input parameter cutoff_injection must be of class 'numeric'")
  }
  
  # cutoff_injection must also be between 0 and 1
  if (cutoff_injection < 0 | cutoff_injection > 1 & version == "WaveICA2.0") {
    stop("Input parameter cutoff_injection must be between 0 and 1 inclusive, the default is 0.1")
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
  
  # cutoff_batch
  # must only be one number
  if (length(cutoff_batch) != 1 & version == "WaveICA") {
    stop("Input parameter cutoff_batch must be of length 1 (e.g. vector containing a single element")
  }
  # cutoff_batch must be numeric
  if (class(cutoff_batch) != "numeric" & version == "WaveICA") {
    stop("Input parameter cutoff_batch must be of class 'numeric'")
  }
  
  # cutoff_batch must also be between 0 and 1
  if (cutoff_batch < 0 | cutoff_batch > 1 & version == "WaveICA") {
    stop("Input parameter cutoff_batch must be between 0 and 1 inclusive, the default is 0.1")
  }
  
  # cutoff_group
  # must only be one number
  if (length(cutoff_group) != 1 & version == "WaveICA") {
    stop("Input parameter cutoff_group must be of length 1 (e.g. vector containing a single element")
  }
  # cutoff_group must be numeric
  if (class(cutoff_group) != "numeric" & version == "WaveICA") {
    stop("Input parameter cutoff_group must be of class 'numeric'")
  }
  
  # cutoff_group must also be between 0 and 1
  if (cutoff_group < 0 | cutoff_group > 1 & version == "WaveICA") {
    stop("Input parameter cutoff_group must be between 0 and 1 inclusive, the default is 0.1")
  }
  
  # WaveICA requires no missing values
  # which column has the edata cname
  cnameCol <- which(colnames(omicsData$e_data) == pmartR::get_edata_cname(omicsData))
  if (sum(is.na(omicsData$e_data[,-cnameCol])) != 0) {
    stop ("WaveICA requires no missing observations. Remove molecules with missing samples.")
  }
  
  # run the WaveICA2.0 calculations -------------------------------------------------
  # retain seed after  running code
  if (!exists(".Random.seed")) runif(1)
  old_seed <- .Random.seed
  on.exit(.Random.seed <- old_seed)
  
  # find the values needed for the normalize function
  # find the parameter Y (the e_data abundance values)
  edat <- as.matrix(omicsData$e_data[,-cnameCol]) %>%
    t() %>%
    data.frame()
  molecules <- omicsData$e_data[,cnameCol]
  colnames(edat) <- molecules
  
  
  if(version == "WaveICA"){
    batch_valuesCol <- which(colnames(omicsData$f_data) == batch_cname)
    batch_values <- omicsData$f_data[,batch_cname]
    set.seed(1)
    waveica <- WaveICA(data = edat, wf = "haar",batch = batch_values,K = K, t = cutoff_batch, t2 = cutoff_group, alpha = alpha)
  } else {
    # make sure fdata sampleid and edat are in the same order
    # find the injection order
    injection_orderCol <- which(colnames(omicsData$f_data) == injection_cname)
    injection_order <- omicsData$f_data[,injection_orderCol]
    # run WaveICA 2.0
    set.seed(1)
    waveica <- WaveICA2.0::WaveICA_2.0(data = edat, Injection_Order = injection_order,K = K, alpha = alpha, Cutoff = cutoff_injection)
  }
  
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
    is_normalized = pmartR::get_data_info(omicsData)$norm_info$is_normalized,
    batch_info = pmartR::get_data_info(omicsData)$batch_info,
    is_bc = pmartR::get_data_info(omicsData)$batch_info$is_bc
  )
  
  # Add the group information to the group_DF attribute in the omicsData object.
  attr(pmartObj, "group_DF") = attr(omicsData,"group_DF")
  
  # Update the data_info attribute for batch
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "bc_waveica",
    params = list(version = version,
                  injection_cname = injection_cname,
                  batch_cname = batch_cname,
                  cutoff_injection = cutoff_injection,
                  cutoff_batch = cutoff_batch,
                  cutoff_group = cutoff_group,
                  alpha = alpha,
                  K = K)
  )
  
  # update normalization as well 
  attributes(pmartObj)$data_info$norm_info <- list(
    is_normalized = TRUE,
    norm_type = "bc_waveica"
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
