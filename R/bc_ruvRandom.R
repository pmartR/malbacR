#' RUV-random normalized adjusted object
#' 
#' This function returns a pmart object that has been undergone for RUV-random
#' normalization regarding batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param nc_cname character string giving the name of the column in omicsData$e_meta
#'  that contains the factor variable indicating whether the sample is a negative control or not
#' @param nc_val value giving the name/value in the nc_cname column that indicates
#'  that the molecule is a negative control or not
#' @param k number of factors of unwanted variation, the default is 3, if the only
#'  factor of unwanted variation is the batch effect, we recommend setting k to be the number
#'  of batches in the data if visual inspection is not an option
#'  
#' @return Object of same class as omicsData that has been undergone
#'   RUV-random normalization
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_mix")
#' pmart_mix <- edata_transform(pmart_mix,"log2")
#' pmart_mix <- group_designation(pmart_mix,main_effects = "BatchNum",batch_id = "BatchNum")
#' pmart_mix <- normalize_global(pmart_mix,subset_fn = "all",norm_fn = "median",
#'                               apply_norm = TRUE,backtransform = TRUE)
#' mix_ruv <- bc_ruvRandom(omicsData = pmart_mix, nc_cname = "tag",nc_val = "IS", k = 3)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_ruvRandom <- function(omicsData,nc_cname,nc_val,k = 3) {
  
  # run checks -----------------------------------------------------------------
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # nc_cname
  if (class(nc_cname) != "character") {
    stop("Input parameter nc_cname must be of class 'character'.")
  }
  
  if (length(nc_cname) != 1) {
    stop("Input parameter nc_cname must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!any(names(omicsData$e_meta) == nc_cname)) {
    stop("Input parameter nc_cname must be a column found in e_meta of omicsData.")
  }
  
  if (length(nc_val) != 1) {
    stop("Input parameter nc_val must be of length 1 (e.g. vector containing a single element")
  }
  
  if (!(nc_val %in% omicsData$e_meta[, nc_cname])) {
    stop(
      "Input parameter nc_val must be a value found in the nc_cname column from the e_meta of omicsData."
    )
  }
  
  # which column has the edata cname
  cnameCol <- which(colnames(omicsData$e_data) == pmartR::get_edata_cname(omicsData))
  if (sum(is.na(omicsData$e_data[,-cnameCol]))!=0) {
    stop ("RUV-random requires no missing observations. Remove molecules with missing samples.")
  }

  # run ruv-random calculations ------------------------------------------------
  
  # find the values needed for the normalize function
  # find the parameter Y (the e_data abundance values)
  edat <- as.matrix(omicsData$e_data[,-1]) %>%
    t()
  molecules <- omicsData$e_data[,1]
  
  # find the parameter ctl (the negative controls)
  edat_cname = pmartR::get_edata_cname(omicsData)
  proper_order_ctl <- omicsData$e_data %>%
    dplyr::select(dplyr::all_of(edat_cname)) %>%
    dplyr::left_join(omicsData$e_meta)
  
  nc_cnameNum <- which(colnames(proper_order_ctl) == nc_cname)
  ctlRUV <- proper_order_ctl[,nc_cnameNum] == nc_val
  
  # now create the ruv-random abundance values!
  # we set lambda to be null, plotk to be FALSE
  ruvObj <- RUVRand(Y = edat,
                             ctl = ctlRUV,
                             k = k,
                             plotk = FALSE,
                             lambda = NULL)
  
  # we only need the newY values for abundance
  edatRUV <- ruvObj$newY
  
  # tranpose back and make it a data frame for pmart
  edatRUV <- edatRUV %>%
    t() %>%
    data.frame()
  edatRUV <- cbind(molecules,edatRUV)
  edat_colNames <- colnames(omicsData$e_data)
  colnames(edatRUV) <- edat_colNames

  # create pmart object --------------------------------------------------------
  
  # find the important values for pmart creation
  
  # find the individual data sets
  fdat = omicsData$f_data[omicsData$f_data$SampleID %in% colnames(edatRUV),]
  emet <- omicsData$e_meta
  
  # get the variables we need to create pmart object
  edat_cname = pmartR::get_edata_cname(omicsData)
  fdat_cname = pmartR::get_fdata_cname(omicsData)
  emet_cname = pmartR::get_emeta_cname(omicsData)
  
  # need to remove the negative control molecules
  # remove the negative controls from emeta
  bad_eggs <- emet %>%
    dplyr::filter(dplyr::across(dplyr::all_of(nc_cname)) == nc_val)
  emet <- emet %>%
    dplyr::filter(dplyr::across(dplyr::all_of(nc_cname)) != nc_val)

  # remove the negative controls from the edata
  to_remove = which(edatRUV[,edat_cname] %in% bad_eggs[,emet_cname])
  edatRUV <- edatRUV[-to_remove,]
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edatRUV,
                                          edata_cname = edat_cname,
                                          f_data = fdat,
                                          fdata_cname = fdat_cname,
                                          e_meta = emet,
                                          emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edatRUV,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edatRUV,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edatRUV,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"lipidDa}ta")){
    pmartObj = pmartR::as.lipidData(e_data = edatRUV,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edatRUV,
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
  
  # Update the data_info attribute.
  attributes(pmartObj)$data_info$batch_info <- list(
    is_bc = TRUE,
    bc_method = "ruv_random",
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
