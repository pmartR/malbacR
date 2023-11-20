#' EigenMS adjusted object
#' 
#' This function returns a pmart object that has been undergone for EigenMS
#' normalization regarding batch effect correction
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @return Object of same class as omicsData that has been undergone
#'   EigenMS normalization
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data("pmart_amide")
#' pmart_amide <- edata_transform(pmart_amide,"log2")
#' pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
#' amide_eigen <- bc_eigenMS(omicsData = pmart_amide)
#' 
#' @author Damon Leach
#' 
#' @export
#' 
bc_eigenMS <- function(omicsData){
  
  # run original data checks ---------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  # check that omicsData has batch id information
  if (is.null(attributes(attr(omicsData,"group_DF")))){
    stop (paste("omicsData must have group_designation ran",
                sep = ' '))
  }
  
  # find the eigenMS values ----------------------------------------------------
  # set up the m parameter (the e_data)
  cnameCol = which(colnames(omicsData$e_data) == pmartR::get_edata_cname(omicsData))
  orig_edat_rows = omicsData$e_data[,cnameCol]
  edat = omicsData$e_data[,-cnameCol]
  # maintain original colname order to compare emeta in the future
  orig_edat_cols = colnames(omicsData$e_data)
  
  # set up the treatment parameter as the main effects of group_designation()
  replicateCol = attr(omicsData,"group_DF")$Group
  replicate = as.factor(replicateCol)
  
  # need to make sure the fdata and edata are aligned
  fdat_cname_num = which(colnames(omicsData$f_data) == pmartR::get_fdata_cname(omicsData))
  rep_fdata = omicsData$f_data[,fdat_cname_num]
  edat = edat[,rep_fdata]

  # check to see if we even have enough data for all the molecules
  edat_transpose <- data.frame(t(edat))
  edat_transpose$replicate <- replicate
  # split the data by replicate
  split_dat <- split(edat_transpose,edat_transpose$replicate)
  # find the number of bad molecules
  column_sum_info = purrr::map(split_dat,function(sd){
    numnum = which(colnames(sd) == "replicate")
    sd <- sd[-numnum]
    good_data = colSums(!is.na(sd))
    return(sum(good_data == 0))
  })
  # find the total number of bad molecules, if we have more than 0 we send out error
  total_bad = sum(unlist(column_sum_info,use.names = FALSE))
  if(total_bad > 0){
    stop("omicsData must have molecule filter applied with use_groups = TRUE and min_num = 1")
  }
  
  # set up the prot.info parameter
  # use emeta data or create a meta data frame if there is none
  # if(!is.null(omicsData$e_meta)) {
  #   meta = omicsData$e_meta
  #   # need to know if the edata and emeta are in the same order
  #   emeta_cname = pmartR::get_emeta_cname(omicsData)
  #   emeta_cname_num = which(colnames(omicsData$e_meta)== pmartR::get_emeta_cname(omicsData))
  #   edat_emet_ordering = match(meta[,emeta_cname_num],orig_edat_rows)
  #   meta = meta[edat_emet_ordering,]
  #   meta <- meta %>%
  #     dplyr::relocate(dplyr::all_of(emeta_cname))
  #     
  # } else {
  #   meta = data.frame(omicsData$e_data[,cnameCol])
  # }
  meta = data.frame(omicsData$e_data[,cnameCol])
  
  # we need to set a seed so results remain constant each time we run it
  set.seed(12345)
  # find the bias trends
  # supress the warning that we need to set seed (it does it automatically unfortunately)
  bias_trends = R.devices::suppressGraphics(suppressWarnings(ProteoMM::eig_norm1(m = edat,
                   treatment = replicate,
                   prot.info = meta)))
  # run the eigenMS normalization
  eigen_ms = R.devices::suppressGraphics(ProteoMM::eig_norm2(rv = bias_trends))
  
  # put everything back into pmart object --------------------------------------
  
  # put edata into proper format
  edat_eigen = data.frame(eigen_ms$norm_m)
  edat_eigen <- tibble::rownames_to_column(edat_eigen, var = pmartR::get_edata_cname(omicsData))
  # reorder back to original ordering from user pmart object
  colEdatNames = colnames(edat)
  colnames(edat_eigen) <- c(pmartR::get_edata_cname(omicsData),colEdatNames)
  edat_eigen <- edat_eigen[,orig_edat_cols]

  # values needed for creating pmart object
  edat_cname = pmartR::get_edata_cname(omicsData)
  fdat = omicsData$f_data
  fdat_cname = pmartR::get_fdata_cname(omicsData)
  emet =  NULL
  emet_cname = NULL
  # need to update emeta if the data since we may have removed samples in the process
  # that do not have enough molecules
  if(!is.null(omicsData$e_meta)){
    emet = omicsData$e_meta
    emet_cname = pmartR::get_emeta_cname(omicsData)
    # emet = eigen_ms$normalized[,colnames(eigen_ms$normalized) %in% colnames(omicsData$e_meta)]
    # rownames(emet) = NULL
    # emet_cname = pmartR::get_emeta_cname(omicsData)
    # emeta_cname_num = which(colnames(omicsData$e_meta)== pmartR::get_emeta_cname(omicsData))
    # nrow_emet = nrow(emet)
    # # put back in the right order
    # emet = emet[rank(match(emet[,emet_cname],omicsData$e_meta[,emet_cname])),]
    # emet <- emet %>%
    #   dplyr::select(dplyr::all_of(colnames(omicsData$e_meta)))
  }
  
  # create the pmart object #
  if(inherits(omicsData,"isobaricpepData") & inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.isobaricpepData(e_data = edat_eigen,
                                          edata_cname = edat_cname,
                                          f_data = fdat,
                                          fdata_cname = fdat_cname,
                                          e_meta = emet,
                                          emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"pepData")){
    pmartObj = pmartR::as.pepData(e_data = edat_eigen,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"proData")){
    pmartObj = pmartR::as.proData(e_data = edat_eigen,
                                  edata_cname = edat_cname,
                                  f_data = fdat,
                                  fdata_cname = fdat_cname,
                                  e_meta = emet,
                                  emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"metabData")){
    pmartObj = pmartR::as.metabData(e_data = edat_eigen,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"lipidData")){
    pmartObj = pmartR::as.lipidData(e_data = edat_eigen,
                                    edata_cname = edat_cname,
                                    f_data = fdat,
                                    fdata_cname = fdat_cname,
                                    e_meta = emet,
                                    emeta_cname = emet_cname)
  }
  else if(inherits(omicsData,"nmrData")){
    pmartObj = pmartR::as.nmrData(e_data = edat_eigen,
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
    bc_method = "bc_eigenMS",
    params = list()
  )
  # update normalization as well 
  attributes(pmartObj)$data_info$norm_info <- list(
    is_normalized = TRUE,
    norm_type = "bc_eigenMS"
  )
  # Update the meta_info attribute.
  attr(pmartObj, 'meta_info') <- pmartR:::set_meta_info(
    e_meta = omicsData$e_meta,
    emeta_cname = pmartR::get_emeta_cname(omicsData)
  )
  
  # Update the filtesr attribute
  attr(pmartObj, 'filters') <- attr(omicsData,'filters')
  
  # return pmart object
  return(pmartObj)
}
