#'Get Parameters for QC-RLSC Normalization
#'
#'Obtain parameters for quality control-based robust LOESS (locally estimated
#'scatterplot smoothing) signal correction
#'
#' @param omicsData an omicsData object (metabData, lipidData, pepData, or
#'  proData) created using the pmartR package, where any zero values have
#'  already been replaced with NAs and a log transformation has already been
#'  applied
#' @param block_cname character string giving name of column in omicsData$f_Data
#'  that contains block (or batch) information
#' @param qc_cname character string giving name of column in omicsData$f_data
#'  that contains the factor variable indicating whether sample is QC or not
#' @param qc_ind character string giving the value from the qc_cname column that
#'  indicates a QC sample
#' @param order_cname character string giving name of column in omicsData$f_data
#'  that contains the run order
#' @param missing_thresh numeric threshold, between 0 and 1, used for filtering
#'  out biomolecules. See details for more information. A value of 0.5 is
#'  reasonable.
#' @param rsd_thresh numeric threshold used for filtering metabolites. See
#'  details for more information. A value of 0.3 is reasonable.
#'
#' @author Lisa Bramer, Kelly Stratton
#' @references Dunn,W.B., Broadhurst,D., Begley,P., Zelena,E., Francis-McIntyre,S., Anderson,N., Brown,M.,
#' Knowles,J.D., Halsall,A., Haselden,J.N. et al. (2011) Procedures for
#' large-scale metabolic profiling of serum and plasma using gas chromatography and
#' liquid chromatography coupled to mass spectrometry. Nat. Protoc., 6, 1060-1083
#' @details Use this function to get the optimal parameter values to use in a
#'  subsequent call to \code{qcrlsc}, as well as a list of biomolecules that
#'  should be removed from the dataset prior to QC-RLSC normalization
#' @return list with elements final_ests (used as \code{optimal_params} input
#'  argument to \code{qcrlsc}) and bad_feats (features, or biomolecules, that
#'  need to be filtered out prior to using \code{qcrlsc})
get_params <-
  function(omicsData,
           block_cname,
           qc_cname,
           qc_ind,
           order_cname,
           missing_thresh,
           rsd_thresh) {
    ## initial checks ##
    
    # make sure data is of appropriate class #
    if (!class(omicsData) %in% c("metabData", "lipidData", "pepData", "proData")) {
      stop(
        "omicsData must be of class 'metabData', 'lipidData', 'pepData', or 'proData'. See package pmartR for details."
      )
    }
    
    # check that input parameters that should be character strings are indeed character strings, and each is a vector of length 1#
    if (class(block_cname) != "character") {
      stop("Input parameter block_cname must be of class 'character'.")
    }
    if (class(qc_cname) != "character") {
      stop("Input parameter qc_cname must be of class 'character'.")
    }
    if (class(qc_ind) != "character") {
      stop("Input parameter qc_ind must be of class 'character'.")
    }
    if (class(order_cname) != "character") {
      stop("Input parameter order_cname must be of class 'character'.")
    }
    
    if (length(block_cname) != 1) {
      stop(
        "Input parameter block_cname must be of length 1 (e.g. vector containing a single element)."
      )
    }
    if (length(qc_cname) != 1) {
      stop(
        "Input parameter qc_cname must be of length 1 (e.g. vector containing a single element)."
      )
    }
    if (length(qc_ind) != 1) {
      stop("Input parameter qc_ind must be of length 1 (e.g. vector containing a single element).")
    }
    if (length(order_cname) != 1) {
      stop(
        "Input parameter order_cname must be of length 1 (e.g. vector containing a single element)."
      )
    }
    
    # check that input parameters that should be numeric are indeed numeric and length 1 #
    if (class(missing_thresh) != "numeric") {
      stop("Input parameter missing_thresh must be numeric.")
    }
    if (class(rsd_thresh) != "numeric") {
      stop("Input parameter rsd_thres must be numeric.")
    }
    
    if (length(missing_thresh) != 1) {
      stop("Input parameter missing_thresh must be length 1.")
    }
    if (length(rsd_thresh) != 1) {
      stop("Input parameter rsd_thresh must be length 1.")
    }
    
    # check that missing_thresh is between 0 and 1 #
    if (missing_thresh < 0 |
        missing_thresh >= 1) {
      stop("Input parameter missing_thresh must be greater than or equal to 0 and less than 1.")
    }
    
    # check that the column names specified are in the f_data #
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
    
    # check that qc_ind is present in qc_cname #
    if (!(qc_ind %in% omicsData$f_data[, qc_cname])) {
      stop(
        "Input parameter qc_ind must be a value found in the qc_cname column from the f_data of omicsData."
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
    
    data <- omicsData$e_data
    f_data <- omicsData$f_data
    cnames <- attributes(omicsData)$cnames
    molecule_cname <- cnames$edata_cname
    samp_cname <- cnames$fdata_cname
    
    # make sure there are no 0's in e_data #
    edata <- data[, -which(colnames(data) == molecule_cname)]
    if (any(!is.na(edata) &
            edata == 0)) {
      stop(
        "Zero values in omicsData$e_data should be replaced with NAs prior to applying the Broadhurst normalization."
      )
    }
    
    
    # make sure the data has already been log transformed #
    if (!attributes(omicsData)$data_info$data_scale %in% c("log2", "log", "log10")) {
      stop(
        "omicsData$e_data should be transformed to the log2, log10, or natural log scale prior to applying the Broadhurst normalization."
      )
    }
    
    ## end of initial checks ##
    
    # separate out QC samples only from each data source #
    qc_cols = as.character(f_data[which(f_data[, qc_cname] == qc_ind), samp_cname])
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
    
    
    ## pull meta info ##
    
    # identify columns of f_data which contain relevant info #
    meta_col_ids = which(names(f_data) %in% c(samp_cname, block_cname, order_cname))
    
    # subset to columns #
    meta_data1 = f_data[, meta_col_ids]
    
    # subset to qc samples only #
    meta_data_qc = meta_data1[which(meta_data1[, samp_cname] %in% qc_cols),]
    
    final_params = list()
    
    ### parameter estimates ###
    
    # create a list of the molecules
    molCol = which(colnames(filt_qc_data1) == molecule_cname)
    molecules <- tibble::tibble(molecule = filt_qc_data1[,molCol])
    molecules <- molecules %>%
      dplyr::group_by(molecule, .add = TRUE) %>%
      dplyr::group_split()
    
    # create a tibble to hold all the data
    # create something for molecule name
    listostuff <- tibble::tibble(molecule = molecules)
    
    # create data frame for each molecules abundance
    listostuff <- listostuff %>%
      dplyr::mutate(
        # create data frame for each molecule abundance
        qc_list = split(filt_qc_data1,seq(nrow(filt_qc_data1))),
        # create data frame for long data
        temp_data = purrr::map(qc_list,function(pd){
          temp <- data.frame(Sample_ID = names(pd[, -which(names(pd) == molecule_cname)]),
                             Abundance = unlist(pd[, -which(names(pd) == molecule_cname)]))
          names(temp)[1] = samp_cname
          return(temp)
          }),
        # temp_data2 is for merge the meta_data_qc and temp_data by samp_cname
        temp_data2 = purrr::map(temp_data,function(td){
          temp = merge(x = meta_data_qc,
                       y = td,
                       by = samp_cname,
                       all.x = F,
                       all.y = T)
          return(temp)
          }),
        # temp_data3 is temp_data2 but now broken out by batch num
        temp_data3 = purrr::map(temp_data2,function(td2){
          temp = td2 %>%
            dplyr::group_by(!!as.name(block_cname)) %>%
            dplyr::select(tidyselect::all_of(samp_cname),tidyselect::all_of(block_cname),tidyselect::all_of(order_cname),Abundance) %>%
            dplyr::group_split()
          return(temp)
        }),
        # data2run is used for optimal span
        data2run = purrr::map(temp_data3,function(td3){
          purrr::map(td3, function(td3mini){
            td3mini %>%
              dplyr::select(tidyselect::all_of(order_cname),Abundance,tidyselect::all_of(block_cname)) %>%
              data.frame()
          })
        }))
    listostuff <- listostuff %>%
      dplyr::mutate(
        # cur_params uses the optimal spans function
        cur_params = purrr::map(data2run,function(d2r){
          purrr::map(d2r, function(d2rmini){
            #if(d2rmini[,2]) need to update it so that we only have min of 6 NA observations?
            if(sum(!is.na(d2rmini[,2])) >= 4) {
              cur_params = optimal_span(d2rmini[,-3])
            }
            else {
              cur_params = NA
            }
            return(cur_params)
          })
        }))

# for(i in 1:length(listostuff$data2run)){
#   for(j in 1:length(listostuff$data2run[[i]])){
#     if(sum(is.na(listostuff$data2run[[i]][[j]])) > 0){
#       print(c(i,j,sum(is.na(listostuff$data2run[[i]][[j]]))))
#     }
#   }
# }
    
    # 
    # listostuff <- listostuff %>%
    #   dplyr::mutate(
    #     naVal = purrr::map(data2run,function(d2r){
    #       purrr:::map(d2r,function(d2rmini){
    #         if(sum(!is.na(d2rmini[,2])) >= 4){return(TRUE)}
    #         else{return(FALSE)}
    #       })
    #     }))
    
    listostuff <- listostuff %>%
      dplyr::mutate(
        # param ests
        param_ests = purrr::pmap(list(data2run,qc_list,cur_params),function(d2r,qc,cp){
          purrr::map2(d2r,cp,function(d2rmini,cpmini){
            if(!is.na(cpmini$Lambda)){
              param_est = data.frame(
                Block = d2rmini[1,3],
                Metabolite = qc[,1],
                Span = cpmini$Alpha,
                Poly_Degree = cpmini$Lambda,
                MSE = cpmini$MSE
              )
            }
            else {
              param_est = data.frame(
                Block = d2rmini[1,3],
                Metabolite = qc[,1],
                Span = NA,
                Poly_Degree = NA,
                MSE = NA
              )
            }
            names(param_est)[which(names(param_est) == "Metabolite")] <- molecule_cname
            return(param_est)
          })
        }),
        # final params
        final_params = purrr::map(param_ests,function(pe){
          return(do.call(rbind,pe))
        })
        )

    # now we consolidate all of the final_ests together
    final_ests = do.call(rbind,listostuff$final_params)
    return(list(final_ests = final_ests, bad_feats = bad_feats))
  }
