#' QC-RLSC Normalization
#'
#' Quality control-based robust LOESS (locally estimated scatterplot smoothing)
#' signal correction
#'
#' @param omicsData_filt an omicsData object (metabData, lipidData,
#'   pepData, or proData) created using the pmartR package, where any zero
#'   values have already been replaced with NAs and a log transformation has
#'   already been applied, where any biomolecules that were identified by
#'   \code{get_params()} as needing to be removed have been removed
#' @param optimal_params final_ests element of the output from \code{get_params()} function (this is
#'   needed for Span and Poly_Degree values)
#' @param block_cname character string giving name of column in omicsData_filt$f_Data
#'  that contains block (or batch) information. Values in this column must be numeric.
#' @param qc_cname character string giving name of column in omicsData_filt$f_data
#'  that contains the factor variable indicating whether sample is QC or not
#' @param qc_ind character string giving the value from the qc_cname column that
#'  indicates a QC sample
#' @param backtransform logical value indicated whether or not to backtransform
#'   the data to put the normalized data back on the same scale as the original
#'   data. Defaults to FALSE. If TRUE, the median of the y-values of the loess
#'   curve is added to the normalized value for each biomolecule for each batch.
#' @author Kelly Stratton, Lisa Bramer
#' @references Dunn,W.B., Broadhurst,D., Begley,P., Zelena,E., Francis-McIntyre,S., Anderson,N., Brown,M.,
#' Knowles,J.D., Halsall,A., Haselden,J.N. et al. (2011) Procedures for
#' large-scale metabolic profiling of serum and plasma using gas chromatography and
#' liquid chromatography coupled to mass spectrometry. Nat. Protoc., 6, 1060-1083
#' @details This function applies the QC-RLSC normalization to each batch of input data. Must use \code{get_params} function first in order to get optimal_params input for \code{qcrlsc}.
#' @return omicsData object of same class as omicsData_filt, where e_data contains the QC-RLSC normalized values
#' @rdname normalize_qcrlsc
#' @name normalize_qcrlsc
normalize_qcrlsc <-
  function(omicsData_filt,
           optimal_params,
           block_cname,
           qc_cname,
           qc_ind,
           backtransform = FALSE) {
    ### initial checks ###

    if (!class(omicsData_filt) %in% c("metabData", "lipidData", "pepData", "proData")) {
      stop(
        "omicsData_filt must be an S3 object of class 'metabData', 'lipidData', 'pepData', or 'proData'. See pmartR package for more details."
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

    # check that column names specified exist in f_data #
    if (!any(names(omicsData_filt$f_data) == block_cname)) {
      stop("The f_data component of omicsData_filt must contain a column for 'Batch'.")
    }
    if (!any(names(omicsData_filt$f_data) == qc_cname)) {
      stop("The f_data component of omicsData_filt must contain a column for 'QCtype'.")
    }

    # check that qc_ind is present in qc_cname #
    if (!(qc_ind %in% omicsData_filt$f_data[, qc_cname])) {
      stop(
        "The value for qc_ind is not present in the column specified by qc_cname (in f_data component of omicsData_filt)."
      )
    }

    # check that optimal params is data.frame with expected elements #
    if (!is.data.frame(optimal_params)) {
      stop(
        "Input parameter optimal_params must be a data.frame. See 'get_params' function documentation for more detail."
      )
    }
    if (!(all(c("Span", "Poly_Degree", "MSE") %in% names(optimal_params)))) {
      stop(
        "Input parameter optimal_params must be the final_ests element from the output from 'get_params' function. See documentation for more detail."
      )
    }

    # check that values in block_cname are numeric #
    if (!is.numeric(omicsData_filt$f_data[, block_cname])) {
      stop(
        "Values in block_cname column (in f_data component of omicsData_filt) must be numeric."
      )
    }

    # check that backtransform is logical #
    if (!is.logical(backtransform)) {
      stop("Input parameter backtransform must be of class 'logical'.")
    }

    ### end of initial checks ###


    e_data <- omicsData_filt$e_data
    f_data <- omicsData_filt$f_data

    cnames <- attributes(omicsData_filt)$cnames
    molecule_cname <- cnames$edata_cname
    samp_cname <- cnames$fdata_cname

    # get medians of each metabolite across all samples (across all batches) #
    metab_medians <- apply(e_data[, -1], 1, median, na.rm = TRUE)

    # initialize list #
    gold_std_batch_i_norm <- list()

    for (i in 1:length(unique(f_data[, block_cname]))) {
      # restrict to batch i #
      samp_id_ind <- which(names(f_data) == samp_cname)

      ids <-
        as.character(f_data[which(f_data[, block_cname] == i), samp_id_ind])

      batch_i <- e_data[, which(names(e_data) %in% ids)] # note: we lose the biomolecule identifier column here, need to add it back in later

      biomolecules_i <- e_data[, which(names(e_data) == molecule_cname)]

      ids <-
        as.character(f_data[which(f_data[, qc_cname] == qc_ind &
                                    f_data[, block_cname] == i), samp_id_ind])
      qcs_i <- e_data[, which(names(e_data) %in% ids)]

      optimal_params_i <- optimal_params[optimal_params$Block == i, ]

      # initialize gold_std_batch_i_norm #
      gold_std_batch_i_norm[[i]] <-
        matrix(NA, nrow = nrow(batch_i), ncol = ncol(batch_i))
      # dim(gold_std_batch_i_norm[[i]])

      for (j in 1:nrow(batch_i)) {
        # fit loess curve for each metab #
        x <-
          which(f_data[, qc_cname] == qc_ind &
                  f_data[, block_cname] == i)  # isolate the QC indices for the ith batch
        y.qcs <- as.numeric(qcs_i[j,])
        if (!is.na(optimal_params_i$Poly_Degree[j]) &
            sum(is.na(y.qcs)) < length(y.qcs) / 2) {
          y.loess <-
            loess(
              y.qcs ~ x,
              span = optimal_params_i$Span[j],
              degree = optimal_params_i$Poly_Degree[j],
              family = "gaussian"
            )
          y <- batch_i[j, ]
          y.predict <-
            predict(y.loess, data.frame(x = which(f_data[, block_cname] == i)))

          ## apply the normalization ##
          if (backtransform == TRUE) {
            # add back the median of y.loess, in order to get the data back on the original scale of the data
            y.norm <- as.numeric(y - y.predict + metab_medians[j])
          } else{
            y.norm <- as.numeric(y - y.predict)
          }

        } else{
          # catch cases where there were not at least 4 nonmissing values in the data
          y.norm <- NA
        }


        gold_std_batch_i_norm[[i]][j,] <- y.norm
      }
      gold_std_batch_i_norm[[i]] <-
        data.frame((e_data[which(names(e_data) == molecule_cname)]), gold_std_batch_i_norm[[i]])
      gold_std_batch_i_norm[[i]] <-
        as.data.frame(gold_std_batch_i_norm[[i]])
      names(gold_std_batch_i_norm[[i]]) <-
        c(molecule_cname, names(batch_i))

    }

    norm_data <- gold_std_batch_i_norm
    norm_data2 <- norm_data
    # make sure norm_data is numeric (it is defaulting to factor...) #
    for(i in 1:length(norm_data)){
      norm_data2[[i]][, -1] <- apply(norm_data[[i]][, -1], 2, as.numeric)
    }

    # put normalized data into metabData object #
    output_data <- omicsData_filt # original omicsData object

    # if > 1 batch, then the molecule cname will be repeated #
    if(length(unique(f_data[, block_cname])) > 1){
      inds <- grep(molecule_cname, names(as.data.frame(norm_data2)))
      # keep 1st occurrence but remove the remaining

      inds_rm <- inds[2:length(inds)]
      output_data$e_data <- as.data.frame(norm_data2)[, -inds_rm]
    }else{
      # single batch
      output_data$e_data <- as.data.frame(norm_data2)
    }


    attributes(output_data)$data_info$data_norm <- TRUE

    # filter out QC samples #
    inds <- which(output_data$f_data[, qc_cname] == qc_ind)
    to_remove <- output_data$f_data[inds, samp_cname]

    myfilter <- pmartR::custom_filter(output_data, f_data_remove = to_remove)
    output_data2 <- pmartR::applyFilt(myfilter, output_data)

    return(output_data2)
  }
