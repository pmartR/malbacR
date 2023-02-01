#' TIGER filtering
#' 
#' This function returns a object of class 'tigerFilt' that has been filtered to ensure 
#' TIGER batch correction can run successfully
#' 
#' @param omicsData an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param sampletype_cname character string giving the name of the column in omicsData$f_data
#'  that contains the sample type information (such as quality control samples)
#' @param test_val character string giving the name of the value within the column sampletype_cname to be used
#'  as the testing value for TIGER
#'  
#' @return Object of tigerFilt class
#' 
#' @examples
#' \dontrun{
#' library(pmartRdata)
#' data(pmart_amide)
#' tigerFilt <- tiger_filter(omicsData = pmart_amide,sampletype_cname = "group", test_val = "QC")
#' pmart_amide_filt <- pmartR::applyFilt(filter_object = tigerFilt,omicsData = pmart_amide_filt)
#' amide_tiger <- bc_tiger(omicsData = pmart_amide_filt, sampletype_cname = "group", test_val = "QC")
#' }
#' 
#' @author Damon Leach
#' 
#' @export
#' 
tiger_filter <- function(omicsData,sampletype_cname,test_val){
  
  ## some initial checks ## ----------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData", "seqData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', or 'seqData'",
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
  
  # Extricate the names of the different ID columns.
  edata = omicsData$e_data
  fdata = omicsData$f_data
  edata_cname <- pmartR::get_edata_cname(omicsData)
  fdata_cname = pmartR::get_fdata_cname(omicsData)
  edat_cname_col <- which(colnames(edata) == edata_cname)
  if(!is.null(omicsData$e_meta)){
    emeta = omicsData$e_meta
    emeta_cname = pmartR::get_emeta_cname(omicsData)
  }else{emeta_cname = NULL}
  
  # we do not need the molecule column for this portion
  edata_abundance <- edata[,-edat_cname_col]
  
  # make sure that if fdata and edata are not in the same order, we at least
  # are ordering them correctly
  fdata_ordering <- fdata[,fdata_cname]
  fdata_qc_info <- fdata[,sampletype_cname]
  # set up the proper order of QC vs not QC samples
  proper_ordering <- match(colnames(edata_abundance),fdata_ordering)
  proper_qc_order <- fdata_qc_info[proper_ordering]
  
  # select QC only columns
  qc_samples <- edata_abundance[,proper_qc_order == test_val]
  # select nonQC only columns
  non_qc_samples <- edata_abundance[,proper_qc_order != test_val]
  
  # set up which columns will be used
  ij_combo <- c(NA,NA)
  # find the smallest amount of missingness based on that combination
  smallest_num <- nrow(qc_samples)
  for(i in 1:ncol(qc_samples)){
    for(j in 1:ncol(non_qc_samples)){
      missing_num <- sum(is.na(rowSums(data.frame(qc_samples[,i],non_qc_samples[,j]))))
      if(missing_num < smallest_num){
        smallest_num = missing_num
        ij_combo = c(i,j)
      }
    }
  }
  
  # qc sample name
  qc_sample_remover <- colnames(qc_samples)[ij_combo[1]]
  # non qc sample name
  non_qc_sample_remover <- colnames(non_qc_samples)[ij_combo[2]]
  
  # find how many molecules need to be removed across the two samples with least
  # amount of missing variables
  mini_edata <- edata[,c(edata_cname,qc_sample_remover,non_qc_sample_remover)]
  # mini_edata$QC16[4] <- NA
  # mini_edata$A35[c(5,10)] <- NA
  bad_mols <- mini_edata[which(rowSums(is.na(mini_edata)) != 0),edata_cname]
  
  filter_object <- list(e_data_remove = bad_mols)
  
  # output #
  output <- filter_object
  
  orig_class <- class(output)
  class(output) <- c("tigerFilt", orig_class)
  attr(output, "min_mol_removed") <- smallest_num
  
  attr(output,"samples_used") <- data.frame(mini_edata)
  
  # attributes #
  attr(output, "num_samples") = length(unique(omicsData$f_data[, fdata_cname]))
  attr(output, "num_edata") = length(unique(omicsData$e_data[, edata_cname]))
  attr(output, "num_emeta") = if(!is.null(emeta_cname)) length(unique(omicsData$e_meta[, emeta_cname]))
  
  attr(output, "cnames")$edata_cname = edata_cname
  attr(output, "cnames")$emeta_cname = emeta_cname
  attr(output, "cnames")$fdata_cname = fdata_cname
  
  attr(output, "omicsData") = omicsData # added 12/5/2017 by KS #
  
  return(output)
  
  # look at each possible pair of qc vs non qc and find one with biggest union observed
  # summary statistics if its own distinct filter
  # will need an applyFilt for TIGER
}

#' print.tigerFilt
#' 
#' For printing an S3 object of type 'tigerFilt':
#' 
#'@rdname print-tigerFilt
#'@export
#'
print.tigerFilt <- function(filter_object){
  if(!inherits(filter_object, "tigerFilt")) stop("filter object must be of the class 'tigerFilt'")
  edata_rm<- filter_object$e_data_remove
  
  if(length(edata_rm)< 8){
    out<- cbind(edata_rm)
    colnames(out)<- 'e_data_remove'
    cat(capture.output(out), sep ='\n')
    cat('\n')
  }else{
    edata_rm_hd<- cbind(head(edata_rm, 4))
    edata_rm_tl<- cbind(tail(edata_rm, 4))
    blank_row<- '---'
    
    out<-rbind(edata_rm_hd, blank_row, edata_rm_tl)
    colnames(out)<- 'e_data_remove'
    rownames(out)<-NULL
    cat(capture.output(out), sep = '\n')
    cat('\n')
  }
  
}

#' malbacR Filter Worker
#'@rdname danceR_filter_worker
#'
danceR_filter_worker <- function(filter_object, omicsData){
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  samp_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname
  
  # pull group_DF attribute #
  group_DF = attr(omicsData, "group_DF")
  
  # initialize the new omicsData parts #
  temp.edata <- omicsData$e_data
  temp.fdata <- omicsData$f_data
  temp.emeta <- omicsData$e_meta
  
  #check if filter object contains remove arguments
  if(!is.null(filter_object$edata_filt)|!is.null(filter_object$emeta_filt)|!is.null(filter_object$samples_filt)){
    
    # remove any samples from f_data and e_data #
    if(!is.null(filter_object$samples_filt)){
      inds <- which(temp.fdata[, which(names(temp.fdata) == samp_cname)] %in% filter_object$samples_filt)
      temp.fdata <- temp.fdata[-inds, ]
      
      inds <- which(names(temp.edata) %in% filter_object$samples_filt)
      temp.edata <- temp.edata[ ,-inds]
    }
    
    # remove any edata molecules from e_data and e_meta #
    if(!is.null(filter_object$edata_filt)){
      inds <- which(temp.edata[ , which(names(temp.edata) == edata_cname)] %in% filter_object$edata_filt)
      temp.edata <- temp.edata[-inds, ]
      
      # also remove these from e_meta, if it is present #
      if(!is.null(temp.emeta)){
        inds <- which(temp.emeta[ , which(names(temp.emeta) == edata_cname)] %in% filter_object$edata_filt)
        temp.emeta <- temp.emeta[-inds, ]
      }
    }
    
    # remove any emeta molecules from e_meta and e_data #
    if(!is.null(filter_object$emeta_filt)){
      inds <- which(temp.emeta[ , which(names(temp.emeta) == emeta_cname)] %in% filter_object$emeta_filt)
      if(length(inds) > 0){
        temp.emeta <- temp.emeta[-inds, ] 
      }
      
      # subset to the intersection of the edata_molecules in both e_data and e_meta, in case more were removed in one than the other #
      mols <- intersect(temp.edata[, which(names(temp.edata) == edata_cname)], temp.emeta[, which(names(temp.emeta) == edata_cname)])
      inds <- which(temp.edata[, which(names(temp.edata) == edata_cname)] %in% mols)
      temp.edata <- temp.edata[inds, ]
      inds <- which(temp.emeta[, which(names(temp.emeta) == edata_cname)] %in% mols)
      temp.emeta <- temp.emeta[inds, ]
    }
    
  }else{ # filter object contains keep arguments #
    
    # keep samples in f_data and e_data #
    if(!is.null(filter_object$samples_keep)){
      inds <- which(temp.fdata[, which(names(temp.fdata) == samp_cname)] %in% filter_object$samples_keep)
      temp.fdata <- temp.fdata[inds, ]
      
      inds <- c(which(names(temp.edata) == edata_cname), which(names(temp.edata) %in% filter_object$samples_keep))
      temp.edata <- temp.edata[ , inds]
    }
    
    # keep edata molecules in e_data #
    if(!is.null(filter_object$edata_keep)){
      inds <- which(temp.edata[ , which(names(temp.edata) == edata_cname)] %in% filter_object$edata_keep)
      temp.edata <- temp.edata[inds, ]
      
      # if e_meta is present and we aren't explicitly specifying to keep anything in it, also keep these e_data molecules in e_meta #
      if(!is.null(temp.emeta) & is.null(filter_object$emeta_keep)){
        inds <- which(temp.emeta[ , which(names(temp.emeta) == edata_cname)] %in% filter_object$edata_keep)
        temp.emeta <- temp.emeta[inds, ]
      }
    }
    
    # keep emeta molecules in e_meta (here, we are explicitly specifying things to keep) #
    if(!is.null(filter_object$emeta_keep)){
      inds <- which(temp.emeta[ , which(names(temp.emeta) == emeta_cname)] %in% filter_object$emeta_keep)
      temp.emeta <- temp.emeta[inds, ]
      
      # keep the union of the edata_molecules in both e_data and e_meta, in case more were kept in one than the other #
      if(is.null(filter_object$edata_keep)){
        # use intersection here, since nothing was explicitly specified to keep from edata, and if we use the union, then edata doesn't actually get filtered at all #
        mols <- intersect(temp.edata[, which(names(temp.edata) == edata_cname)], temp.emeta[, which(names(temp.emeta) == edata_cname)])
        inds <- which(temp.edata[, which(names(temp.edata) == edata_cname)] %in% mols)
        temp.edata <- temp.edata[inds, ]
        inds <- which(temp.emeta[, which(names(temp.emeta) == edata_cname)] %in% mols)
        temp.emeta <- temp.emeta[inds, ]
      }else{
        # use union here, since there WERE things explicitly specified to keep from edata #
        mols <- union(temp.edata[, which(names(temp.edata) == edata_cname)], temp.emeta[, which(names(temp.emeta) == edata_cname)])
        inds <- which(temp.edata[, which(names(temp.edata) == edata_cname)] %in% mols)
        temp.edata <- temp.edata[inds, ]
        inds <- which(temp.emeta[, which(names(temp.emeta) == edata_cname)] %in% mols)
        temp.emeta <- temp.emeta[inds, ]
      }
    }
  }
  
  # return the pieces needed to assemble a proData/pepData/lipidData/metabData object #
  output <- list(temp.pep2 = temp.edata, temp.samp2 = temp.fdata, temp.meta1 = temp.emeta, edata_cname = edata_cname, emeta_cname = emeta_cname, samp_cname = samp_cname)
  
  return(output)
}

#' Apply a S3 tiger filter object to a pmartR S3 object
#'
#' This function takes a filter object of class 'tigerFilt' and 
#'   applies the filter to a dataset of \code{pepData}, \code{proData}, 
#'   \code{lipidData}, \code{metabData}, or \code{nmrData}.
#'   
#' @param filter_object an object of the class 'tigerFilt' created by
#'   \code{tiger_filter}
#' @param omicsData an object of the class \code{pepData}, \code{proData},
#'   \code{lipidData}, \code{metabData}, or \code{nmrData} usually created by
#'   \code{\link{as.pepData}}, \code{\link{as.proData}},
#'   \code{\link{as.lipidData}}, \code{\link{as.metabData}}, or \code{\link{as.nmrData}}, respectively.
#'
#' @return An object of the class \code{pepData}, \code{proData},
#'   \code{lipidData}, \code{metabData}, or \code{nmrData} with specified cname_ids,
#'   edata_cnames, and emeta_cnames filtered out of the appropriate datasets.
#'
#'   
#' @examples
#' dontrun{
#' library(pmartRdata)
#' data("pmart_amide")
#' to_filter <- tiger_filter(omicsData = pmart_amide,sampletype_cname = "group", test_val = "QC")
#' amideFilt <- apply_tigerFilt(filter_object = to_filter, omicsData = pep_object)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton,Damon Leach
#'
#' @export
apply_tigerFilt <- function(filter_object, omicsData){
  
  
  ## some initial checks ## ----------------------------------------------------
  
  # check that omicsData is of appropriate class #
  if (!inherits(omicsData, c("pepData", "proData", "metabData", "lipidData",
                             "nmrData", "seqData"))) {
    
    stop (paste("omicsData must be of class 'pepData', 'proData', 'metabData',",
                "'lipidData', 'nmrData', or 'seqData'",
                sep = ' '))
  }
  
  # sampletype_cname - type of each sample
  if (!"tigerFilt" %in% class(filter_object)) {
    stop("filter_object must be of class tigerFilt.")
  }
  
  if(length(filter_object$e_data_remove) == 0){
    message("There are no molecules that need to be removed. Returning the original omicsData")
    return(omicsData)
  }else{
    edata_cname <- attributes(omicsData)$cnames$edata_cname
    fdata_cname <- attributes(omicsData)$cnames$fdata_cname
    emeta_cname <- attributes(omicsData)$cnames$emeta_cname
    
    # if filter_object contains 'removes' #
    filter_object_new = list(edata_filt = filter_object$e_data_remove, emeta_filt = filter_object$e_meta_remove, samples_filt = filter_object$f_data_remove)
    
    # check that edata_filt doesn't specify ALL the items in omicsData #
    if(all(omicsData$e_data[, edata_cname] %in% filter_object_new$edata_filt)){stop("e_data_remove specifies all the items in the data")}
    
    # call the function that does the filter application
    results_pieces <- danceR_filter_worker(omicsData = omicsData, filter_object = filter_object_new)
    
    # return filtered data object #
    results <- omicsData
    results$e_data <- results_pieces$temp.pep2
    results$f_data <- results_pieces$temp.samp2
    if(!is.null(omicsData$e_meta)){ # if-statement added by Kelly 3/24/2017 #
      results$e_meta <- data.frame(results_pieces$temp.meta1)
      names(results$e_meta)[which(names(omicsData$e_meta) == emeta_cname)] <- emeta_cname
    }else{
      # e_meta is null
      results$e_meta <- NULL
    }
    
    # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure #
    if(!is.null(attr(results, "group_DF"))){
      results <- pmartR::group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    }else{
      # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
      attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
      attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$prop_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
      attributes(results)$data_info$num_samps = ncol(results$e_data) - 1
      
      if(!is.null(results$e_meta)){
        # number of unique proteins that map to a peptide in e_data #
        if(!is.null(emeta_cname)){
          num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
        }else{num_emeta = NULL}
      }else{
        num_emeta = NULL
      }
      attr(results, "data_info")$num_emeta = num_emeta
      ## end of update attributes (7/11/2016 by KS)
    }
    
    # set attributes for which filters were run
    attr(results, "filters")$tigerFilt <- list(report_text = "", threshold = c(), filtered = c())
    n_edata_filtered <- nrow(omicsData$e_data) - nrow(results$e_data)
    n_fdata_filtered <- nrow(omicsData$f_data) - nrow(results$f_data)
    if(!is.null(omicsData$e_meta)){
      n_emeta_filtered <- nrow(omicsData$e_meta) - nrow(results$e_meta)
    }else{
      n_emeta_filtered = NA
    }
    if(!is.null(omicsData$e_meta)){
      attr(results, "filters")$tigerFilt$report_text <- paste("A tiger filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data, ", n_emeta_filtered, " ", emeta_cname, "s from e_meta, and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }else{
      attr(results, "filters")$tigerFilt$report_text <- paste("A tiger filter was applied to the data, removing ", n_edata_filtered, " ", edata_cname, "s from e_data and ", n_fdata_filtered, " ", fdata_cname, "s from f_data.", sep="")
    }
    attr(results, "filters")$tigerFilt$filtered <- filter_object
    return(results)
  }
}
