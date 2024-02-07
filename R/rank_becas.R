#' Rank Batch Corrected Data
#' 
#' This function returns a data.frame that ranks the batch corrected data supplied
#' from first to last based on one of the three different metrics, difference in 
#' R2m/R2c, coefficient of variation, and median distance of centroids of batch clusters
#' 
#' @param omicsData_beca_list an list containing at least 2 objects
#'  of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.
#' @param comparison_method a character string that can take on one of the following
#' three options: r2_diff, cv, or distance_pca corresponding to the metrics
#' difference in R2m/R2c, coefficient of variation, and median distance of centroids of batch
#' clusters respectively
#' @param batch_effect_cname a character string giving the name of the column in f_data of omics objects
#'  that contians the batch information
#' @param main_effect_cname a character string giving the name of the column in f_data of omics objects
#'  that contians the treatment (main effect) information
#' @param omicsData_unnormalized an object of the class 'pepData', 'proData', 'metabData',
#'   'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}},
#'   \code{\link{as.proData}}, \code{\link{as.metabData}},
#'   \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively containing
#'   the unnormalized data (but log2 transformed data), required for the metric difference 
#'   in R2m/R2c, otherwise this value can be set to the default as NULL
#'
#' @return data.frame with each corresponding to a different BECA from omicsData_beca_list
#' as well as a corresponding rank comparing that method with others
#' 
#' @examples
#' library(malbacR)
#' library(pmartR)
#' data(pmart_amide)
#' pmart_amide <- group_designation(pmart_amide, main_effects = "group",batch_id= "batch")
#' pmart_amide_log <- edata_transform(pmart_amide,"log2")
#' impObj <- imputation(pmart_amide_log)
#' pmart_amide_imp <- apply_imputation(impObj,pmart_amide_log)
#' pmart_amide_norm <- normalize_global(pmart_amide_imp, subset_fn = "all", norm_fn = "median",
#'                                     apply_norm = TRUE, backtransform = TRUE)
#' pmart_amide_imp_raw <- edata_transform(pmart_amide_imp,"abundance")
#' pmart_combat <- bc_combat(pmart_amide_norm)
#' pmart_serrf <- bc_serrf(pmart_amide_imp_raw,"group","QC","group")
#' pmart_serrf <- edata_transform(pmart_serrf,"log2")
#' pmart_power <- bc_power(pmart_amide_imp)
#' pmart_qcrfsc <- bc_qcrfsc(pmart_amide_imp_raw,qc_cname = "group",qc_val = "QC",
#'                           order_cname = "Injection_order",group_cname = "group")
#' pmart_qcrfsc <- edata_transform(pmart_qcrfsc,"log2")
#' pmart_range <- bc_range(pmart_amide_imp)
#' becas = list(ComBat = pmart_combat,SERRF = pmart_serrf,Power = pmart_power,
#'              QCRFSC = pmart_qcrfsc,Range = pmart_range)
#' r2_diff_ranking = rank_becas(omicsData_beca_list = becas,comparison_method = "r2_diff",
#'                              omicsData_unnormalized = pmart_amide_log,
#'                              main_effect_cname = "group",
#'                             batch_effect_cname = "batch")
#' @author Damon Leach
#' 
#' @export
#' 
rank_becas <- function(omicsData_beca_list,comparison_method = "r2_diff",
                       batch_effect_cname,main_effect_cname,
                       omicsData_unnormalized = NULL){
  
  # CHECKS REGARDING OMICSDATA_beca_list
  # must be a list with length of 2 or more
  if(class(omicsData_beca_list) != "list" || length(omicsData_beca_list) < 2){
    stop (paste("omicsData_beca_list must be a list with at least 2 elements"))
  }
  
  # check that list has names
  if (is.null(names(omicsData_beca_list))){
    stop("The parameter omicsData_beca_list must be a named list")
  }
  
  # check that omicsData_beca_list is of appropriate class #
  if(sum(unlist(lapply(omicsData_beca_list,function(x) !inherits(x,c("pepData", "proData", "metabData", "lipidData",
                                                                   "nmrData"))))) > 0){
    stop (paste("omicsData_beca_list must be a list containing only objects of class 'pepData', 'proData', 'metabData',",
                "'lipidData', or 'nmrData'",
                sep = ' '))
  }
  
  #  batch effect cname must be a character
  if (class(batch_effect_cname) != "character") {
    stop("Input parameter batch_effect_cname must be of class 'character'.")
  }
  if (length(batch_effect_cname) != 1) {
    stop("Input parameter batch_effect_cname must be of length 1 (e.g. vector containing a single element)")
  }
  if(sum(unlist(purrr::map(omicsData_beca_list,function(pmartObj){!any(names(pmartObj$f_data) == batch_effect_cname)}))) > 0){
    stop (paste("Input parameter batch_effect_cname must be a column found in f_data of omicsData objects in omicsData_beca_list"))
  }
  
  #  fixed effect cname must be a character
  if (class(main_effect_cname) != "character") {
    stop("Input parameter main_effect_cname must be of class 'character'.")
  }
  if (length(main_effect_cname) != 1) {
    stop("Input parameter main_effect_cname must be of length 1 (e.g. vector containing a single element)")
  }
  if(sum(unlist(purrr::map(omicsData_beca_list,function(pmartObj){!any(names(pmartObj$f_data) == main_effect_cname)}))) > 0){
    stop (paste("Input parameter main_effect_cname must be a column found in f_data of omicsData objects in omicsData_beca_list"))
  }
  
  # CHECKS REGARDING COMPARISON METHOD
  # must be a character
  if (class(comparison_method) != "character") {
    stop("Input parameter comparison_method must be of class 'character'.")
  }
  if (length(comparison_method) != 1) {
    stop("Input parameter comparison_method must be of length 1 (e.g. vector containing a single element)")
  }
  if (!comparison_method %in% c("cv","distance_pca","r2_diff")){
    stop("Input parameter comparison_method must be either 'cv','distance_pca', or 'r2_diff'.")
  }
  
  # check that group designation is stated for list of data 
  if(sum(unlist(purrr::map(omicsData_beca_list,function(pmartObj){is.null(attr(pmartObj,"group_DF"))}))) > 0){
    stop (paste("omicsData in omicsData_beca_list must have group designation applied"))
  }
  
  # check that batch designation is stated for list of data 
  if(sum(unlist(purrr::map(omicsData_beca_list,function(pmartObj){is.null(attributes(attr(pmartObj,"group_DF"))$batch_id)}))) > 0){
    stop (paste("omicsData in omicsData_beca_list must have batch designation applied"))
  }
  
  # CHECKS FOR OTHER ARGUMENTS IF COMPARISON METHOD IS R2_DIFF
  if(comparison_method == "r2_diff"){
    
    # omicsData_unnormalized cannot be NULL
    if(is.null(omicsData_unnormalized)){
      stop (paste("For the difference in R2, omicsData_unnormalized cannot be NULL"))
    }
    # check that omicsData is of appropriate class #
    if (!inherits(omicsData_unnormalized, c("pepData", "proData", "metabData", "lipidData",
                                            "nmrData"))) {
      
      stop (paste("omicsData_unnormalized must be of class 'pepData', 'proData', 'metabData',",
                  "'lipidData', or 'nmrData'",
                  sep = ' '))
    }
    if (attributes(omicsData_unnormalized)$data_info$norm_info$is_norm) {
      stop (paste("Input parameter omicsData_unnormalized must be unnormalized data"))
    }
    if (attributes(omicsData_unnormalized)$data_info$data_scale != "log2") {
      stop (paste("Input parameter omicsData_unnormalized must be on log2 scale"))
    }
    if (!any(names(omicsData_unnormalized$f_data) == batch_effect_cname)){
      stop("Input parameter batch_effect_cname must be a column found in f_data of omicsData_unnormalized")
    }
    if (!any(names(omicsData_unnormalized$f_data) == main_effect_cname)){
      stop("Input parameter main_effect_cname must be a column found in f_data of omicsData_unnormalized")
    }
    
    # check that all molecules in bc data are in unnormalized data
    numDifferent = purrr::map(omicsData_beca_list,function(pmartObj){
      edat_cnameList = pmartR::get_edata_cname(pmartObj)
      edat_cnameUnNorm = pmartR::get_edata_cname(omicsData_unnormalized)
      moleculesList = pmartObj$e_data[,which(colnames(pmartObj$e_data) == edat_cnameList)]
      moleculesUnNorm = omicsData_unnormalized$e_data[,which(colnames(omicsData_unnormalized$e_data) == edat_cnameUnNorm)]
      sum(!moleculesList %in% moleculesUnNorm)
    })
    if(sum(unlist(numDifferent)) > 0){
      stop (paste("At least one omics dataset from omicsData_beca_list has molecules that are not found in omicsData_unnormalized"))
    }
  }
  
  # if the method is cv
  if(comparison_method == "cv"){
    cv_filt_scores = purrr::map(omicsData_beca_list,function(omicsData){
      cv_obj = pmartR::cv_filter(omicsData,use_groups = TRUE)
      median(cv_obj$CV,na.rm=T)
    })
    
    ranking_df = data.frame(cv_filt_scores) %>% t() %>% data.frame() %>%
      dplyr::rename(Value = ".") %>%
      dplyr::arrange(Value) %>%
      dplyr::mutate(Value = round(Value,3)) %>%
      dplyr::mutate(Ranking = seq(from = 1, to = length(omicsData_beca_list), by = 1)) %>%
      tibble::rownames_to_column(var = "BECA")
  }
  # if the method is distance pca
  else if(comparison_method == "distance_pca"){
    distance_scores = purrr::map(omicsData_beca_list,function(omicsData){
      dim_pmart = pmartR::dim_reduction(omicsData)
      # plot(dim_pmart,omicsData = omicsData, color_by = "Batch")
      
      centroid_pmart = data.frame(SampleID = dim_pmart$SampleID,PC1 = dim_pmart$PC1,PC2 = dim_pmart$PC2) %>%
        dplyr::left_join(omicsData$f_data, by = "SampleID") %>%
        dplyr::group_by(!!as.symbol(batch_effect_cname)) %>%
        dplyr::summarise(meanPC1 = mean(PC1,na.rm=T),
                         meanPC2 = mean(PC2,na.rm=T))
      median(dist(centroid_pmart[,-1]))
    })
    ranking_df = data.frame(distance_scores) %>% t() %>% data.frame() %>%
      dplyr::rename(Value = ".") %>%
      dplyr::arrange(Value) %>%
      dplyr::mutate(Value = round(Value,3)) %>%
      dplyr::mutate(Ranking = seq(from = 1, to = length(omicsData_beca_list), by = 1)) %>%
      tibble::rownames_to_column(var = "BECA")
  }
  # if the method is r2_diff
  else if(comparison_method == "r2_diff"){
    r2_scores = purrr::map(omicsData_beca_list,function(omicsData_bc){
      edata_cname = pmartR::get_edata_cname(omicsData_unnormalized)
      fdata_cname = pmartR::get_fdata_cname(omicsData_unnormalized)
      edata_cnameNum = which(colnames(omicsData_unnormalized$e_data) == edata_cname)
      
      ball_nest <- omicsData_unnormalized$e_data %>%
        tidyr::pivot_longer(-dplyr::all_of(edata_cname),
                            names_to = fdata_cname,
                            values_to = "abundance") %>%
        dplyr::left_join(omicsData_unnormalized$f_data, by = fdata_cname) %>%
        dplyr::group_by(!!as.symbol(edata_cname)) %>%
        tidyr::nest()
      
      ball_nest <- ball_nest %>%
        dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
          #dat = ball_nest$data[[1]] %>% dplyr::filter(Batch != "B1")
          lmer_res = lme4::lmer(abundance ~ !!as.symbol(main_effect_cname) + (1|!!as.symbol(batch_effect_cname)), data = dat,
                                control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
          MuMIn::r.squaredGLMM(lmer_res)
        }))
      
      ball_r2_unnormalized <- ball_nest %>% 
        dplyr::mutate(R2m_unnorm = sapply(glmm_r2,'[',1),
                      R2c_unnorm = sapply(glmm_r2,'[',2)) %>%
        dplyr::select(!!as.symbol(edata_cname),R2m_unnorm,R2c_unnorm) %>%
        dplyr::mutate(r2_ratio_unnorm = R2m_unnorm/R2c_unnorm)
      
      # combat version
      ball_nest <- omicsData_bc$e_data %>%
        tidyr::pivot_longer(-dplyr::all_of(edata_cname),
                            names_to = fdata_cname,
                            values_to = "abundance") %>%
        dplyr::left_join(omicsData_bc$f_data, by = fdata_cname) %>%
        dplyr::group_by(!!as.symbol(edata_cname)) %>%
        tidyr::nest()
      
      ball_nest <- ball_nest %>%
        dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
          lmer_res = lme4::lmer(abundance ~ !!as.symbol(main_effect_cname) + (1|!!as.symbol(batch_effect_cname)), data = dat,
                                control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
          MuMIn::r.squaredGLMM(lmer_res)
        }))
      
      ball_r2_bc <- ball_nest %>% 
        dplyr::mutate(R2m = sapply(glmm_r2,'[',1),
                      R2c = sapply(glmm_r2,'[',2)) %>%
        dplyr::select(!!as.symbol(edata_cname),R2m,R2c) %>%
        dplyr::mutate(r2_ratio = R2m/R2c)
      
      ball_all <- ball_r2_bc %>%
        dplyr::left_join(ball_r2_unnormalized,by = edata_cname)
      
      abund_df = omicsData_unnormalized$e_data %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = edata_cname) %>%
        as.matrix() %>%
        matrixStats::rowMedians(na.rm=T) %>%
        data.frame() %>%
        tibble::rownames_to_column(var = edata_cname) %>%
        dplyr::rename("medAbund" = ".")
      
      quartile_info = quantile(abund_df$medAbund)
      abund_df <- abund_df %>%
        dplyr::mutate(quantileNum = ifelse(medAbund < quartile_info[2],"Q1",
                                           ifelse(medAbund > quartile_info[4],"Q4",
                                                  ifelse(medAbund > quartile_info[2] & medAbund < 
                                                           quartile_info[3],"Q2","Q3")))) %>%
        dplyr::filter(quantileNum %in% c("Q3","Q4"))
      
      abundant_molecules = abund_df[,which(colnames(abund_df) == edata_cname)]
      
      ball_all_med <- ball_all %>%
        dplyr::filter(!!as.symbol(edata_cname) %in% abundant_molecules) %>%
        dplyr::mutate(r2_ratio_diff = r2_ratio - r2_ratio_unnorm) %>%
        dplyr::ungroup() %>%
        dplyr::summarise(medDiff = median(r2_ratio_diff))
      
      ball_all_med$medDiff
    })
    ranking_df = data.frame(r2_scores) %>% t() %>% data.frame() %>%
      dplyr::rename(Value = ".") %>%
      dplyr::arrange(desc(Value)) %>%
      dplyr::mutate(Ranking = seq(from = 1, to = length(omicsData_beca_list), by = 1)) %>%
      tibble::rownames_to_column(var = "BECA")
    # update for ties
    if(length(ranking_df$Value) != length(unique(ranking_df$Value))){
      for(i in 2:nrow(ranking_df)){
        if(ranking_df$Value[i] == ranking_df$Value[i-1]){
          ranking_df$Ranking[i] = ranking_df$Ranking[i-1]
        }
      }
    }
    ranking_df <- ranking_df %>% dplyr::mutate(Value = round(Value,3))
  }
  ranking_df
}
