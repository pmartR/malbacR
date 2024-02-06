context('Rank BECAs')

test_that('bc_tiger returns the correct data frame and attributes',{
  
  # Load the reduced peptide data frames ---------------------------------------
  load(system.file('testdata',
                   'mini_udn.Rdata',
                   package = 'malbacR'))
  
  # Run as.pepData with agreeable data frames ----------------------------------
  
  # Construct a pepData object with the edata, fdata, and emeta data frames.
  mdata <- pmartR::as.metabData(e_data = edata,
                                f_data = fdata,
                                e_meta = emeta,
                                edata_cname = 'Metabolite',
                                fdata_cname = 'SampleID',
                                emeta_cname = 'Metabolite',
                                data_scale = "log2")
  
  # run batch correction on five different methods
  mdata <- pmartR::group_designation(mdata,main_effects = "Sex",batch_id = "BatchName")
  molfilt <- pmartR::molecule_filter(mdata,use_groups = TRUE, use_batch = TRUE)
  mdata <- pmartR::applyFilt(molfilt, mdata)
  mdata_norm <- pmartR::normalize_global(mdata,subset_fn = "all", norm_fn = "median",
                                         apply_norm = TRUE, backtransform = TRUE)
  # imputation
  impObj <- imputation(mdata)
  mdata_imp <- apply_imputation(impObj,mdata)
  mdata_imp_raw <- pmartR::edata_transform(mdata_imp,"abundance")
  
  # run batch correction based on standard procedure for each BECA
  # combat
  mdata_combat <- bc_combat(mdata_norm)
  # eigenms
  mdata_eigen <- bc_eigenMS(mdata)
  # serrf 
  mdata_serrf <- bc_serrf(mdata_imp_raw, sampletype_cname = "Sex",test_val = "QC.NIST","Sex")
  mdata_serrf <- pmartR::edata_transform(mdata_serrf, "log2")
  # qcrfsc
  mdata_qcrfsc <- bc_qcrfsc(mdata_imp_raw, qc_cname = "Sex", qc_val = "QC.NIST",
                            order_cname = "RunOrderOverall",group_cname = "Sex")
  mdata_qcrfsc <- pmartR::edata_transform(mdata_qcrfsc,"log2")
  # power
  mdata_power <- bc_power(mdata)
  
  # Run through the potential error messages -----------------------------------
  
  # omicsData_beca_list needs to be a list with at least two elements
  expect_error(rank_becas(omicsData_beca_list = mdata_combat,comparison_method = "cv"),
               "omicsData_beca_list must be a list with at least 2 elements")
  expect_error(rank_becas(omicsData_beca_list = list(mdata_combat),comparison_method = "cv"),
               "omicsData_beca_list must be a list with at least 2 elements")
  # must be a named list
  expect_error(rank_becas(omicsData_beca_list = list(mdata_combat,mdata_eigen),comparison_method = "cv"),
               "The parameter omicsData_beca_list must be a named list")
  # objects must all be omicsData
  expect_error(rank_becas(omicsData_beca_list = list(ComBat = mdata_combat$e_data,EigenMS = mdata_eigen),comparison_method = "cv"),
               "omicsData_beca_list must be a list containing only objects of class")
  becas = list(ComBat = mdata_combat, EigenMS = mdata_eigen,
               SERRF = mdata_serrf, QCRFSC = mdata_qcrfsc, Power = mdata_power)

  # batch effect cname must be character
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method ="cv",batch_effect_cname = 2),
               "Input parameter batch_effect_cname must be of class 'character'")
  # batch effect cname must be character in fdata columns
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method ="cv",batch_effect_cname = "batch"),
               "Input parameter batch_effect_cname must be a column found in f_data")
  # batch effect cname must be character in fdata columns length of 1
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method ="cv",batch_effect_cname = c("Batch","BatchNum")),
               "Input parameter batch_effect_cname must be of length 1")
  
  # main effect cname must be character
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method ="cv",batch_effect_cname = "BatchNum",
                          main_effect_cname = 2),
               "Input parameter main_effect_cname must be of class 'character'")
  # main effect cname must be character in fdata columns
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method ="cv",batch_effect_cname = "BatchNum",
                          main_effect_cname = "group"),
               "Input parameter main_effect_cname must be a column found in f_data")
  # main effect cname must be character in fdata columns length of 1
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method ="cv",batch_effect_cname = "BatchNum",
                          main_effect_cname = c("Sex","Age")),
               "Input parameter main_effect_cname must be of length 1")
  
  # comparison method must be character
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = 2,batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex"),
               "Input parameter comparison_method must be of class 'character'")
  # comparison method must be one of three options cv, distance_pca, r2_diff
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "PCA",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex"),
               "Input parameter comparison_method must be either 'cv'")
  # comparison method must be of length 1
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = c("cv","r2_diff"),batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex"),
               "Input parameter comparison_method must be of length 1")
  
  # extra checks if comparison method is r2_diff
  # we need unnormalized data
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex"),
               "For the difference in R2, omicsData_unnormalized cannot be NULL")
  # unnormalized data must be omics object
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex", omicsData_unnormalized = mdata$e_data),
               "omicsData_unnormalized must be of class")
  # data cannot be normalized
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex", omicsData_unnormalized = mdata_norm),
               "Input parameter omicsData_unnormalized must be unnormalized data")
  # unnormalized must be log2 transformed
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex", omicsData_unnormalized = pmartR::edata_transform(mdata,"abundance")),
               "Input parameter omicsData_unnormalized must be on log2 scale")
  # no main effect or batch effect in unnormalized is bad
  mdata_noInfo = mdata
  mdata_noInfo$f_data <- mdata_noInfo$f_data[,-c(2,6)]
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex", omicsData_unnormalized = mdata_noInfo),
               "Input parameter batch_effect_cname must be a column found in f_data of omicsData_unnormalized")
  mdata_noInfo = mdata
  mdata_noInfo$f_data <- mdata_noInfo$f_data[,-c(6)]
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex", omicsData_unnormalized = mdata_noInfo),
               "Input parameter main_effect_cname must be a column found in f_data of omicsData_unnormalized")
  
  # every molecule in omics list should be in unnormalized data (does not have to be the other way per se)
  unknown = mdata$e_data$Metabolite[which(stringr::str_detect(mdata$e_data$Metabolite,"Unidentified"))]
  molfilt <- pmartR::custom_filter(mdata, e_data_remove = unknown)
  mdata_reduced <- pmartR::applyFilt(molfilt,mdata)
  expect_error(rank_becas(omicsData_beca_list = becas,
                          comparison_method = "r2_diff",batch_effect_cname = "BatchNum",
                          main_effect_cname = "Sex", omicsData_unnormalized = mdata_reduced),
               "At least one omics dataset from omicsData_beca_list has molecules that are not found in")

  # Check the dimensions of results --------------------------------------------
  # coefficient of variation
  # function version
  cv_rankings = rank_becas(omicsData_beca_list = becas,
                            comparison_method = "cv",batch_effect_cname = "BatchNum",
                            main_effect_cname = "Sex")
  # manual version
  cv_combat = median(pmartR::cv_filter(mdata_combat,use_groups = T)$CV,na.rm=T)
  cv_eigen = median(pmartR::cv_filter(mdata_eigen,use_groups = T)$CV,na.rm=T)
  cv_serrf = median(pmartR::cv_filter(mdata_serrf,use_groups = T)$CV,na.rm=T)
  cv_qcrfsc = median(pmartR::cv_filter(mdata_qcrfsc,use_groups = T)$CV,na.rm=T)
  cv_power = median(pmartR::cv_filter(mdata_power,use_groups = T)$CV,na.rm=T)
  cv_manual_rankings = data.frame(BECA = c("ComBat","EigenMS","SERRF","QCRFSC","Power"),
             Value = c(cv_combat,cv_eigen,cv_serrf,cv_qcrfsc,cv_power)) %>%
    dplyr::arrange(Value) %>%
    dplyr::mutate(Value = round(Value,3)) %>%
    dplyr::mutate(Ranking = seq(from = 1, to = 5, by = 1))
  # compare versions
  expect_equal(cv_rankings, cv_manual_rankings)
  
  # distance_pca
  pca_rankings = rank_becas(omicsData_beca_list = becas,
                            comparison_method = "distance_pca", batch_effect_cname = "BatchNum",
                            main_effect_cname = "Sex")
  # manual version
  dim_combat = pmartR::dim_reduction(mdata_combat)
  centroid_combat = data.frame(SampleID = dim_combat$SampleID, PC1 = dim_combat$PC1, PC2 = dim_combat$PC2) %>%
    dplyr::left_join(mdata_combat$f_data, by = "SampleID") %>%
    dplyr::group_by(BatchNum) %>%
    dplyr::summarise(meanPC1 = mean(PC1, na.rm = T),
                     meanPC2 = mean(PC2, na.rm = T))
  dist_combat = as.numeric(sqrt((centroid_combat[1,2] - centroid_combat[2,2])^2 + (centroid_combat[1,3] - centroid_combat[2,3])^2))
  dim_eigen = pmartR::dim_reduction(mdata_eigen)
  centroid_eigen = data.frame(SampleID = dim_eigen$SampleID, PC1 = dim_eigen$PC1, PC2 = dim_eigen$PC2) %>%
    dplyr::left_join(mdata_eigen$f_data, by = "SampleID") %>%
    dplyr::group_by(BatchNum) %>%
    dplyr::summarise(meanPC1 = mean(PC1, na.rm = T),
                     meanPC2 = mean(PC2, na.rm = T))
  dist_eigen = as.numeric(sqrt((centroid_eigen[1,2] - centroid_eigen[2,2])^2 + (centroid_eigen[1,3] - centroid_eigen[2,3])^2))
  dim_qcrfsc = pmartR::dim_reduction(mdata_qcrfsc)
  centroid_qcrfsc = data.frame(SampleID = dim_qcrfsc$SampleID, PC1 = dim_qcrfsc$PC1, PC2 = dim_qcrfsc$PC2) %>%
    dplyr::left_join(mdata_qcrfsc$f_data, by = "SampleID") %>%
    dplyr::group_by(BatchNum) %>%
    dplyr::summarise(meanPC1 = mean(PC1, na.rm = T),
                     meanPC2 = mean(PC2, na.rm = T))
  dist_qcrfsc = as.numeric(sqrt((centroid_qcrfsc[1,2] - centroid_qcrfsc[2,2])^2 + (centroid_qcrfsc[1,3] - centroid_qcrfsc[2,3])^2))
  dim_serrf = pmartR::dim_reduction(mdata_serrf)
  centroid_serrf = data.frame(SampleID = dim_serrf$SampleID, PC1 = dim_serrf$PC1, PC2 = dim_serrf$PC2) %>%
    dplyr::left_join(mdata_serrf$f_data, by = "SampleID") %>%
    dplyr::group_by(BatchNum) %>%
    dplyr::summarise(meanPC1 = mean(PC1, na.rm = T),
                     meanPC2 = mean(PC2, na.rm = T))
  dist_serrf = as.numeric(sqrt((centroid_serrf[1,2] - centroid_serrf[2,2])^2 + (centroid_serrf[1,3] - centroid_serrf[2,3])^2))
  dim_power = pmartR::dim_reduction(mdata_power)
  centroid_power = data.frame(SampleID = dim_power$SampleID, PC1 = dim_power$PC1, PC2 = dim_power$PC2) %>%
    dplyr::left_join(mdata_power$f_data, by = "SampleID") %>%
    dplyr::group_by(BatchNum) %>%
    dplyr::summarise(meanPC1 = mean(PC1, na.rm = T),
                     meanPC2 = mean(PC2, na.rm = T))
  dist_power = as.numeric(sqrt((centroid_power[1,2] - centroid_power[2,2])^2 + (centroid_power[1,3] - centroid_power[2,3])^2))

  pca_manual_rankings = data.frame(BECA = c("ComBat","EigenMS","SERRF","QCRFSC","Power"),
                                   Value = c(dist_combat,dist_eigen,dist_serrf,dist_qcrfsc,dist_power)) %>%
    dplyr::arrange(Value) %>%
    dplyr::mutate(Value = round(Value,3)) %>%
    dplyr::mutate(Ranking = seq(from = 1, to = 5, by = 1))
  expect_equal(pca_rankings,pca_manual_rankings)
  
  # r2_diff
  r2_rankings = rank_becas(omicsData_beca_list = becas,
                                          comparison_method = "r2_diff", batch_effect_cname = "BatchNum",
                                          main_effect_cname = "Sex", omicsData_unnormalized = mdata)
  # manual version
  # run the data on unnormalized data
  unnorm_nest <- mdata$e_data %>%
    tidyr::pivot_longer(-Metabolite,
                        names_to = "SampleID",
                        values_to = "Abundance") %>%
    dplyr::left_join(mdata$f_data, by = "SampleID") %>%
    dplyr::group_by(Metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
      lmer_res = lme4::lmer(Abundance ~ Sex + (1|BatchNum), data = dat,
                            control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
      MuMIn::r.squaredGLMM(lmer_res)
    })) %>% 
    dplyr::mutate(R2m_unnorm = sapply(glmm_r2,'[',1),
                  R2c_unnorm = sapply(glmm_r2,'[',2)) %>%
    dplyr::select(Metabolite,R2m_unnorm,R2c_unnorm) %>%
    dplyr::mutate(r2_ratio_unnorm = R2m_unnorm/R2c_unnorm)
    
  # find the 50% most abundant molecules
  abund_df = mdata$e_data %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Metabolite") %>%
    as.matrix() %>%
    matrixStats::rowMedians(na.rm=T) %>%
    data.frame() %>%
    tibble::rownames_to_column(var = "Metabolite") %>%
    dplyr::rename("medAbund" = ".")
  quartile_info = quantile(abund_df$medAbund)
  abund_df <- abund_df %>%
    dplyr::mutate(quantileNum = ifelse(medAbund < quartile_info[2],"Q1",
                                       ifelse(medAbund > quartile_info[4],"Q4",
                                              ifelse(medAbund > quartile_info[2] & medAbund < 
                                                       quartile_info[3],"Q2","Q3")))) %>%
    dplyr::filter(quantileNum %in% c("Q3","Q4"))
  abundant_molecules = abund_df[,which(colnames(abund_df) == "Metabolite")]
  
  # run the data on BC combat version
  combat_nest <- mdata_combat$e_data %>%
    tidyr::pivot_longer(-Metabolite,
                        names_to = "SampleID",
                        values_to = "Abundance") %>%
    dplyr::left_join(mdata_combat$f_data, by = "SampleID") %>%
    dplyr::group_by(Metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
      lmer_res = lme4::lmer(Abundance ~ Sex + (1|BatchNum), data = dat,
                            control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
      MuMIn::r.squaredGLMM(lmer_res)
    })) %>% 
    dplyr::mutate(R2m = sapply(glmm_r2,'[',1),
                  R2c = sapply(glmm_r2,'[',2)) %>%
    dplyr::select(Metabolite,R2m,R2c) %>%
    dplyr::mutate(r2_ratio = R2m/R2c)
  
  combat_all <- combat_nest %>%
    dplyr::left_join(unnorm_nest,by = "Metabolite")
  combat_all_med <- combat_all %>%
    dplyr::filter(Metabolite %in% abundant_molecules) %>%
    dplyr::mutate(r2_ratio_diff = r2_ratio - r2_ratio_unnorm) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(medDiff = median(r2_ratio_diff))
  combat_diff = combat_all_med$medDiff
  
  eigen_nest <- mdata_eigen$e_data %>%
    tidyr::pivot_longer(-Metabolite,
                        names_to = "SampleID",
                        values_to = "Abundance") %>%
    dplyr::left_join(mdata_eigen$f_data, by = "SampleID") %>%
    dplyr::group_by(Metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
      lmer_res = lme4::lmer(Abundance ~ Sex + (1|BatchNum), data = dat,
                            control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
      MuMIn::r.squaredGLMM(lmer_res)
    })) %>% 
    dplyr::mutate(R2m = sapply(glmm_r2,'[',1),
                  R2c = sapply(glmm_r2,'[',2)) %>%
    dplyr::select(Metabolite,R2m,R2c) %>%
    dplyr::mutate(r2_ratio = R2m/R2c)
  
  eigen_all <- eigen_nest %>%
    dplyr::left_join(unnorm_nest,by = "Metabolite")
  eigen_all_med <- eigen_all %>%
    dplyr::filter(Metabolite %in% abundant_molecules) %>%
    dplyr::mutate(r2_ratio_diff = r2_ratio - r2_ratio_unnorm) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(medDiff = median(r2_ratio_diff))
  eigen_diff = eigen_all_med$medDiff
  
  power_nest <- mdata_power$e_data %>%
    tidyr::pivot_longer(-Metabolite,
                        names_to = "SampleID",
                        values_to = "Abundance") %>%
    dplyr::left_join(mdata_power$f_data, by = "SampleID") %>%
    dplyr::group_by(Metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
      lmer_res = lme4::lmer(Abundance ~ Sex + (1|BatchNum), data = dat,
                            control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
      MuMIn::r.squaredGLMM(lmer_res)
    })) %>% 
    dplyr::mutate(R2m = sapply(glmm_r2,'[',1),
                  R2c = sapply(glmm_r2,'[',2)) %>%
    dplyr::select(Metabolite,R2m,R2c) %>%
    dplyr::mutate(r2_ratio = R2m/R2c)
  
  power_all <- power_nest %>%
    dplyr::left_join(unnorm_nest,by = "Metabolite")
  power_all_med <- power_all %>%
    dplyr::filter(Metabolite %in% abundant_molecules) %>%
    dplyr::mutate(r2_ratio_diff = r2_ratio - r2_ratio_unnorm) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(medDiff = median(r2_ratio_diff))
  power_diff = power_all_med$medDiff
  
  serrf_nest <- mdata_serrf$e_data %>%
    tidyr::pivot_longer(-Metabolite,
                        names_to = "SampleID",
                        values_to = "Abundance") %>%
    dplyr::left_join(mdata_serrf$f_data, by = "SampleID") %>%
    dplyr::group_by(Metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
      lmer_res = lme4::lmer(Abundance ~ Sex + (1|BatchNum), data = dat,
                            control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
      MuMIn::r.squaredGLMM(lmer_res)
    })) %>% 
    dplyr::mutate(R2m = sapply(glmm_r2,'[',1),
                  R2c = sapply(glmm_r2,'[',2)) %>%
    dplyr::select(Metabolite,R2m,R2c) %>%
    dplyr::mutate(r2_ratio = R2m/R2c)
  
  serrf_all <- serrf_nest %>%
    dplyr::left_join(unnorm_nest,by = "Metabolite")
  serrf_all_med <- serrf_all %>%
    dplyr::filter(Metabolite %in% abundant_molecules) %>%
    dplyr::mutate(r2_ratio_diff = r2_ratio - r2_ratio_unnorm) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(medDiff = median(r2_ratio_diff))
  serrf_diff = serrf_all_med$medDiff
  
  qcrfsc_nest <- mdata_qcrfsc$e_data %>%
    tidyr::pivot_longer(-Metabolite,
                        names_to = "SampleID",
                        values_to = "Abundance") %>%
    dplyr::left_join(mdata_qcrfsc$f_data, by = "SampleID") %>%
    dplyr::group_by(Metabolite) %>%
    tidyr::nest() %>%
    dplyr::mutate(glmm_r2 = purrr::map(data,function(dat){
      lmer_res = lme4::lmer(Abundance ~ Sex + (1|BatchNum), data = dat,
                            control = lme4::lmerControl(check.conv.singular = lme4::.makeCC(action = "ignore", tol = 1e-4)))
      MuMIn::r.squaredGLMM(lmer_res)
    })) %>% 
    dplyr::mutate(R2m = sapply(glmm_r2,'[',1),
                  R2c = sapply(glmm_r2,'[',2)) %>%
    dplyr::select(Metabolite,R2m,R2c) %>%
    dplyr::mutate(r2_ratio = R2m/R2c)
  
  qcrfsc_all <- qcrfsc_nest %>%
    dplyr::left_join(unnorm_nest,by = "Metabolite")
  qcrfsc_all_med <- qcrfsc_all %>%
    dplyr::filter(Metabolite %in% abundant_molecules) %>%
    dplyr::mutate(r2_ratio_diff = r2_ratio - r2_ratio_unnorm) %>%
    dplyr::ungroup() %>%
    dplyr::summarise(medDiff = median(r2_ratio_diff))
  qcrfsc_diff = qcrfsc_all_med$medDiff
  
  r2_manual_rankings = data.frame(BECA = c("ComBat","EigenMS","SERRF","QCRFSC","Power"),
                                   Value = c(combat_diff,eigen_diff,serrf_diff,qcrfsc_diff,power_diff)) %>%
    dplyr::arrange(Value) %>%
    dplyr::mutate(Value = round(Value,3)) %>%
    dplyr::mutate(Ranking = seq(from = 1, to = 5, by = 1))
  for(i in 2:nrow(r2_manual_rankings)){
    if(r2_manual_rankings$Value[i] == r2_manual_rankings$Value[i-1]){
      r2_manual_rankings$Ranking[i] = r2_manual_rankings$Ranking[i-1]
      }
  }
  expect_equal(r2_rankings,r2_manual_rankings)
  
})
