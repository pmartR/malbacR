## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE)

## -----------------------------------------------------------------------------
library(malbacR)
library(pmartR)
library(ggplot2)

## -----------------------------------------------------------------------------
data(pmart_amide)

## -----------------------------------------------------------------------------
pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
pmart_amide_log <- edata_transform(pmart_amide,"log2")
pmart_amide_norm <- normalize_global(pmart_amide_log,subset_fn = "all",norm_fn = "median",
                                apply_norm = TRUE,backtransform = TRUE)

## -----------------------------------------------------------------------------
# range scaling
amide_range <- bc_range(omicsData = pmart_amide_log)

# power scaling
amide_power <- bc_power(omicsData = pmart_amide_log)

# pareto scaling
amide_pareto <- bc_pareto(omicsData = pmart_amide_log)

## -----------------------------------------------------------------------------
tigerFilt <- tiger_filter(pmart_amide,sampletype_cname = "group",test_val = "QC")
pmart_amideFilt <- apply_tigerFilt(tigerFilt,pmart_amide)
amide_tiger_abundance <- bc_tiger(omicsData = pmart_amideFilt,sampletype_cname = "group",
                        test_val = "QC",injection_cname = "Injection_order",group_cname = "group")
amide_tiger <- edata_transform(omicsData = amide_tiger_abundance, data_scale = "log2")

## -----------------------------------------------------------------------------
amide_qcrlsc <- bc_qcrlsc(omicsData = pmart_amide_log,block_cname = "batch",
                          qc_cname = "group", qc_val = "QC", order_cname = "Injection_order",
                          missing_thresh = 0.5, rsd_thresh = 0.3, backtransform  = FALSE,keep_qc = FALSE)

## -----------------------------------------------------------------------------
impObj <- imputation(omicsData = pmart_amide_log)
amide_imp_log <- apply_imputation(impObj,pmart_amide_log)
amide_imp_raw <- edata_transform(amide_imp_log,"abundance")

## -----------------------------------------------------------------------------
amide_serrf_abundance <- bc_serrf(omicsData = amide_imp_raw,sampletype_cname = "group",test_val = "QC",group_cname = "group")
amide_serrf <- edata_transform(omicsData = amide_serrf_abundance,data_scale = "log2")

## -----------------------------------------------------------------------------
amide_qcrfsc_abundance <- bc_qcrfsc(omicsData = amide_imp_raw,qc_cname = "group",qc_val = "QC",
                                    order_cname = "Injection_order",group_cname = "group",
                                    ntree = 500, keep_qc = FALSE)
amide_qcrfsc <- edata_transform(omicsData = amide_qcrfsc_abundance,data_scale = "log2")

## -----------------------------------------------------------------------------
# combat batch correction
amide_combat <- bc_combat(omicsData = pmart_amide_norm, use_groups = FALSE)

## -----------------------------------------------------------------------------
amide_eigen <- bc_eigenMS(omicsData = pmart_amide_log)

## ----results = 'hide'---------------------------------------------------------
amide_wave_abundance <- bc_waveica(omicsData = amide_imp_raw, batch_cname = "batch",
                         version = "WaveICA2.0",
                         injection_cname = "Injection_order",alpha = 0,
                         cutoff_injection = 0.1, K = 10,
                         negative_to_na = TRUE)
amide_wave <- edata_transform(omicsData = amide_wave_abundance, data_scale = "log2")

## ----out.width = "33%"--------------------------------------------------------
p1 <- plot(dim_reduction(omicsData = pmart_amide_log),omicsData = pmart_amide,color_by = "batch") + labs(title = "Amide: Unadjusted")
p2 <- plot(dim_reduction(omicsData = amide_range),omicsData = amide_range,color_by = "batch") + labs(title = "Amide: Range")
p3 <- plot(dim_reduction(omicsData = amide_power),omicsData = amide_power,color_by = "batch") + labs(title = "Amide: Power")
p4 <- plot(dim_reduction(omicsData = amide_pareto),omicsData = amide_pareto,color_by = "batch") + labs(title = "Amide: Pareto")
p5 <- plot(dim_reduction(omicsData = amide_tiger),omicsData = amide_tiger,color_by = "batch") + labs(title = "Amide: TIGER")
p6 <- plot(dim_reduction(omicsData = amide_qcrlsc),omicsData = amide_qcrlsc,color_by = "batch") + labs(title = "Amide: QC-RLSC")
p7 <- plot(dim_reduction(omicsData = amide_combat),omicsData = amide_combat,color_by = "batch") + labs(title = "Amide: ComBat")
p8 <- plot(dim_reduction(omicsData = amide_wave),omicsData = amide_wave,color_by = "batch") + labs(title = "Amide: WaveICA2.0")
p9 <- plot(dim_reduction(omicsData = amide_serrf),omicsData = amide_serrf,color_by = "batch") + labs(title = "Amide: SERRF")
p10 <- plot(dim_reduction(omicsData = amide_qcrfsc),omicsData = amide_qcrfsc,color_by = "batch") + labs(title = "Amide: QC-RFSC")
p11 <- plot(dim_reduction(omicsData = amide_eigen),omicsData = amide_eigen,color_by = "batch") + labs(title = "Amide: EigenMS")

p1;p2;p3;p4;p5;p6;p7;p8;p9;p10;p11

## ----include = FALSE, eval = FALSE--------------------------------------------
#  # unadjusted
#  pmart_amide_res <- dim_reduction(omicsData = pmart_amide_log)
#  amide1 <- data.frame(sample = pmart_amide_res$SampleID, Batch = attr(pmart_amide_res,"group_DF")$Group,
#             pc1 = pmart_amide_res$PC1, pc2 = pmart_amide_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(pmart_amide_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(pmart_amide_res,"R2")[2],3),")"),
#         title = "Amide: Unadjusted")
#  # combat
#  amide_combat_res <- dim_reduction(omicsData = amide_combat)
#  amide2 <- data.frame(sample = amide_combat_res$SampleID, Batch = attr(amide_combat_res,"group_DF")$Group,
#             pc1 = amide_combat_res$PC1, pc2 = amide_combat_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(amide_combat_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(amide_combat_res,"R2")[2],3),")"),
#         title = "Amide: ComBat")
#  # waveica
#  amide_wave_res <- dim_reduction(omicsData = amide_wave)
#  amide3 <- data.frame(sample = amide_wave_res$SampleID, Batch = attr(amide_wave_res,"group_DF")$Group,
#             pc1 = amide_wave_res$PC1, pc2 = amide_wave_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.50,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(amide_wave_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(amide_wave_res,"R2")[2],3),")"),
#         title = "Amide: WaveICA2.0")
#  # serrf
#  amide_serrf_res <- dim_reduction(omicsData = amide_serrf)
#  amide4 <- data.frame(sample = amide_serrf_res$SampleID, Batch = attr(amide_serrf_res,"group_DF")$Group,
#             pc1 = amide_serrf_res$PC1, pc2 = amide_serrf_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(amide_serrf_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(amide_serrf_res,"R2")[2],3),")"),
#         title = "Amide: SERRF")
#  
#  (amide1|amide2)/(amide3|amide4)

## -----------------------------------------------------------------------------
data(pmart_mix)

## -----------------------------------------------------------------------------
pmart_mix <- group_designation(pmart_mix,main_effects = "BatchNum",batch_id = "BatchNum")
pmart_mix_log <- edata_transform(pmart_mix,"log2")
pmart_mix_norm <- normalize_global(pmart_mix_log,subset_fn = "all",norm_fn = "median",
                                apply_norm = TRUE,backtransform = TRUE)

## -----------------------------------------------------------------------------
mix_ruv <- bc_ruvRandom(omicsData = pmart_mix_log, nc_cname = "tag",nc_val = "IS", k = 3)

## -----------------------------------------------------------------------------
mix_nomis_abundance <- bc_nomis(omicsData = pmart_mix, is_cname = "tag", is_val = "IS", num_pc = 2)
mix_nomis <- edata_transform(omicsData = mix_nomis_abundance, data_scale = "log2")

## -----------------------------------------------------------------------------
mix_combat <- bc_combat(omicsData = pmart_mix_norm)
mix_range <- bc_range(omicsData = pmart_mix_log)
mix_pareto <- bc_pareto(omicsData = pmart_mix_log)
mix_power <- bc_power(omicsData = pmart_mix_log)
mix_waveica_abundance <- bc_waveica(omicsData = pmart_mix,batch_cname = "BatchNum",
                          version = "WaveICA",cutoff_batch = 0.05, cutoff_group = 0.05,
                          negative_to_na = TRUE)
mix_waveica <- edata_transform(omicsData = mix_waveica_abundance, data_scale = "log2")

## ----out.width = "33%"--------------------------------------------------------
p1 <- plot(dim_reduction(omicsData = pmart_mix_log),omicsData = pmart_mixFilt,color_by = "BatchNum") + labs(title = "Mix: Unadjusted")
p2 <- plot(dim_reduction(omicsData = mix_ruv),omicsData = mix_ruv,color_by = "BatchNum") + labs(title = "Mix: ruv-Random")
p3 <- plot(dim_reduction(omicsData = mix_nomis),omicsData = mix_nomis,color_by = "BatchNum") + labs(title = "Mix: NOMIS")
p4 <- plot(dim_reduction(omicsData = mix_combat),omicsData = mix_combat,color_by = "BatchNum") + labs(title = "Mix: ComBat")
p5 <- plot(dim_reduction(omicsData = mix_range),omicsData = mix_range,color_by = "BatchNum") + labs(title = "Mix: Range")
p6 <- plot(dim_reduction(omicsData = mix_power),omicsData = mix_power,color_by = "BatchNum") + labs(title = "Mix: Power")
p7 <- plot(dim_reduction(omicsData = mix_pareto),omicsData = mix_pareto,color_by = "BatchNum") + labs(title = "Mix: Pareto")
p8 <- plot(dim_reduction(omicsData = mix_waveica),omicsData = mix_pareto,color_by = "BatchNum") + labs(title = "Mix: WaveICA")

p1;p2;p3;p4;p5;p6;p7;p8

## ----include = FALSE, eval = FALSE--------------------------------------------
#  # unadjusted
#  pmart_mix_res <- dim_reduction(omicsData = pmart_mix)
#  mix1 <- data.frame(sample = pmart_mix_res$SampleID, Batch = attr(pmart_mix_res,"group_DF")$Group,
#             pc1 = pmart_mix_res$PC1, pc2 = pmart_mix_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(pmart_mix_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(pmart_mix_res,"R2")[2],3),")"),
#         title = "Mix: Unadjusted")
#  # ruv random
#  mix_ruv_res <- dim_reduction(omicsData = mix_ruv)
#  mix2 <- data.frame(sample = mix_ruv_res$SampleID, Batch = attr(mix_ruv_res,"group_DF")$Group,
#             pc1 = mix_ruv_res$PC1, pc2 = mix_ruv_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(mix_ruv_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(mix_ruv_res,"R2")[2],3),")"),
#         title = "Mix: RUV-random")
#  # nomis
#  mix_nomis_res <- dim_reduction(omicsData = mix_nomis)
#  mix3 <- data.frame(sample = mix_nomis_res$SampleID, Batch = attr(mix_nomis_res,"group_DF")$Group,
#             pc1 = mix_nomis_res$PC1, pc2 = mix_nomis_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(mix_nomis_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(mix_nomis_res,"R2")[2],3),")"),
#         title = "Mix: NOMIS")
#  # combat
#  mix_combat_res <- dim_reduction(omicsData = mix_combat)
#  mix4 <- data.frame(sample = mix_combat_res$SampleID, Batch = attr(mix_combat_res,"group_DF")$Group,
#             pc1 = mix_combat_res$PC1, pc2 = mix_combat_res$PC2) %>%
#    ggplot(aes(x = pc1, y = pc2,color = Batch)) +
#    geom_point(alpha = 0.5,size = 4) +
#    theme_bw() +
#    labs(x = paste0("PC1 (R^2 = ",round(attr(mix_combat_res,"R2")[1],3),")"),
#         y = paste0("PC1 (R^2 = ",round(attr(mix_combat_res,"R2")[2],3),")"),
#         title = "Mix: ComBat")
#  
#  (mix1|mix2)/(mix3|mix4)

