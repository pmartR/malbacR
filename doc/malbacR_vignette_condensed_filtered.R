## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE)

## -----------------------------------------------------------------------------
library(malbacR)
library(pmartR)
library(ggplot2)

## -----------------------------------------------------------------------------
data(pmart_amideFilt)

## ----results = 'hide'---------------------------------------------------------
# SCALING METHODS
# range scaling
amide_range <- bc_range(omicsData = pmart_amideFilt)

# power scaling
amide_power <- bc_power(omicsData = pmart_amideFilt)

# pareto scaling
amide_pareto <- bc_pareto(omicsData = pmart_amideFilt)

# QUALITY CONTROL METHODS
# TIGER
pmart_amideFilt_abundance <- edata_transform(omicsData = pmart_amideFilt, data_scale = "abundance")
amide_tiger_abundance <- bc_tiger(omicsData = pmart_amideFilt_abundance,sampletype_cname = "group",
                        test_val = "QC",injection_cname = "Injection_order",group_cname = "group")
amide_tiger <- edata_transform(omicsData = amide_tiger_abundance, data_scale = "log2")

# QC-RLSC
amide_qcrlsc <- bc_qcrlsc(omicsData = pmart_amideFilt,block_cname = "batch",
                          qc_cname = "group", qc_val = "QC", order_cname = "Injection_order",
                          missing_thresh = 0.5, rsd_thresh = 0.3, backtransform  = FALSE,keep_qc = FALSE)

# SERRF
amide_serrf_abundance <- bc_serrf(omicsData = pmart_amideFilt_abundance,sampletype_cname = "group",test_val = "QC",group_cname = "group")
amide_serrf <- edata_transform(omicsData = amide_serrf_abundance, data_scale = "log2")

# OTHER METHODS
# ComBat
amide_combat <- bc_combat(omicsData = pmart_amideFilt, use_groups = FALSE)

# EigenMS
amide_eigen <- bc_eigenMS(omicsData = pmart_amideFilt)

# WaveICA2.0
amide_wave_abundance <- bc_waveica(omicsData = pmart_amideFilt_abundance, injection_cname = "Injection_order",
                         version = "WaveICA2.0", alpha = 0, cutoff_injection = 0.1, K = 10,
                         negative_to_na = TRUE)
amide_wave <- edata_transform(omicsData = amide_wave_abundance, data_scale = "log2")

# QC-RFSC
amide_qcrfsc_abundance <- bc_qcrfsc(omicsData = pmart_amideFilt_abundance,qc_cname = "group",qc_val = "QC",
                          order_cname = "Injection_order",ntree = 500, keep_qc = FALSE,group_cname = "group")
amide_qcrfsc <- edata_transform(amide_qcrfsc_abundance, data_scale = "log2")

## ----out.width = "33%"--------------------------------------------------------

p1 <- plot(dim_reduction(omicsData = pmart_amideFilt),omicsData = pmart_amide,color_by = "batch") + labs(title = "Amide: Unadjusted")
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

## -----------------------------------------------------------------------------
data(pmart_mixFilt)

## ----warning = FALSE----------------------------------------------------------
# SCALING METHODS
# range scaling
mix_range <- bc_range(omicsData = pmart_mixFilt)

# pareto scaling
mix_pareto <- bc_pareto(omicsData = pmart_mixFilt)

# power scaling
mix_power <- bc_power(omicsData = pmart_mixFilt)

# INTERNAL STANDARDS/NEGATIVE CONTROLS
# RUV-random
mix_ruv <- bc_ruvRandom(omicsData = pmart_mixFilt, nc_cname = "tag",nc_val = "IS", k = 3)

# NOMIS
pmart_mixFilt_abundance <- edata_transform(omicsData = pmart_mixFilt, data_scale = "abundance")
mix_nomis_abundance <- bc_nomis(omicsData = pmart_mixFilt_abundance, is_cname = "tag", is_val = "IS")
mix_nomis <- edata_transform(omicsData = mix_nomis_abundance, data_scale = "log2")

# OTHER METHODS
mix_combat <- bc_combat(omicsData = pmart_mixFilt)

mix_waveica_abundance <- bc_waveica(omicsData = pmart_mixFilt_abundance,batch_cname = "BatchNum",
                                    version = "WaveICA",alpha = 0, K = 20, cutoff_batch = 0.05, 
                                    cutoff_group = 0.05,negative_to_na = TRUE)
mix_waveica <- edata_transform(omicsData = mix_waveica_abundance, data_scale = "log2")

## ----out.width = "33%"--------------------------------------------------------
p1 <- plot(dim_reduction(omicsData = pmart_mixFilt),omicsData = pmart_mixFilt,color_by = "BatchNum") + labs(title = "Mix: Unadjusted")
p2 <- plot(dim_reduction(omicsData = mix_ruv),omicsData = mix_ruv,color_by = "BatchNum") + labs(title = "Mix: ruv-Random")
p3 <- plot(dim_reduction(omicsData = mix_nomis),omicsData = mix_nomis,color_by = "BatchNum") + labs(title = "Mix: NOMIS")
p4 <- plot(dim_reduction(omicsData = mix_combat),omicsData = mix_combat,color_by = "BatchNum") + labs(title = "Mix: ComBat")
p5 <- plot(dim_reduction(omicsData = mix_range),omicsData = mix_range,color_by = "BatchNum") + labs(title = "Mix: Range")
p6 <- plot(dim_reduction(omicsData = mix_power),omicsData = mix_power,color_by = "BatchNum") + labs(title = "Mix: Power")
p7 <- plot(dim_reduction(omicsData = mix_pareto),omicsData = mix_pareto,color_by = "BatchNum") + labs(title = "Mix: Pareto")
p8 <- plot(dim_reduction(omicsData = mix_waveica),omicsData = mix_waveica,color_by = "BatchNum") + labs(title = "Mix: WaveICA")

p1;p2;p3;p4;p5;p6;p7;p8

