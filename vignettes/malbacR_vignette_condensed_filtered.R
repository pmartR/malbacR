## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE)

## -----------------------------------------------------------------------------
library(malbacR)
library(pmartR)
library(ggplot2)

## -----------------------------------------------------------------------------
data(pmart_amideFilt)

## ---- results = 'hide'--------------------------------------------------------
# SCALING METHODS
# range scaling
amide_range <- bc_range(omicsData = pmart_amideFilt)

# power scaling
amide_power <- bc_power(omicsData = pmart_amideFilt)

# pareto scaling
amide_pareto <- bc_pareto(omicsData = pmart_amideFilt)

# QUALITY CONTROL METHODS
# TIGER
amide_tiger <- bc_tiger(omicsData = pmart_amideFilt,sampletype_cname = "group",
                        test_val = "QC",injection_cname = "Injection_order")

# QC-RLSC
amide_qcrlsc <- bc_qcrlsc(omicsData = pmart_amideFilt,block_cname = "batch",
                          qc_cname = "group", qc_val = "QC", order_cname = "Injection_order",
                          missing_thresh = 0.5, rsd_thresh = 0.3, backtransform  = FALSE)

# SERRF
amide_serrf <- bc_serrf(omicsData = pmart_amideFilt,sampletype_cname = "group",test_val = "QC")

# OTHER METHODS
# ComBat
amide_combat <- bc_combat(omicsData = pmart_amideFilt, use_groups = FALSE)

# EigenMS
amide_eigen <- bc_eigenMS(omicsData = pmart_amideFilt)

# WaveICA2.0
amide_wave <- bc_waveica(omicsData = pmart_amideFilt, injection_cname = "Injection_order",
                         alpha = 0, cutoff = 0.1, K = 10)

## ---- out.width = "33%"-------------------------------------------------------
pmart_amide <- group_designation(pmart_amideFilt,main_effects = "batch")
amide_range <- group_designation(amide_range,main_effects = "batch")
amide_power <- group_designation(amide_power,main_effects = "batch")
amide_pareto <- group_designation(amide_pareto,main_effects = "batch")
amide_tiger <- group_designation(amide_tiger,main_effects = "batch")
amide_qcrlsc <- group_designation(amide_qcrlsc,main_effects = "batch")
amide_combat <- group_designation(amide_combat,main_effects = "batch")
amide_wave <- group_designation(amide_wave,main_effects = "batch")
amide_serrf <- group_designation(amide_serrf,main_effects = "batch")

p1 <- plot(dim_reduction(omicsData = pmart_amide))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: Unadjusted")
p2 <- plot(dim_reduction(omicsData = amide_range))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: Range")
p3 <- plot(dim_reduction(omicsData = amide_power))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: Power")
p4 <- plot(dim_reduction(omicsData = amide_pareto))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: Pareto")
p5 <- plot(dim_reduction(omicsData = amide_tiger))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: TIGER")
p6 <- plot(dim_reduction(omicsData = amide_qcrlsc))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: QC-RLSC")
p7 <- plot(dim_reduction(omicsData = amide_combat))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: ComBat")
p8 <- plot(dim_reduction(omicsData = amide_wave))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: WaveICA2.0")
p9 <- plot(dim_reduction(omicsData = amide_serrf))+ scale_colour_discrete(name="Batch") + labs(title = "Amide: SERRF")

p1;p2;p3;p4;p5;p6;p7;p8;p9

## -----------------------------------------------------------------------------
data(pmart_mixFilt)

## ---- warning = FALSE---------------------------------------------------------
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
mix_nomis <- bc_nomis(omicsData = pmart_mixFilt, is_cname = "tag", is_val = "IS", num_pc = 2)

# OTHER METHODS
mix_combat <- bc_combat(omicsData = pmart_mixFilt)

## ----out.width = "33%"--------------------------------------------------------
pmart_mix <- group_designation(pmart_mixFilt,main_effects = "BatchNum")
mix_range <- group_designation(mix_range,main_effects = "BatchNum")
mix_power <- group_designation(mix_power,main_effects = "BatchNum")
mix_pareto <- group_designation(mix_pareto,main_effects = "BatchNum")
mix_combat <- group_designation(mix_combat,main_effects = "BatchNum")
mix_ruv <- group_designation(mix_ruv,main_effects = "BatchNum")
mix_nomis <- group_designation(mix_nomis,main_effects = "BatchNum")

p1 <- plot(dim_reduction(omicsData = pmart_mix)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: Unadjusted")
p2 <- plot(dim_reduction(omicsData = mix_ruv)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: RUV-random")
p3 <- plot(dim_reduction(omicsData = mix_nomis)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: NOMIS")
p4 <- plot(dim_reduction(omicsData = mix_combat)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: ComBat")
p5 <- plot(dim_reduction(omicsData = mix_range)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: Range")
p6 <- plot(dim_reduction(omicsData = mix_power)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: Power")
p7 <- plot(dim_reduction(omicsData = mix_pareto)) + scale_colour_discrete(name="Batch") + labs(title = "Mix: Pareto")

p1;p2;p3;p4;p5;p6;p7

