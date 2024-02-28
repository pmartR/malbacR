# damon leach
# march 14, 2023
# malbac figures
  
# load in necessary libraries
library(malbacR)
library(pmartR)
library(patchwork)
library(ggplot2)

################################### figure 1 ###################################
# load in data
data(pmart_amide)

# log2 transformation
pmart_amide <- edata_transform(pmart_amide,"log2")
# group designation
pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
# normalize data
pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
                                apply_norm = TRUE,backtransform = TRUE)

# batch correction methods
# serrf (including imputation)
impObj <- imputation(omicsData = pmart_amide)
amide_imp <- apply_imputation(imputeData = impObj, omicsData = pmart_amide)
amide_serrf <- bc_serrf(omicsData = amide_imp,sampletype_cname = "group",test_val = "QC")
# combat
amide_combat <- bc_combat(omicsData = pmart_amide, use_groups = FALSE)
# waveica (using imputed dataset created during serrf)
amide_wave <- bc_waveica(omicsData = amide_imp, injection_cname = "Injection_order",
                         alpha = 0, cutoff = 0.1, K = 10)

# data visualization
# change group designation to be batch for plotting purposes
pmart_amide <- group_designation(pmart_amide,main_effects = "batch")
amide_combat <- group_designation(amide_combat,main_effects = "batch")
amide_wave <- group_designation(amide_wave,main_effects = "batch")
amide_serrf <- group_designation(amide_serrf,main_effects = "batch")
# unadjusted
pmart_amide_res <- dim_reduction(omicsData = pmart_amide)
amide1 <- data.frame(sample = pmart_amide_res$SampleID, Batch = attr(pmart_amide_res,"group_DF")$Group,
                     pc1 = pmart_amide_res$PC1, pc2 = pmart_amide_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(pmart_amide_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(pmart_amide_res,"R2")[2],3),")"),
       title = "Amide: Unadjusted")
# combat
amide_combat_res <- dim_reduction(omicsData = amide_combat)
amide2 <- data.frame(sample = amide_combat_res$SampleID, Batch = attr(amide_combat_res,"group_DF")$Group,
                     pc1 = amide_combat_res$PC1, pc2 = amide_combat_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(amide_combat_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(amide_combat_res,"R2")[2],3),")"),
       title = "Amide: ComBat")
# waveica
amide_wave_res <- dim_reduction(omicsData = amide_wave)
amide3 <- data.frame(sample = amide_wave_res$SampleID, Batch = attr(amide_wave_res,"group_DF")$Group,
                     pc1 = amide_wave_res$PC1, pc2 = amide_wave_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.50,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(amide_wave_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(amide_wave_res,"R2")[2],3),")"),
       title = "Amide: WaveICA2.0")
# serrf
amide_serrf_res <- dim_reduction(omicsData = amide_serrf)
amide4 <- data.frame(sample = amide_serrf_res$SampleID, Batch = attr(amide_serrf_res,"group_DF")$Group,
                     pc1 = amide_serrf_res$PC1, pc2 = amide_serrf_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(amide_serrf_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(amide_serrf_res,"R2")[2],3),")"),
       title = "Amide: SERRF")
# patchwork to show plot the four figures
(amide1|amide2)/(amide3|amide4)

################################### figure 2 ###################################
# load in data
data(pmart_mix)

# log2 transformation
pmart_mix <- edata_transform(pmart_mix,"log2")
# group designation
pmart_mix <- group_designation(pmart_mix,main_effects = "BatchNum",batch_id = "BatchNum")
# normalize data
pmart_mix <- normalize_global(pmart_mix,subset_fn = "all",norm_fn = "median",
                              apply_norm = TRUE,backtransform = TRUE)

# batch correction methods
# ruv-random
mix_ruv <- bc_ruvRandom(omicsData = pmart_mix, nc_cname = "tag",nc_val = "IS", k = 3)
# nomis
mix_nomis <- bc_nomis(omicsData = pmart_mix, is_cname = "tag", is_val = "IS")
# combat
mix_combat <- bc_combat(omicsData = pmart_mix)

# data visualization
# change group designation to be batch for plotting purposes
pmart_mix <- group_designation(pmart_mix,main_effects = "BatchNum")
mix_combat <- group_designation(mix_combat,main_effects = "BatchNum")
mix_ruv <- group_designation(mix_ruv,main_effects = "BatchNum")
mix_nomis <- group_designation(mix_nomis,main_effects = "BatchNum")

# unadjusted
pmart_mix_res <- dim_reduction(omicsData = pmart_mix)
mix1 <- data.frame(sample = pmart_mix_res$SampleID, Batch = attr(pmart_mix_res,"group_DF")$Group,
                   pc1 = pmart_mix_res$PC1, pc2 = pmart_mix_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(pmart_mix_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(pmart_mix_res,"R2")[2],3),")"),
       title = "Mix: Unadjusted")
# ruv random
mix_ruv_res <- dim_reduction(omicsData = mix_ruv)
mix2 <- data.frame(sample = mix_ruv_res$SampleID, Batch = attr(mix_ruv_res,"group_DF")$Group,
                   pc1 = mix_ruv_res$PC1, pc2 = mix_ruv_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(mix_ruv_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(mix_ruv_res,"R2")[2],3),")"),
       title = "Mix: RUV-random")
# nomis
mix_nomis_res <- dim_reduction(omicsData = mix_nomis)
mix3 <- data.frame(sample = mix_nomis_res$SampleID, Batch = attr(mix_nomis_res,"group_DF")$Group,
                   pc1 = mix_nomis_res$PC1, pc2 = mix_nomis_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(mix_nomis_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(mix_nomis_res,"R2")[2],3),")"),
       title = "Mix: NOMIS")
# combat
mix_combat_res <- dim_reduction(omicsData = mix_combat)
mix4 <- data.frame(sample = mix_combat_res$SampleID, Batch = attr(mix_combat_res,"group_DF")$Group,
                   pc1 = mix_combat_res$PC1, pc2 = mix_combat_res$PC2) %>%
  ggplot(aes(x = pc1, y = pc2,color = Batch)) + 
  geom_point(alpha = 0.5,size = 4) + 
  theme_bw() + 
  labs(x = paste0("PC1 (R^2 = ",round(attr(mix_combat_res,"R2")[1],3),")"),
       y = paste0("PC1 (R^2 = ",round(attr(mix_combat_res,"R2")[2],3),")"),
       title = "Mix: ComBat")
# patchwork to show plot the four figures
(mix1|mix2)/(mix3|mix4)
