## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(message = FALSE)

## -----------------------------------------------------------------------------
library(malbacR)
library(pmartR)
library(ggplot2)

## -----------------------------------------------------------------------------
data(pmart_amide)

## -----------------------------------------------------------------------------
pmart_amide <- edata_transform(pmart_amide,"log2")
pmart_amide <- group_designation(pmart_amide,main_effects = "group",batch_id = "batch")
pmart_amide <- normalize_global(pmart_amide,subset_fn = "all",norm_fn = "median",
                                apply_norm = TRUE,backtransform = TRUE)

## -----------------------------------------------------------------------------
# range scaling
amide_range <- bc_range(omicsData = pmart_amide)

# power scaling
amide_power <- bc_power(omicsData = pmart_amide)

# pareto scaling
amide_pareto <- bc_pareto(omicsData = pmart_amide)

## -----------------------------------------------------------------------------
tigerFilt <- tiger_filter(pmart_amide,sampletype_cname = "group",test_val = "QC")
pmart_amideFilt <- apply_tigerFilt(tigerFilt,pmart_amide)
amide_tiger <- bc_tiger(omicsData = pmart_amideFilt,sampletype_cname = "group",
                        test_val = "QC",injection_cname = "Injection_order")

