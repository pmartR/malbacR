# damon leach
# march 14, 2023
# data manipulation code

################################### dataset 1 ##################################
# retain seed after  running code
old_seed <- .Random.seed
on.exit(.Random.seed <- old_seed)

# load in library with original data
library(WaveICA2.0)
# load in dataset
data(Amide_data)

# create fdata
amide_fdata <- Amide_data %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  dplyr::select(SampleID,Injection_order,group,batch)

# create edata
amide_edata <- Amide_data %>%
  dplyr::select(-c(Injection_order,group,batch)) %>%
  t() %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "Molecule")

# there is no emeta in this object

# create metabObject
pmart_amide <- pmartR::as.metabData(e_data = amide_edata,f_data = amide_fdata,
                                    edata_cname = "Molecule",fdata_cname = "SampleID")
# sample only 500 of the molecules for easier user speed
set.seed(1)
good_mol <- sample(pmart_amide$e_data$Molecule,500)
# run custom filter with that sample of 500
cfilt2 <- custom_filter(pmart_amide,e_data_keep = good_mol)
# apply the custom filter
pmart_amide <- applyFilt(cfilt2,pmart_amide)

################################### dataset 2 ##################################
# load in library with original data
library(crmn)
# load in dataset
data(mix)

# run through their code to get the expression data
Y <- exprs(mix)
G <- with(pData(mix), model.matrix(~-1 + type))
isIS <- with(fData(mix),tag == "IS")

# create emeta
emeta <- fData(mix) %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "moleculeNum")
# find abundance values
abundance <- Y %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "moleculeNum")
# create edata
edata <- emeta %>%
  dplyr::left_join(abundance) %>%
  dplyr::select(-c(moleculeNum,mark,tag,synonym,RI,query,known)) %>%
  dplyr::rename("Metabolite" = "preferred")
# update emeta names
emeta <- emeta %>%
  dplyr::select(-moleculeNum) %>%
  dplyr::rename("Metabolite" = "preferred")
# create fdata
fdata <- data.frame(SampleID = colnames(edata)[-1],
                    BatchNum = rep(NA,42))
# update batch information in code
fdata <- fdata %>%
  dplyr::mutate(BatchNum = ifelse(stringr::str_detect(SampleID,"STDs_1")|stringr::str_detect(SampleID,"STDs1"),1,BatchNum),
                BatchNum = ifelse(stringr::str_detect(SampleID,"STDs_2")|stringr::str_detect(SampleID,"STDs2"),2,BatchNum),
                BatchNum = ifelse(stringr::str_detect(SampleID,"STDs_3")|stringr::str_detect(SampleID,"STDs3"),3,BatchNum))
colnames(edata) <- c("Metabolite",paste0("sample",seq(from = 1, to = 42, by = 1)))
# update sampleID to match edata
fdata$SampleID <- paste0("sample",seq(from = 1, to = 42, by = 1))

# create metabObject
pmart_mix <- as.metabData(e_data = edata,f_data = fdata,e_meta = emeta,
                          edata_cname = "Metabolite",fdata_cname = "SampleID",emeta_cname = "Metabolite")