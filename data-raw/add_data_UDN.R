### UDN Data Formatting ###

## Background batches were all run consecutively over the course of a week (June 2016)
## Pilot1 was standalone but only a couple days after the background (June 2016)
## Pilot2 was standalone a couple months later (Oct 2016)
## Project01 was standalone (Dec 2016)
## Project02 was standalone (Feb 2017)
## Project04 & Project06 were consecutive (Sep 2017)
## Project03 & Project05 were consecutive after a couple days from 4 & 6 (Sep 2017)
## Project07 was standalone (Oct 2017)
## Project08 was standalone (Dec 2017)
## Project09 was standalone (Feb 2018)
## Project10 & Project11 were consecutive (May 2018)
## Project12 was standalone (Oct 2018)

# So, overall run order was:
# BG, Pilot1, Pilot2, Project01, Project02, Project04 & Project06, Project03 &
# Project05, Project07, Project08, Project09, Project10 & Project11, Project12

rm(list = ls())

library(stringr)
library(pmartR)

# read in UDN data files (from Jen) #
mydir <- "/Volumes/mq_DANCE/Example_Data/UDN_Data/Plasma_Metabolomics_Raw/Raw data/"
all_files <- list.files(mydir)
rds_files <- all_files[grep(".rds", all_files)]

rds_names <- unlist(lapply(rds_files, function(x) str_replace(x, ".rds", "")))
rds_names <- unlist(lapply(rds_names, function(x) str_replace(x, "_v2", "")))
rds_names <- unlist(lapply(rds_names, function(x) str_replace(x, "_v3b", "")))
rds_names <- unlist(lapply(rds_names, function(x) str_replace(x, "_v4", "")))
rds_names <- unlist(lapply(rds_names, function(x) str_replace(x, "_v3", "")))

# read in the files (metabData objects, but from older version of pmartR)
for(i in 1:length(rds_files)){
  assign(rds_names[i], readRDS(paste0(mydir, rds_files[i])))
}

# attributes(bg_metab)$data_info
## log2
## not normalized


### re-create the metabData objects ###
update_metabData <- function(cur_data, fdata_cname = "Names"){
  newdata <- as.metabData(e_data = cur_data$e_data,
                          f_data = cur_data$f_data,
                          e_meta = cur_data$e_meta,
                          edata_cname = get_edata_cname(cur_data),
                          fdata_cname = fdata_cname,
                          emeta_cname = get_emeta_cname(cur_data),
                          data_scale = "log2",
                          check.names = FALSE)

  qc_inds <- grep("QC", newdata$f_data[, fdata_cname])
  not_qc_inds <- setdiff(1:nrow(newdata$f_data), qc_inds)

  if(any(is.na(newdata$f_data$Age))){
    na_inds <- which(is.na(newdata$f_data$Age))
    not_na_inds <- setdiff(1:nrow(newdata$f_data), na_inds)
    # if it's a QC sample that's missing Age, give it the value "QC"
    if(length(intersect(na_inds, qc_inds)) > 0){
      cur_inds <- intersect(na_inds, qc_inds)
      newdata$f_data$Age[cur_inds] <- "QC"
    }
    # if it's a non-QC sample that's missing Age, give it the value "Unknown"
    if(length(intersect(na_inds, not_qc_inds)) > 0){
      cur_inds <- intersect(na_inds, not_qc_inds)
      newdata$f_data$Age[cur_inds] <- "Unknown"
    }
  }

  if(any(is.na(newdata$f_data$Sex))){
    na_inds <- which(is.na(newdata$f_data$Sex))
    not_na_inds <- setdiff(1:nrow(newdata$f_data), na_inds)
    # if it's a QC sample that's missing Sex, give it the value "QC"
    if(length(intersect(na_inds, qc_inds)) > 0){
      cur_inds <- intersect(na_inds, qc_inds)
      newdata$f_data$Sex[cur_inds] <- "QC"
    }
    # if it's a non-QC sample that's missing Sex, give it the value "Unknown"
    if(length(intersect(na_inds, not_qc_inds)) > 0){
      cur_inds <- intersect(na_inds, not_qc_inds)
      newdata$f_data$Sex[cur_inds] <- "Unknown"
    }
  }

  newdata <- group_designation(newdata, main_effects = c("Sex", "Age"))

  if(length(newdata$f_data$RunOrder) == 0 & length(newdata$f_data$OverallRunOrder) > 0){
    newdata$f_data$RunOrder <- newdata$f_data$OverallRunOrder
  }else{
    if(length(newdata$f_data$RunOrder) == 0 & length(newdata$f_data$Run_Order_M) > 0){
      newdata$f_data$RunOrder <- newdata$f_data$Run_Order_M
    }
  }

  return(newdata)
}

bg <- update_metabData(cur_data = bg_metab)
summary(bg)
plot(bg, order_by = "RunOrder", color_by = "Batch")
rm(bg_metab)

pilot01 <- update_metabData(cur_data = pilot01_metab)
summary(pilot01)
plot(pilot01, order_by = "RunOrder", color_by = "Batch")
rm(pilot01_metab)

pilot02 <- update_metabData(cur_data = pilot02_metab)
summary(pilot02)
plot(pilot02, order_by = "RunOrder", color_by = "Batch")
rm(pilot02_metab)

project01 <- update_metabData(cur_data = project01_metab)
summary(project01)
plot(project01, order_by = "RunOrder", color_by = "Batch")
rm(project01_metab)

project02 <- update_metabData(cur_data = project02_metab)
summary(project02)
plot(project02, order_by = "RunOrder", color_by = "Batch")
rm(project02_metab)

project03 <- update_metabData(cur_data = project03_metab)
summary(project03)
plot(project03, order_by = "RunOrder", color_by = "Batch")
rm(project03_metab)

project04 <- update_metabData(cur_data = project04_metab)
summary(project04)
plot(project04, order_by = "RunOrder", color_by = "Batch")
rm(project04_metab)

project05 <- update_metabData(cur_data = project05_metab)
summary(project05)
plot(project05, order_by = "RunOrder", color_by = "Batch")
rm(project05_metab)

project06 <- update_metabData(cur_data = project06_metab)
summary(project06)
plot(project06, order_by = "RunOrder", color_by = "Batch")
rm(project06_metab)

project07 <- update_metabData(cur_data = p7_metab)
summary(project07)
plot(project07, order_by = "RunOrder", color_by = "Batch")
rm(p7_metab)

project08 <- update_metabData(cur_data = p8_metab)
summary(project08)
plot(project08, order_by = "RunOrder", color_by = "Batch")
rm(p8_metab)

project09 <- update_metabData(cur_data = p9_metab)
summary(project09)
plot(project09, order_by = "RunOrder", color_by = "Batch")
rm(p9_metab)

project10 <- update_metabData(cur_data = p10_metab)
summary(project10)
plot(project10, order_by = "RunOrder", color_by = "Batch")
rm(p10_metab)

project11 <- update_metabData(cur_data = p11_metab)
summary(project11)
plot(project11, order_by = "RunOrder", color_by = "Batch")
rm(p11_metab)

project12 <- update_metabData(cur_data = p12_metab)
summary(project12)
plot(project12, order_by = "RunOrder", color_by = "Batch")
rm(p12_metab)


### Now let's get Overall Run Order Updated in each of these objects ###
# Overall run order was:
# BG, Pilot1, Pilot2, Project01, Project02, Project04 & Project06, Project03 &
# Project05, Project07, Project08, Project09, Project10 & Project11, Project12

run_orders <- data.frame(SampleID = c(bg$f_data$Names,
                                      pilot01$f_data$Names,
                                      pilot02$f_data$Names,
                                      project01$f_data$Names,
                                      project02$f_data$Names,
                                      project04$f_data$Names,
                                      project06$f_data$Names,
                                      project03$f_data$Names,
                                      project05$f_data$Names,
                                      project07$f_data$Names,
                                      project08$f_data$Names,
                                      project09$f_data$Names,
                                      project10$f_data$Names,
                                      project11$f_data$Names,
                                      project12$f_data$Names),
                         BatchNum = c(bg$f_data$Batch, # 10 BG batches
                                      rep(11, nrow(pilot01$f_data)), # pilot 1
                                      rep(12, nrow(pilot02$f_data)), # pilot 2
                                      rep(13, nrow(project01$f_data)), # project 1
                                      rep(14, nrow(project02$f_data)), # project 2
                                      rep(15, nrow(project04$f_data)), # project 4
                                      rep(16, nrow(project06$f_data)), # project 6
                                      rep(17, nrow(project03$f_data)), # project 3
                                      rep(18, nrow(project05$f_data)), # project 5
                                      rep(19, nrow(project07$f_data)), # project 7
                                      rep(20, nrow(project08$f_data)), # project 8
                                      rep(21, nrow(project09$f_data)), # project 9
                                      rep(22, nrow(project10$f_data)), # project 10
                                      rep(23, nrow(project11$f_data)), # project 11
                                      rep(24, nrow(project12$f_data))), # project 12
                         BatchName = c(rep("BG01", length(which(bg$f_data$Batch == 1))),
                                       rep("BG02", length(which(bg$f_data$Batch == 2))),
                                       rep("BG03", length(which(bg$f_data$Batch == 3))),
                                       rep("BG04", length(which(bg$f_data$Batch == 4))),
                                       rep("BG05", length(which(bg$f_data$Batch == 5))),
                                       rep("BG06", length(which(bg$f_data$Batch == 6))),
                                       rep("BG07", length(which(bg$f_data$Batch == 7))),
                                       rep("BG08", length(which(bg$f_data$Batch == 8))),
                                       rep("BG09", length(which(bg$f_data$Batch == 9))),
                                       rep("BG010", length(which(bg$f_data$Batch == 10))),
                                       rep("Pi01", nrow(pilot01$f_data)),
                                       rep("Pi02", nrow(pilot02$f_data)),
                                       rep("Pr01", nrow(project01$f_data)),
                                       rep("Pr02", nrow(project02$f_data)),
                                       rep("Pr04", nrow(project04$f_data)),
                                       rep("Pr06", nrow(project06$f_data)),
                                       rep("Pr03", nrow(project03$f_data)),
                                       rep("Pr05", nrow(project05$f_data)),
                                       rep("Pr07", nrow(project07$f_data)),
                                       rep("Pr08", nrow(project08$f_data)),
                                       rep("Pr09", nrow(project09$f_data)),
                                       rep("Pr10", nrow(project10$f_data)),
                                       rep("Pr11", nrow(project11$f_data)),
                                       rep("Pr12", nrow(project12$f_data))),
                         RunOrderOverall = NA,
                         RunOrderBatch = c(bg$f_data$RunOrder,
                                           pilot01$f_data$RunOrder,
                                           pilot02$f_data$RunOrder,
                                           project01$f_data$RunOrder,
                                           project02$f_data$RunOrder,
                                           project04$f_data$RunOrder,
                                           project06$f_data$RunOrder,
                                           project03$f_data$RunOrder,
                                           project05$f_data$RunOrder,
                                           project07$f_data$RunOrder,
                                           project08$f_data$RunOrder,
                                           project09$f_data$RunOrder,
                                           project10$f_data$RunOrder,
                                           project11$f_data$RunOrder,
                                           project12$f_data$RunOrder),
                         Sex = c(bg$f_data$Sex,
                                 pilot01$f_data$Sex,
                                 pilot02$f_data$Sex,
                                 project01$f_data$Sex,
                                 project02$f_data$Sex,
                                 project04$f_data$Sex,
                                 project06$f_data$Sex,
                                 project03$f_data$Sex,
                                 project05$f_data$Sex,
                                 project07$f_data$Sex,
                                 project08$f_data$Sex,
                                 project09$f_data$Sex,
                                 project10$f_data$Sex,
                                 project11$f_data$Sex,
                                 project12$f_data$Sex),
                         Age = c(bg$f_data$Age,
                                 pilot01$f_data$Age,
                                 pilot02$f_data$Age,
                                 project01$f_data$Age,
                                 project02$f_data$Age,
                                 project04$f_data$Age,
                                 project06$f_data$Age,
                                 project03$f_data$Age,
                                 project05$f_data$Age,
                                 project07$f_data$Age,
                                 project08$f_data$Age,
                                 project09$f_data$Age,
                                 project10$f_data$Age,
                                 project11$f_data$Age,
                                 project12$f_data$Age))
run_orders$RunOrderOverall <- 1:nrow(run_orders)
# now run_orders is my massive f_data

### create massive metabData object :)
f_data <- run_orders

e_data <- merge(x = bg$e_data, y = pilot01$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = pilot02$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project01$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project02$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project04$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project06$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project03$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project05$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project07$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project08$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project09$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project10$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project11$e_data, by = "Metabolite", all = TRUE)
e_data <- merge(x = e_data, y = project12$e_data, by = "Metabolite", all = TRUE)

all(names(e_data)[-1] %in% f_data$SampleID)

# e_meta <- merge(x = bg$e_meta, y = pilot01$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = pilot02$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project01$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project02$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project04$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project06$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project03$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project05$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project07$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project08$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project09$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project10$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project11$e_meta, by = "Metabolite", all = TRUE)
# e_meta <- merge(x = e_meta, y = project12$e_meta, by = "Metabolite", all = TRUE)

alldata <- as.metabData(e_data = e_data,
                        f_data = f_data,
                        edata_cname = "Metabolite",
                        fdata_cname = "SampleID",
                        data_scale = "log2",
                        check.names = FALSE)
alldata <- group_designation(alldata, main_effects = c("BatchName"))

# library(ggplot2)
# library(plotly)
# png("/Volumes/mq_DANCE/Example_Data/UDN_Data/boxplots_alldata.png", width = 2000)
# plot(alldata, color_by = "BatchName", order_by = "RunOrderOverall", interactive = TRUE)
#   # theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
# dev.off()



# subset to just background data to understand why the plot isn't ordering the way I expect
myfilt <- custom_filter(alldata, f_data_keep = alldata$f_data$SampleID[grep("BG", alldata$f_data$SampleID)])
bgdata <- applyFilt(myfilt, alldata)
saveRDS(bgdata, "/Volumes/mq_DANCE/Example_Data/UDN_Data/bgdata.RDS")
bgdata <- readRDS("/Volumes/mq_DANCE/Example_Data/UDN_Data/bgdata.RDS")
plot(bgdata, color_by = "BatchName", order_by = "RunOrderOverall", interactive = TRUE)
plot(bgdata$f_data$BatchNum, bgdata$f_data$RunOrderOverall)


plot_omicsData_local(bgdata, order_by = "RunOrderOverall", color_by = "BatchName",
                     facet_by = NULL, facet_cols = NULL,
                     interactive = TRUE, x_lab = NULL, y_lab = NULL, x_lab_size = NULL, y_lab_size = NULL,
                     x_lab_angle = NULL, title_lab = NULL, title_lab_size = NULL, legend_lab = NULL,
                     legend_position = NULL, ylimit = NULL, bw_theme = TRUE, palette = NULL,
                     use_VizSampNames = FALSE)


# alldata <- readRDS("/Volumes/mq_DANCE/Example_Data/UDN_Data/alldata.RDS")

# subset to just QC NIST samples
tokeep <- grep("NIST", alldata$f_data$SampleID)
myfilt <- custom_filter(alldata, f_data_keep = alldata$f_data$SampleID[tokeep])
allnist <- applyFilt(myfilt, alldata)
summary(allnist)
View(allnist$f_data)
saveRDS(allnist, "/Volumes/mq_DANCE/Example_Data/UDN_Data/allnist.RDS")

qcdata <- readRDS("/Volumes/mq_DANCE/Example_Data/UDN_Data/allnist.RDS")


#### add the data to the package
udn_metab_all <- alldata
usethis::use_data(udn_metab_all, overwrite = TRUE)

udn_metab_qc <- qcdata
usethis::use_data(udn_metab_qc, overwrite = TRUE)
