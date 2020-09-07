rm(list = ls(all = TRUE)) 
setwd("Z:/Experiments_of_Maria/20200505_strain_extra")
ukbbdata<-readRDS("Z:/UKBB_40616/Phenotypes/41354/post_processed/my_ukb_data_41354.rda")

folderImagingDataPath <- c("Z:/UKBB_New_Data/data")

foldersNames <- list.dirs(folderImagingDataPath,
                          full.names = FALSE,
                          recursive = F)

ukbb_bridge <-data.frame(read.csv("Z:/UKBB_40616/Bridging file/Bridge40616_18545.csv"))

p<-match(foldersNames,ukbb_bridge$eid_18545)
names<-ukbb_bridge[p,]
pdata<-match(names$eid_40616,ukbbdata$eid)
age<-ukbbdata$age_when_attended_assessment_centre_f21003_2_0[pdata]
sex<-ukbbdata$sex_f31_0_0[pdata]
library(plyr)

sex<-as.matrix(mapvalues(sex, from = c("Female", "Male"), to = c("0","1")))
sex<-as.numeric(as.character(sex))

bsa<-ukbbdata$body_surface_area_f22427_2_0[pdata]


mSBP<-cbind(rowMeans(cbind(ukbbdata$systolic_blood_pressure_manual_reading_f93_2_0[pdata],ukbbdata$systolic_blood_pressure_manual_reading_f93_2_1[pdata]), na.rm = TRUE, dims = 1),
            rowMeans(cbind(ukbbdata$systolic_blood_pressure_automated_reading_f4080_2_0[pdata], ukbbdata$systolic_blood_pressure_automated_reading_f4080_2_1[pdata]), na.rm = TRUE, dims = 1))
sbp<-rowMeans(mSBP, na.rm = TRUE, dims = 1)
mDBP<-cbind(rowMeans(cbind(ukbbdata$diastolic_blood_pressure_manual_reading_f94_2_0[pdata], ukbbdata$diastolic_blood_pressure_manual_reading_f94_2_1[pdata]), na.rm = TRUE, dims = 1),
            rowMeans(cbind(ukbbdata$diastolic_blood_pressure_automated_reading_f4079_2_0[pdata], ukbbdata$diastolic_blood_pressure_automated_reading_f4079_2_1[pdata]), na.rm = TRUE, dims = 1))
dbp<-rowMeans(mDBP, na.rm = TRUE, dims = 1)

ukb_datatble<-as.data.frame(cbind(names$eid_18545,age, sex, bsa,sbp,dbp))
colnames(ukb_datatble)<-c("Encoded_ID","Age","Sex","BSA","SBP","DBP")
write.csv(ukb_datatble,"phenotype_ukbb_1800.csv", col.names = TRUE, row.names = FALSE)
# n<-which(is.na(ukb_datatble$Age))
