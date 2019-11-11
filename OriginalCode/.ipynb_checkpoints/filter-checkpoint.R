# Make the plots of mean and variance vs. theory, with error bars from bootstrapped confidence
# intervals

# adjust as appropriate
setwd("/media/everyan/NewVolume/files/workspace/Overdose")
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(ggmap)

block1 = read.xlsx("EMSdata/KEMSIS-KY-Lexington-1.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
block2 = read.xlsx("EMSdata/KEMSIS-KY-Lexington-2.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
block3 = read.xlsx("EMSdata/KEMSIS-KY-Lexington-3.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
block4 = read.xlsx("EMSdata/KEMSIS-KY-Lexington-4.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
block5 = read.xlsx("EMSdata/KEMSIS-KY-Lexington-5.xlsx", sheet = 1, startRow = 1, colNames = TRUE)
block = do.call("rbind", list(block1, block2, block3, block4, block5))

myvars <- c("Scene.Incident.Postal.Code.(eScene.19)", "Situation.Provider.Primary.Impression.(eSituation.11)", "Situation.Provider.Secondary.Impression.List.(eSituation.12)", "Situation.Symptom.Onset.Date.Time.(eSituation.01)", "Medication.Administered.Date.Time.(eMedications.01)", "Incident.Unit.Notified.By.Dispatch.Date.Time.(eTimes.03)")
data <- block[myvars]
names(data) <- c("ZIP", "L1", "L2", "T1","T2","T3")
data <- data %>% filter(str_detect(L1, "F11.|T40.1|T40.2")|str_detect(L2, "F11.|T40.1|T40.2"))
    
data[data == ""] <- NA
#data$T1 <- ifelse(is.na(data$T1),data$T2,data$T1)
data$T3 <-convertToDate(data$T3)
data$month <- as.POSIXct(as.character(data$T3), format="%Y-%m-%d")
data$month <- format(data$month, "%Y/%m")
data$ZIP <- substr(data$ZIP, 0,5)
table1 <- table(data$ZIP, data$month)
data %>% group_by(ZIP, month) %>% 
  summarise(Count = n(), Variance=var(n()))
write.csv(table1, file = "lexingtonOpioid.csv")
# block1 = read.csv("cincinnati-ems-overdose-Geo.csv")
# myvars <- c("Create.Time.Incident", "Latitude.X", "Longitude.X")
# data1 <- block1[myvars]
# data1$month <- as.POSIXct(as.character(data1$"Create.Time.Incident"), format="%m/%d/%Y %I:%M:%S %p")
# data1$month <- format(data1$month, "%Y/%m")
# coords <- as.matrix(data1[,2:3])
# data1$ZIP <- revgeocode(coords, output="more")

block1 = read.csv("shp/LZipcsv.csv", header=TRUE, sep=",")
block2 = read.csv("lexingtonOpioid.csv", header=TRUE, sep=",")
colnames(block2)[1] <- "ZIP"
block2$opiodTotal <- apply(block2[2:31], 1, sum)
block3<-merge(x = block1, y = block2, by = "ZIP", all.x = TRUE)
write.csv(block3, file = "lexingtonOpioidData.csv")

block4 = read.csv("EMSdata/cincinnati_ems_data_tractID.csv", header=TRUE, sep=",")
myvars <- c("CREATE_TIME_INCIDENT", "INCIDENT_TYPE_ID", "tract_fipscode")
data2 <- block4[myvars]
names(data2) <- c("T1","L1","ZIP")
data2 <- data2 %>% filter(str_detect(L1, "HERO|23C5"))
data2[data2 == ""] <- NA
#data2$T1 <-convertToDate(data2$T1)
data2$month <- as.POSIXct(as.character(data2$"T1"), format="%m/%d/%Y %I:%M:%S %p")
data2$month <- format(data2$month, "%Y/%m")
#data2$ZIP <- substr(data2$ZIP, 0,5)
table2 <- table(data2$ZIP, data2$month)
data2 %>% group_by(ZIP, month) %>% 
summarise(Count = n(), Variance=var(n()))
write.csv(table2, file = "cincinnatiOpioid.csv")

block1 = read.csv("shp/CincinnatiCityTracts.csv", header=TRUE, sep=",")
block2 = read.csv("cincinnatiOpioid.csv", header=TRUE, sep=",")
colnames(block2)[1] <- "Tract"
block2$opiodTotal <- apply(block2[2:31], 1, sum)
block3<-merge(x = block1, y = block2, by = "Tract", all.x = TRUE)
write.csv(block3, file = "cincinnatiOpioidData.csv")

block1 = read.csv("shp/CZipcsv.csv", header=TRUE, sep=",")
block2 = read.csv("lexingtonOpioid.csv", header=TRUE, sep=",")
block3 = read.csv("cincinnatiOpioid.csv", header=TRUE, sep=",")
block3 <- block3[ -c(2:13) ]  # delete columns 2 through 13
block3 <- block3[ -c(26:28) ]  # delete columns 2 through 13
block2 <- block2[ -c(26:31) ]  # delete columns 2 through 13
colnames(block2)[1] <- "ZIP"
colnames(block3)[1] <- "ZIP"
block2$opiodTotal <- apply(block2[2:25], 1, sum)
block3$opiodTotal <- apply(block3[2:25], 1, sum)
block4 = do.call("rbind", list(block2, block3))
block5<-merge(x = block1, y = block4, by = "ZIP", all.x = TRUE)
write.csv(block5, file = "cincinnatiOpioidData.csv")
