############################################

# Title: "Script for processing Luminex files"
# Authors: Bio-plex: "L. Bos", updated by: "R. van Amstel"
#          Luminex200: "R. van Amstel", updated by: T. Neumann
# Date: "05/2026"

####################
#library
library(readxl)
library(reshape2)

######## paths
#input paths
input_path_1 <- snakemake@input[[1]]

input_path_3 <- snakemake@input[[2]]
input_path_4 <- snakemake@input[[3]]

#output paths
outpath_path_1 <- snakemake@output[[1]]
outpath_path_2 <- snakemake@output[[2]]
Processed_files <- c(outpath_path_1, outpath_path_2)

outpath_path_4 <- snakemake@output[[3]]
outpath_path_5 <- snakemake@output[[4]]
CVs_files <- c(outpath_path_4, outpath_path_5)

outpath_path_7 <- snakemake@output[[5]]
outpath_path_8 <- snakemake@output[[6]]
Standard_curve_files <- c(outpath_path_7, outpath_path_8)

outpath_path_10 <- snakemake@output[[7]]
outpath_path_11 <- snakemake@output[[8]]
Standard_expected_files <- c(outpath_path_10, outpath_path_11)


### read plate layout in
df_plate_layout <- read.csv2(input_path_1, sep = ";")

all_files <- c(input_path_3, input_path_4)


for (i in 1:length(all_files)){
  
  tmp <- read_excel(all_files[i], sheet = "Obs Conc")
  tmp <- tmp[-c(3, 200),]
  
  tmp_beads  <- read_excel(all_files[i], sheet = "Bead Count")
  tmp_beads <- tmp_beads[-c(3, 200),]
  
  tmp_dilution  <- read_excel(all_files[i], sheet = "Dilution")
  tmp_dilution <- tmp_dilution[-c(3, 200),]
  
  tmp_expected <- read_excel(all_files[i], sheet = "Exp Conc")
  tmp_expected <- tmp_expected[-c(3, 200),]
  
  tmp_cv <- read_excel(all_files[i], sheet = "%CV")
  tmp_cv <- tmp_cv[-c(3, 200),]
  
  
  dochterplaat <- i#NB: Check the length of file title!
  
  ## Extract CVs 
  if(length(which(tmp_cv[,1]=="B"))>0){data_cvs<-data.frame(tmp_cv[c(
    which(tmp_cv[,1]=="B")[1]: which(tmp_cv[,2]=="C1")[1]),])
  }
  if(length(which(tmp_cv[,1]=="B"))==0){data_cvs<-data.frame(tmp_cv[c(
    which(tmp_cv[,1]=="eS1")[1]: which(tmp_cv[,1]=="C1")[1]),]) 
  } 
  colnames(data_cvs)<-as.character(tmp_cv[7,]) 
  colnames(data_cvs)[1]<-"sample_type"
  for (j in 4:ncol(data_cvs)){
    data_cvs[,j]<-as.numeric(gsub(pattern = ",", ".", data_cvs[,j]))}
  data_cvs_long<-melt(data_cvs[,c(1,4:ncol(data_cvs))], id.vars="sample_type")
  data_cvs_long$Plate<-dochterplaat
  
  ## Extract standard curve data
  if (length(which(tmp[,1]=="eS9"))>0){
    expected_low_concentration1<-tmp[which(tmp[,1]=="eS9"),][1,]
    expected_low_concentration2<-tmp[which(tmp[,1]=="eS8"),][1,]
    expected_low_concentration3<-tmp[which(tmp[,1]=="eS7"),][1,]
    expected_low_concentration4<-tmp[which(tmp[,1]=="eS6"),][1,]
    expected_high_concentration1<-tmp[which(tmp[,1]=="eS1"),][1,]
    expected_high_concentration2<-tmp[which(tmp[,1]=="eS2"),][1,]
    expected_high_concentration3<-tmp[which(tmp[,1]=="eS3"),][1,]}else{
      expected_low_concentration1<-tmp[which(tmp[,1]=="S9"),][1,]
      expected_low_concentration2<-tmp[which(tmp[,1]=="S8"),][1,]
      expected_low_concentration3<-tmp[which(tmp[,1]=="S7"),][1,]
      expected_low_concentration4<-tmp[which(tmp[,1]=="S6"),][1,]
      expected_high_concentration1<-tmp[which(tmp[,1]=="S1"),][1,]
      expected_high_concentration2<-tmp[which(tmp[,1]=="S2"),][1,]
      expected_high_concentration3<-tmp[which(tmp[,1]=="S3"),][1,]}
  
  expected_low_concentration <- rbind(expected_low_concentration1, expected_low_concentration2, expected_low_concentration3, expected_low_concentration4)
  expected_low_concentration[] <- lapply(expected_low_concentration, as.character)
  expected_low_concentration[,-1] <- lapply(expected_low_concentration[,-1], function(x) as.numeric(gsub(pattern = ",", ".", x)))
  expected_low_concentration1 <- as.data.frame(lapply(expected_low_concentration, min, na.rm = T))
  
  expected_low_concentration1 <- do.call(data.frame, lapply(expected_low_concentration1, function(x) replace(x, is.infinite(x),NA)))
  
  expected_high_concentration <- rbind(expected_high_concentration1, expected_high_concentration2, expected_high_concentration3)
  expected_high_concentration[] <- lapply(expected_high_concentration, as.character)
  expected_high_concentration[,-1] <- lapply(expected_high_concentration[,-1], function(x) as.numeric(gsub(pattern = ",", ".", x)))
  expected_high_concentration1 <- as.data.frame(lapply(expected_high_concentration, max, na.rm = T))
  expected_high_concentration1 <- do.call(data.frame, lapply(expected_high_concentration1, function(x) replace(x, is.infinite(x),NA)))
  
  expected_low_concentration1 <- sub("^[^[:alnum:]]", "",as.character(expected_low_concentration1))
  expected_high_concentration1 <- sub("^[^[:alnum:]]", "",as.character(expected_high_concentration1))
  colnames<-as.character(tmp[7,])
  
  
  if (length(which(tmp[,1]=="eS9"))>0){
    expected_low_concentrationA<-tmp[which(tmp[,1]=="eS9"),][1,]
    expected_low_concentrationB<-tmp[which(tmp[,1]=="eS8"),][1,]
    expected_low_concentrationC<-tmp[which(tmp[,1]=="eS7"),][1,]
    expected_low_concentrationD<-tmp[which(tmp[,1]=="eS6"),][1,]
    expected_low_concentrationE<-tmp[which(tmp[,1]=="eS5"),][1,]
    expected_low_concentrationF<-tmp[which(tmp[,1]=="eS4"),][1,]
    expected_high_concentrationG<-tmp[which(tmp[,1]=="eS1"),][1,]
    expected_high_concentrationH<-tmp[which(tmp[,1]=="eS2"),][1,]
    expected_high_concentrationI<-tmp[which(tmp[,1]=="eS3"),][1,]}else{
    expected_low_concentrationA<-tmp[which(tmp[,1]=="S9"),][1,]
    expected_low_concentrationB<-tmp[which(tmp[,1]=="S8"),][1,]
    expected_low_concentrationC<-tmp[which(tmp[,1]=="S7"),][1,]
    expected_low_concentrationD<-tmp[which(tmp[,1]=="S6"),][1,]
    expected_low_concentrationE<-tmp[which(tmp[,1]=="S5"),][1,]
    expected_low_concentrationF<-tmp[which(tmp[,1]=="S4"),][1,]
    expected_high_concentrationG<-tmp[which(tmp[,1]=="S1"),][1,]
    expected_high_concentrationH<-tmp[which(tmp[,1]=="S2"),][1,]
    expected_high_concentrationI<-tmp[which(tmp[,1]=="S3"),][1,]}
  
  Standard_curve_values_save <- rbind(expected_low_concentrationA, expected_low_concentrationB,
                                      expected_low_concentrationC, expected_low_concentrationD,
                                      expected_low_concentrationE, expected_low_concentrationF,
                                      expected_high_concentrationG, expected_high_concentrationH,
                                      expected_high_concentrationI)

  colnames(Standard_curve_values_save) <- colnames
  Standard_curve_values_save <-Standard_curve_values_save[,-c(2:3)]
  colnames(Standard_curve_values_save)[1] <- "Standard"
  Standard_curve_values_save <- Standard_curve_values_save[order(Standard_curve_values_save$Standard,decreasing=F),]
  Standard_curve_values_save[,2:14] <- lapply(Standard_curve_values_save[,2:14], function(x) as.numeric(gsub(pattern = ",", ".", x)))
  Standard_curve_values_save$Plate <- i
  
  colnames(expected_low_concentration) <- colnames
  expected_low_concentration <- melt(expected_low_concentration)
  expected_low_concentration$variable <- as.character(expected_low_concentration$variable)
  expected_low_concentration$id <- 1:nrow(expected_low_concentration)
  lowest_standard_curve_values<-data.frame(cbind(colnames, expected_low_concentration1))[c(3:length(colnames)),]
  lowest_standard_curve_values <- merge(expected_low_concentration, lowest_standard_curve_values, 
                                        by.x = c("variable", "value"), by.y = c("colnames", "expected_low_concentration1"))
  lowest_standard_curve_values <- lowest_standard_curve_values[order(lowest_standard_curve_values$id), ]
  lowest_standard_curve_values$value <- NULL
  
  colnames(expected_high_concentration) <- colnames
  expected_high_concentration <- melt(expected_high_concentration)
  expected_high_concentration$variable <- as.character(expected_high_concentration$variable)
  expected_high_concentration$id <- 1:nrow(expected_high_concentration)
  highest_standard_curve_values<-data.frame(cbind(colnames, expected_high_concentration1))[c(3:length(colnames)),]
  highest_standard_curve_values <- merge(expected_high_concentration, highest_standard_curve_values, 
                                         by.x = c("variable", "value"), by.y = c("colnames", "expected_high_concentration1"))
  highest_standard_curve_values <- highest_standard_curve_values[order(highest_standard_curve_values$id), ]
  highest_standard_curve_values$value <- NULL
  
  ## Extract expected standard values in long format
  if (length(which(tmp_expected[,1]=="eS9"))>0){
    standard_expected <- data.frame(tmp_expected[c(which(tmp_expected[,1]=="eS1")[1]:(grep("eS9", as.matrix(tmp_expected[,1])))),])}else{
      standard_expected <- data.frame(tmp_expected[c(which(tmp_expected[,1]=="S1")[1]:(grep("S9", as.matrix(tmp_expected[,1])))),])}
  colnames(standard_expected) <- colnames
  standard_expected[,2:3] <- NULL
  standard_expected[] <- lapply(standard_expected, as.character)
  standard_expected[,-1] <- lapply(standard_expected[,-1], function(x) as.numeric(gsub(pattern = ",", ".", x)))
  
  standard_expected_save <- standard_expected
  standard_expected_save$Plate <- i
  
  standard_expected <- melt(standard_expected)
  
  # Make list of LLQ and ULQ values for each variable
  lowest_standard_curve_values <- merge(standard_expected, lowest_standard_curve_values, by = c("variable", NA))
  lowest_standard_curve_values <- lowest_standard_curve_values[order(lowest_standard_curve_values$id), ]
  lowest_standard_curve_values$id <- NULL
  lowest_standard_curve_values[,2] <- NULL
  lowest_standard_curve_values[,2]<-as.numeric(gsub(pattern = ",", ".", lowest_standard_curve_values[,2]))
  lowest_standard_curve_values[,3]<-dochterplaat
  
  dilution_factor<-as.numeric(gsub(pattern = ",", ".", data.frame(tmp_dilution[which(tmp_dilution[,1]== "X10"),4])[1,1])) #NB: There should be at least 10 Unknown samples
  colnames(lowest_standard_curve_values)<-c("variable","LLQ", "Plate_number")
  lowest_standard_curve_values$LLQ<-lowest_standard_curve_values$LLQ*dilution_factor
  
  highest_standard_curve_values <- merge(standard_expected, highest_standard_curve_values, by = c("variable", NA))
  highest_standard_curve_values <- highest_standard_curve_values[order(highest_standard_curve_values$id), ]
  highest_standard_curve_values$id <- NULL
  highest_standard_curve_values[,2] <- NULL
  highest_standard_curve_values[,2]<-as.numeric(gsub(pattern = ",", ".", highest_standard_curve_values[,2]))
  highest_standard_curve_values[,3]<-dochterplaat
  colnames(highest_standard_curve_values)<-c("variable","ULQ", "Plate_number")
  highest_standard_curve_values$ULQ<-highest_standard_curve_values$ULQ*dilution_factor 
  
  ## Dataframe for concentration extraction. 
  data<-data.frame(tmp[c(which(tmp[,1]=="X1")[2]:(grep("Sampling Errors", as.matrix(tmp[,1]))-2)),])
  data[,1]<-dochterplaat
  
  ## Loop to: 
  ## 1. Put everything below measurable concentration on the lowest expected value of the standard curve (so at LLQ)
  ## 2. Same for *** 
  ## 3. Remove all stars and other characters so that everything is a number. 
  
  for (j in 3:ncol(data)){
    data[,j] <- unlist(ifelse(data[,j]=="OOR <", lowest_standard_curve_values[j-3,2], data[,j]))
    data[,j] <- unlist(ifelse(data[,j]=="OOR >", highest_standard_curve_values[j-3,2], data[,j]))
    data[,j] <- unlist(ifelse(data[,j]=="***", lowest_standard_curve_values[j-3,2], data[,j]))
    #data[,j] <- sub("^[^[:alnum:]]", "",as.character(data[,j] )) #NB: This will remove * with extrapolated values, advised to NOT use this line
    data[,j] <- gsub(pattern = ",", ".", data[,j])
    #data[,j] <- ifelse(data[,j]<lowest_standard_curve_values[j-2,2], lowest_standard_curve_values[j-2,2], data[,j]) #NB: This will change extrapolated values, advies to NOT use this line
  }
  data[,3] <- NA
  colnames(data)<-c("Plate_number", "Position", colnames[c(3:length(colnames))])
  
  ### Number of beads. 
  data_beads<-data.frame(tmp_beads[c(which(tmp_beads[,1]=="X1")[2]:grep("Sampling Errors", as.matrix(tmp_beads[,1]))),])
  data_beads<-data_beads[c(grep("+X",data_beads[,1]),grep("+C",data_beads[,1])),]
  data_beads[,1]<-dochterplaat
  data_beads[,3]<- NA
  colnames(data_beads)<-c("Plate_number", "Position", colnames[c(3:length(colnames))])
  
  ## Now melt all these dataframes 
  require(reshape2)
  data_long<-reshape2::melt(data, id.vars=c("Position", "Plate_number"))
  colnames(data_long)[4]<-"Concentration"
  data_beads_long<-reshape2::melt(data_beads, id.vars=c("Position", "Plate_number"))
  colnames(data_beads_long)[4]<-"Beads"
  
  data_long2<-merge(data_long, data_beads_long, all.x=T, by=c("Plate_number","Position","variable"))
  data_long2<-merge(data_long2, lowest_standard_curve_values, all.x=T)
  data_long2<-merge(data_long2, highest_standard_curve_values, all.x=T)
  data_long2$Beads<-as.numeric(data_long2$Beads)
  
  #Add Quality check
  #Important to realise that only One Quality check can be given, so the value can be Extrapolated but also have a low bead count,
  #However, you will only see extrapolated
  data_long2$Quality_check_beads<-"Good Bead Count"
  data_long2$Quality_check_beads<-ifelse(data_long2$Beads<50, "Low Acceptable Bead Count", data_long2$Quality_check_beads)
  data_long2$Quality_check_beads<-ifelse(data_long2$Beads<25, "Low Bead Count", data_long2$Quality_check_beads)
  
  data_long2$Quality_check<-"Good quality"
  data_long2$Quality_check<-ifelse(grepl("\\*", data_long2$Concentration), "Extrapolated", data_long2$Quality_check)
  
  
  data_long2$Concentration <- as.numeric(sub("^[^[:alnum:]]", "",as.character(data_long2$Concentration)))
  
  data_long2$Quality_check <-ifelse(data_long2$Quality_check == "Extrapolated" & data_long2$Concentration > data_long2$ULQ*0.1,
                                    "Above ULQ - extrapolated based on standard curve", data_long2$Quality_check) #*0.1 as safety range
  data_long2$Quality_check <-ifelse(data_long2$Quality_check == "Extrapolated" & data_long2$Concentration < data_long2$LLQ*10,
                                    "Under LLQ - extrapolated based on standard curve", data_long2$Quality_check) #*10 as safety range
  data_long2$Quality_check<-ifelse(round(data_long2$Concentration, 4) == round(data_long2$LLQ, 4), "Under LLQ - set to min of standard curve", data_long2$Quality_check)
  data_long2$Quality_check<-ifelse(round(data_long2$Concentration, 4) == round(data_long2$ULQ, 4), "Above ULQ - set to max of standard curve", data_long2$Quality_check)
  data_long2[which(is.na(data_long2$Beads)==T),"Concentration"]<-NA
  data_long2[which(is.na(data_long2$Beads)==T),"Quality_check"]<-"No Beads Counted"
  data_long2[which(is.na(data_long2$Beads)==T),"Quality_check_beads"]<-"No Beads Counted"
  
  #Add IDs/Plate_layout and save files
  data_long2 <- merge(data_long2, df_plate_layout, all.x=T, by = c("Plate_number","Position"))
  data_long2$Sample_material<-"Plasma"

  write.csv(data_long2, Processed_files[i], row.names = F)
  
  write.csv(data_cvs_long, CVs_files[i], row.names = F)
  
  write.csv(Standard_curve_values_save, Standard_curve_files[i], row.names = F)

  write.csv(standard_expected_save, Standard_expected_files[i], row.names = F)
  
  print(i)

}

