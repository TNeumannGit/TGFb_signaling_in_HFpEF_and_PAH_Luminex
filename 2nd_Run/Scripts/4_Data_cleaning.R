#library
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(ggpubr)

####### paths
#input paths
input_path_1 <- snakemake@input[[1]]

#output paths
output_path_1 <- snakemake@output[[1]]
output_path_2 <- snakemake@output[[2]]
output_path_3 <- snakemake@output[[3]]
output_path_4 <- snakemake@output[[4]]
output_path_5 <- snakemake@output[[5]]

output_path_6 <- snakemake@output[[6]]
output_path_7 <- snakemake@output[[7]]
output_path_8 <- snakemake@output[[8]]
output_path_9 <- snakemake@output[[9]]

####### read data in
df_in <- read.csv(input_path_1)

######################################## check for duplicates

ID_vec <- unique(df_in$ID_number)

duplicated_vec <- c()
for(i in 1:length(ID_vec)){
  ID_sub <- subset(df_in, ID_number == ID_vec[i])
  duplicated_vec <- c(duplicated_vec, any(duplicated(ID_sub$variable)))
}

duplicate <- ID_vec[which(duplicated_vec)]
duplicate_sub <- subset(df_in, ID_number == duplicate)
duplicated_values <- df_in[which(duplicated(df_in[,c(1, 8)])), ]

if(length(duplicate) > 0){
  df_no_dupl <- subset(df_in, ID_number != duplicate)
}else{
  df_no_dupl <- df_in
}

######################################## remove not detectable values
##write function to plot bad quality targets and patients IDs
bad_quality_plot_fun <- function(data_in){
  
  ##### Plot bead counts and remove samples with bead count under 25
  
  bead_plot_list <- list()
  
  for(i in 1:length(unique(data_in$variable))){
    
    variable_tmp <- unique(data_in$variable)[[i]]
    
    subset_tmp <- subset(data_in, variable == variable_tmp)
    
    bead_plot_tmp <- ggplot(subset_tmp, aes(x =  Beads,
                                            y = reorder(ID_number, Beads),
                                            fill = Quality_check_beads)) +
      geom_bar(stat = "identity", colour = "black") + theme_bw() + 
      ylab("") + ggtitle(sprintf("%s", variable_tmp))
    
    bead_plot_list[[i]] <- bead_plot_tmp
    
  }
  
  bead_plots_arrange <- ggarrange(plotlist = bead_plot_list,
                                  common.legend = T,
                                  legend="bottom")
  
  bead_plots_arrange <- annotate_figure(bead_plots_arrange,
                                        top = text_grob("Bead counts barplots",
                                                        size = 20))
  
  
  
  ###### Target
  variable_vec <- unique(data_in$variable)
  
  for(i in 1:length(variable_vec)){
    Quality_file_sub <- subset(data_in, variable == variable_vec[i])
    Bad_quality_sub <- subset(Quality_file_sub, Quality_check == "Under LLQ - set to min of standard curve")
    out_df_tmp <-  data.frame(ID_number = variable_vec[i],
                              Percentage_bad_quality =  (nrow(Bad_quality_sub) / nrow(Quality_file_sub) * 100))
    
    if(exists("out_df") && is.data.frame(get("out_df"))){
      out_df <- rbind(out_df, out_df_tmp)
    }else{out_df <- out_df_tmp}
  }
  
  Bad_quality_perc_target_plot <- ggplot(out_df, aes(x = reorder(ID_number, -Percentage_bad_quality), 
                                                     y =  Percentage_bad_quality)) +
    geom_bar(stat = "identity", fill = "skyblue2", colour = "black") + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    xlab("") + ylab("Percentage") + ggtitle("Percentage under LLQ between Targets")
  
  
  rm(out_df)
  ####### Participant ID
  ID_vec <- unique(data_in$ID_number)
  
  for(i in 1:length(ID_vec)){
    Quality_file_sub <- subset(data_in, ID_number == ID_vec[i])
    Bad_quality_sub <- subset(Quality_file_sub, Quality_check == "Under LLQ - set to min of standard curve")
    out_df_tmp <-  data.frame(ID_number = ID_vec[i],
                              Percentage_bad_quality =  (nrow(Bad_quality_sub) / nrow(Quality_file_sub) * 100))
    
    if(exists("out_df") && is.data.frame(get("out_df"))){
      out_df <- rbind(out_df, out_df_tmp)
    }else{out_df <- out_df_tmp}
  }
  
  
  Bad_quality_perc_ID_plot <- ggplot(out_df, aes(x =  Percentage_bad_quality,
                                                 y = reorder(ID_number, Percentage_bad_quality), 
                                                 )) +
    geom_bar(stat = "identity", fill = "skyblue2", colour = "black") + theme_bw() + 
    ylab("") + xlab("Percentage") + ggtitle("Percentage under LLQ between Participant IDs")
  
  rm(out_df)
  #### Plate number
  Plate_vec <- unique(data_in$Plate_number)
  
  for(i in 1:length(Plate_vec)){
    Plate_sub <- subset(data_in, Plate_number == Plate_vec[i])
    Bad_quality_sub <- subset(Plate_sub, Quality_check == "Under LLQ - set to min of standard curve")
    out_df_tmp <-  data.frame(Plate_number = Plate_vec[i],
                              Percentage_bad_quality =  (nrow(Bad_quality_sub) / nrow(Plate_sub) * 100))
    
    if(exists("out_df") && is.data.frame(get("out_df"))){
      out_df <- rbind(out_df, out_df_tmp)
    }else{out_df <- out_df_tmp}
  }
  
  Bad_quality_perc_Plate_plot <- ggplot(out_df, aes(x = Plate_number, y =  Percentage_bad_quality)) +
    geom_bar(stat = "identity", fill = "skyblue2", colour = "black") + theme_bw() + 
    xlab("") + ylab("Percentage") + ggtitle("Percentage under LLQ between Plates")

  
  #### Compare groups
  data_in_groups <- data_in %>%
    mutate(
      Group = "Other",
      Group = ifelse(grepl("CTR", ID_number, fixed = TRUE), "Control", Group),
      Group = ifelse(grepl("HF", ID_number, fixed = TRUE), "HFpEF", Group),
      Group = ifelse(grepl("BMP", ID_number, fixed = TRUE), "BMPR2", Group),
      Group = ifelse(grepl("PAH", ID_number, fixed = TRUE), "PAH", Group),
      Group = ifelse(grepl("J", ID_number, fixed = TRUE), "Jessie", Group),
    )
  
  Bad_quality <- subset(data_in_groups, Quality_check == "Under LLQ - set to min of standard curve")
  
  Group_bad_quality_perc_df <-  as.data.frame(table(Bad_quality$Group) / table(data_in_groups$Group) * 100)
  colnames(Group_bad_quality_perc_df) <- c("Group", "Percentage_under_LLQ")
  
  Group_num_data <- as.data.frame(table(data_in_groups$Group))
  colnames(Group_num_data) <- c("Group", "Full_number")
  
  Bad_num_data <- as.data.frame(table(Bad_quality$Group))
  colnames(Bad_num_data) <- c("Group", "Bad_quality")
  
  Group_num_data_2 <- merge(Group_num_data, Bad_num_data, by = "Group")
  
  Group_num_melt <- melt(Group_num_data_2)
  
  Percentage_plots <- ggplot(Group_bad_quality_perc_df, aes(Group, Percentage_under_LLQ)) +
    geom_col() + theme_bw() + 
    ggtitle("Under LLQ - set to min of standard curve (Percentage)")
  
  
  Number_plots <- ggplot(Group_num_melt, aes(Group, value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") + theme_bw() + 
    ggtitle("Under LLQ - set to min of standard curve (Absolute number)")
  
  Bad_quality_group_plots_allign <- grid.arrange(Percentage_plots, Number_plots)
  
  out_list <- list(bead_plots_arrange, Bad_quality_perc_target_plot, Bad_quality_perc_ID_plot, Bad_quality_perc_Plate_plot, Bad_quality_group_plots_allign) 
  names(out_list) <- c("Targets", "IDs", "Plates")
  return(out_list)
}

Bad_quality_plots <- bad_quality_plot_fun(df_no_dupl)

ggsave(Bad_quality_plots[[1]], file = output_path_1,
       height = 49, width = 49)
ggsave(Bad_quality_plots[[2]], file = output_path_2, height = 8, width = 10)
ggsave(Bad_quality_plots[[3]], file = output_path_3, height = 28, width = 10)
ggsave(Bad_quality_plots[[4]], file = output_path_4, height = 5, width = 6)
ggsave(Bad_quality_plots[[5]], file = output_path_5, height = 12, width = 14)


######## remove all entries with bead count under 25

df_bad_beads <- subset(df_no_dupl, Beads < 25)
unique(df_bad_beads$ID_number)

df_Good_Targets_and_IDs <- subset(df_no_dupl, Beads > 24)

#### No removal of of samples or targets

nrow(df_Good_Targets_and_IDs[which(df_Good_Targets_and_IDs$Quality_check == "Under LLQ - set to min of standard curve"), ])

#### All still contained under LLQ values are set to 0.00001

Under_LLQ_sub <- df_Good_Targets_and_IDs[which(df_Good_Targets_and_IDs$Quality_check == "Under LLQ - set to min of standard curve"), ] 

df_Good_Targets_and_IDs[which(
  df_Good_Targets_and_IDs$Quality_check == "Under LLQ - set to min of standard curve"),
  "Concentration"] <- 0.00001


df_out_long <- df_Good_Targets_and_IDs

#### save output table
write.csv(df_out_long, file = output_path_6,
          row.names = F)

## log transform values
df_out_long_log <- df_out_long
df_out_long_log$Concentration <- log(df_out_long_log$Concentration + 1)

#### save log transformed output table
write.csv(df_out_long_log, file = output_path_7,
          row.names = F)

############# Make table with only concentrations and sample IDs

Conc_table <- df_Good_Targets_and_IDs[,c("ID_number", "variable", "Concentration")]

Conc_table <- dcast(Conc_table, ID_number ~ variable)

##log transform data
Conc_table_log <- Conc_table
Conc_table_log[,c(2:ncol(Conc_table_log))] <- log(Conc_table_log[,c(2:ncol(Conc_table_log))] + 1)

#### save output tables
write.csv(Conc_table, file = output_path_8,
          row.names = F)

write.csv(Conc_table_log, file = output_path_9,
          row.names = F)

