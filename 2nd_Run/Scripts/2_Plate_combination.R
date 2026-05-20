

######## paths
#input paths
input_path_1 <- snakemake@input[[1]]
input_path_2 <- snakemake@input[[2]]


#output paths
output_path_1 <- snakemake@output[[1]]


Plate_1_df <- read.csv(input_path_1)
Plate_2_df <- read.csv(input_path_2)


df_combined <- rbind(Plate_1_df, Plate_2_df)
df_combined <- na.omit(df_combined)

df_combined$variable <-  sub(" [^ ]+$", "", df_combined$variable)
df_combined$variable <- sub(" ", "_", df_combined$variable)

df_combined$Position <- NULL
df_combined$Sample_material <- NULL

write.csv(df_combined, file = output_path_1,
          row.names = F)
