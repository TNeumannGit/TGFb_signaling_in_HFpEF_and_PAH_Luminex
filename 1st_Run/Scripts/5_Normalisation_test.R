#library
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

####### paths
#input paths
input_path_1 <- snakemake@input[[1]]

#output paths
output_path_1 <- snakemake@output[[1]]

###read data in
DF <- read.csv(input_path_1) 
DF$ID_number <- NULL

##write function for normality plots
norm_plot_fun <- function(data_in){
  
  norm_plot_list <- list()
  
  for(i in 1:length(data_in)){
    
    tmp_string <- colnames(data_in)[i]
    tmp_sub_df <- data_in[,tmp_string]
    tmp_sub_df <- na.omit(tmp_sub_df)
    
    norm_plot_tmp <- ggplot(data_in, aes_string(tmp_string)) +
      stat_bin(aes(y=..density..), breaks = seq(min(tmp_sub_df), max(tmp_sub_df), 
                                                by = .1), color="black", fill = "skyblue2")  +
      geom_line(stat="density", size = 0.7) + theme_bw() + 
      ggtitle(tmp_string) + xlab("Abundance")
    
    
    norm_plot_list[[i]] <- norm_plot_tmp
    names(norm_plot_list)[[i]] <- tmp_string
  }
  return(norm_plot_list)
}


Norm_plot_list <- norm_plot_fun(DF)

Norm_plots <- do.call(grid.arrange, c(Norm_plot_list,
                                      ncol = 5, top = "Histograms log transformed concentrations"))

ggsave(Norm_plots, file =  output_path_1,
       height = 18, width = 18)



