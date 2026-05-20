##library
library(ggplot2)
library(gridExtra)


######## paths
#input paths
input_path_1 <- snakemake@input[[1]]
input_path_2 <- snakemake@input[[2]]


input_path_4 <- snakemake@input[[3]]
input_path_5 <- snakemake@input[[4]]


#output paths
output_path_1 <- snakemake@output[[1]]
output_path_2 <- snakemake@output[[2]]


Plate_1_standard_curve <- read.csv(input_path_1)
Plate_2_standard_curve <- read.csv(input_path_2)


Plate_1_standard_curve == Plate_2_standard_curve


##Standard curves are identical between plates

Target_vec <- colnames(Plate_1_standard_curve)[2:14]

Standard_curve_plot_list <- list()

for(i in 1:length(Target_vec)){
  
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  Standard_plot_tmp <- ggplot(Plate_1_standard_curve, aes_string(x = "Standard", y = Target_vec[i], group = 1)) + 
     geom_point() + geom_line() + theme_bw() + scale_y_continuous(trans = "log", labels=scaleFUN) +
      ggtitle(Target_vec[i]) + ylab("log scaled concentration")
  
  Standard_curve_plot_list[[i]] <- Standard_plot_tmp
}

Standard_plots_arrange <- do.call(grid.arrange, c(Standard_curve_plot_list,
                                                  ncol = 6, top = "Luminex standard curve (log)"))

ggsave(Standard_plots_arrange, file = output_path_1,
       height = 16, width = 13)

###############################################################################
### Expected standard curves

Plate_1_Expected <- read.csv(input_path_4)
Plate_2_Expected <- read.csv(input_path_5)





##Expected Standard curves are identical between plates

colnames(Plate_1_Expected)[1] <- "Standard"

Target_vec <- colnames(Plate_1_Expected)[2:47]

Expected_curve_plot_list <- list()

for(i in 1:length(Target_vec)){
  
  scaleFUN <- function(x) sprintf("%.0f", x)
  
  Expected_plot_tmp <- ggplot(Plate_1_Expected, aes_string(x = "Standard", y = Target_vec[i], group = 1)) + 
    geom_point() + geom_line() + theme_bw() + scale_y_continuous(trans = "log", labels=scaleFUN) +
    ggtitle(Target_vec[i]) + ylab("log scaled concentration")
  
  Expected_curve_plot_list[[i]] <- Expected_plot_tmp
}

Expected_plots_arrange <- do.call(grid.arrange, c(Expected_curve_plot_list,
                                                  ncol = 6, top = "Luminex expected standard curve (log)"))

ggsave(Expected_plots_arrange, file = output_path_2,
       height = 16, width = 13)

