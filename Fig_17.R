


##For Visualizations of missed dose analysis for relative percentage changes in AUC across different dosage amounts
#load libraries
library(R.matlab)
library(ggplot2)
library(tidyr)
library(ggpubr)
my_theme <- theme_classic() +
  theme(text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))



# Load the data that was saved from MATLAB
AUc_values500  <- readMat('change_means_500.mat')
AUc_values1000  <- readMat('change_means_1000.mat')
AUc_values750  <- readMat('change_means_750.mat')
AUc_values1500  <- readMat('change_means_1500.mat')
# Convert from readMat() to a data frame
AUc_values500 <- as.data.frame(AUc_values500)
AUc_values1000 <- as.data.frame(AUc_values1000)
AUc_values1500<- as.data.frame(AUc_values1500)
AUc_values750 <- as.data.frame(AUc_values750)

#reformat data into  a dataframe
df <- data.frame(values = c(AUc_values500$change.1,AUc_values750$change.1,AUc_values1000$change.1,AUc_values1500$change.1,
                            AUc_values500$change.2,AUc_values750$change.2,AUc_values1000$change.2,AUc_values1500$change.2,
                            AUc_values500$change.3,AUc_values750$change.3,AUc_values1000$change.3,AUc_values1500$change.3,
                            AUc_values500$change.4,AUc_values750$change.4,AUc_values1000$change.4,AUc_values1500$change.4,
                            AUc_values500$change.5,AUc_values750$change.5,AUc_values1000$change.5,AUc_values1500$change.5),
                 dosage_regimen=c("missed by m/5","missed by m/5","missed by m/5","missed by m/5",
                         "missed by 2m/5","missed by 2m/5","missed by 2m/5","missed by 2m/5",
                         "missed by 3m/5","missed by 3m/5","missed by 3m/5","missed by 3m/5",
                         "missed by 4m/5","missed by 4m/5","missed by 4m/5","missed by 4m/5",
                         "missed completely","missed completely","missed completely","missed completely"),
                 dosage=c(500,750,1000,1500,500,750,1000,1500,500,750,1000,1500,500,750,1000,1500,500,750,1000,1500))
df$dosage_regimen <- factor(df$dosage_regimen, levels = c("missed by m/5","missed by 2m/5","missed by 3m/5","missed by 4m/5","missed completely"))


   
plot1 <- ggplot(df,                                      # Grouped barplot using ggplot2
       aes(fill = dosage_regimen,
           y = values,
           x = dosage)) +
  geom_bar(stat = "identity",
           position = "dodge")+scale_x_continuous(name = 'Different doses', breaks = c(seq(0,2000,500),750)) +
  #scale_y_log10(name = 'Percentage change in AUC wrt normal dosage') + +
  scale_fill_manual(name ='Dosage Regimen',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E"),
                    limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely"))+

  scale_y_sqrt(name = 'AUC values', limits = c(0,14), breaks =
                 seq(0,14,2)) +
  labs(title = 'Missed Dose Analysis for different doses')+my_theme

ggsave(filename = 'Fig_17.png', plot = plot1, width = 8, height = 6) 


