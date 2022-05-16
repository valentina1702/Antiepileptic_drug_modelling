

##For Visualizations of missed dose analysis for changes in AUC for different dosage amounts
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
AUc_values500  <- readMat('AUc_values_500.mat')
AUc_values1000  <- readMat('AUc_values_1000.mat')
AUc_values750  <- readMat('AUc_values_750.mat')
AUc_values1500  <- readMat('AUc_values_1500.mat')
# Convert from readMat() to a data frame
AUc_values500 <- as.data.frame(AUc_values500)
AUc_values750 <- as.data.frame(AUc_values750)
AUc_values1000<- as.data.frame(AUc_values1000)
AUc_values1500 <- as.data.frame(AUc_values1500)
#reformat data into  a dataframe
df <- data.frame(values = c(AUc_values500$auc.1,AUc_values750$auc.1,AUc_values1000$auc.1,AUc_values1500$auc.1,
                            AUc_values500$auc.2,AUc_values750$auc.2,AUc_values1000$auc.2,AUc_values1500$auc.2,
                            AUc_values500$auc.3,AUc_values750$auc.3,AUc_values1000$auc.3,AUc_values1500$auc.3,
                            AUc_values500$auc.4,AUc_values750$auc.4,AUc_values1000$auc.4,AUc_values1500$auc.4,
                            AUc_values500$auc.5,AUc_values750$auc.5,AUc_values1000$auc.5,AUc_values1500$auc.5,
                            AUc_values500$auc.6,AUc_values750$auc.6,AUc_values1000$auc.6,AUc_values1500$auc.6),
                 dosage_regimen=c("normal","normal","normal","normal",
                                  "missed by m/5","missed by m/5","missed by m/5","missed by m/5",
                                  "missed by 2m/5","missed by 2m/5","missed by 2m/5","missed by 2m/5",
                                  "missed by 3m/5","missed by 3m/5","missed by 3m/5","missed by 3m/5",
                                  "missed by 4m/5","missed by 4m/5","missed by 4m/5","missed by 4m/5",
                                  "missed completely","missed completely","missed completely","missed completely"),
                 dosage=c(500,750,1000,1500,500,750,1000,1500,500,750,1000,1500,500,750,1000,1500,500,750,1000,1500,500,750,1000,1500))
df$dosage_regimen <- factor(df$dosage_regimen, levels = c("normal","missed by m/5","missed by 2m/5","missed by 3m/5","missed by 4m/5","missed completely"))



plot1 <- ggplot(df,                                      # Grouped barplot using ggplot2
                aes(fill = dosage_regimen,
                    y = values,
                    x = dosage)) +
  geom_bar(stat = "identity",
           position = "dodge")+
  scale_fill_manual(name ='Dosage Regimen',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  scale_x_continuous(name = 'Different doses', breaks = c(seq(0,2000,500),750)) +
  scale_y_continuous(name = 'AUC values', limits = c(0,5000), breaks =
                       seq(0,5000,500)) +
  labs(title = 'Missed Dose Analysis different doses')+my_theme




ggsave(filename = 'Fig_18.png', plot = plot1, width = 8, height = 6) 