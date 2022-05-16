


##For Visualizations of missed dose analysis for relative percentage changes in AUC for set of different patients
#load libraries
library(R.matlab)
library(ggplot2)
library(tidyr)
library(ggpubr)

my_theme <- theme_classic() +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))


# Load the data that was saved from MATLAB
AUc_values50  <- readMat('change_means50.mat')
AUc_values5  <- readMat('change_means5.mat')
AUc_values10  <- readMat('change_means10.mat')
AUc_values75  <- readMat('change_means75.mat')
AUc_values100  <- readMat('change_means100.mat')
# Convert from readMat() to a data frame
AUc_values5 <- as.data.frame(AUc_values5)
AUc_values10 <- as.data.frame(AUc_values10)
AUc_values50<- as.data.frame(AUc_values50)
AUc_values75 <- as.data.frame(AUc_values75)
AUc_values100 <- as.data.frame(AUc_values100)
#reformat data into  a dataframe
df <- data.frame(values = c(AUc_values5$change.mean.1,AUc_values10$change.mean.1,AUc_values50$change.mean.1,AUc_values75$change.mean.1,AUc_values100$change.mean.1,
                            AUc_values5$change.mean.2,AUc_values10$change.mean.2,AUc_values50$change.mean.2,AUc_values75$change.mean.2,AUc_values100$change.mean.2,
                            AUc_values5$change.mean.3,AUc_values10$change.mean.3,AUc_values50$change.mean.3,AUc_values75$change.mean.3,AUc_values100$change.mean.3,
                            AUc_values5$change.mean.4,AUc_values10$change.mean.4,AUc_values50$change.mean.4,AUc_values75$change.mean.4,AUc_values100$change.mean.4,
                            AUc_values5$change.mean.5,AUc_values10$change.mean.5,AUc_values50$change.mean.5,AUc_values75$change.mean.5,AUc_values100$change.mean.5),
                 dosage_regimen=c("missed by m/5","missed by m/5","missed by m/5","missed by m/5","missed by m/5",
                                  "missed by 2m/5","missed by 2m/5","missed by 2m/5","missed by 2m/5","missed by 2m/5",
                                  "missed by 3m/5","missed by 3m/5","missed by 3m/5","missed by 3m/5","missed by 3m/5",
                                  "missed by 4m/5","missed by 4m/5","missed by 4m/5","missed by 4m/5","missed by 4m/5",
                                  "missed completely","missed completely","missed completely","missed completely","missed completely"),
                 no_of_subjects=c(5,10,50,75,100,5,10,50,75,100,5,10,50,75,100,5,10,50,75,100,5,10,50,75,100))

df$dosage_regimen <- factor(df$dosage_regimen, levels = c("missed by m/5","missed by 2m/5","missed by 3m/5","missed by 4m/5","missed completely"))


plot1 <- ggplot(df,                                      # Grouped barplot using ggplot2
                aes(fill = dosage_regimen,
                    y = values,
                    x = no_of_subjects)) +
  geom_bar(stat = "identity",
           position = "dodge")+
  scale_x_continuous(name = 'No of patients', breaks = c(seq(0,100,25),5,10)) +
  
  #scale_y_log10(name = 'Percentage change in AUC wrt normal dosage') +
  scale_fill_manual(name ='Dosage Regimen',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E"),
                    limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely"))+
  scale_y_sqrt(name = 'AUC values', limits = c(0,14), breaks =
                       seq(0,14,2)) +
  labs(title = 'Missed Dose Analysis for population variation')+ my_theme 


plot1

ggsave(filename = 'Fig_21.png', plot = plot1, width = 12, height = 6) 
