


##For Visualizations of missed dose analysis
#load libraries
library(R.matlab)
library(ggplot2)
library(tidyr)
library(ggpubr)


###PART A:  plots for 500 mg

# Load the data that was saved from MATLAB
c_p_values_delayed <- readMat('c_p_values_delayed_500.mat')
c_p_values_missed <- readMat('c_p_values_missed_500.mat')
c_p_normal <- readMat('c_p_normal_500.mat')
conc_normal <- readMat('conc_normal_500.mat')
conc_delayed <- readMat('conc_delayed_500.mat')
conc_missed <- readMat('conc_missed_500.mat')


# Convert from readMat() to a data frame
c_p_values_delayed<- as.data.frame(c_p_values_delayed)
c_p_values_missed<- as.data.frame(c_p_values_missed)
c_p_normal <- as.data.frame(c_p_normal)
conc_normal <- as.data.frame(conc_normal)
conc_delayed<- as.data.frame(conc_delayed)
conc_missed <- as.data.frame(conc_missed)



#rename columns
names(c_p_values_delayed)[1] <- 'ctimes_m_5'
names(c_p_values_delayed)[2] <- 'ctimes_2m_5'
names(c_p_values_delayed)[3] <- 'ctimes_3m_5'
names(c_p_values_delayed)[4] <- 'ctimes_4m_5'
names(c_p_values_delayed)[5] <- 'ctrough_m_5'
names(c_p_values_delayed)[6] <- 'ctrough_2m_5'
names(c_p_values_delayed)[7] <- 'ctrough_3m_5'
names(c_p_values_delayed)[8] <- 'ctrough_4m_5'
names(c_p_values_delayed)[9] <- 'ptimes_m_5'
names(c_p_values_delayed)[10] <- 'ptimes_2m_5'
names(c_p_values_delayed)[11] <- 'ptimes_3m_5'
names(c_p_values_delayed)[12] <- 'ptimes_4m_5'
names(c_p_values_delayed)[13] <- 'peak_m_5'
names(c_p_values_delayed)[14] <- 'peak_2m_5'
names(c_p_values_delayed)[15] <- 'peak_3m_5'
names(c_p_values_delayed)[16] <- 'peak_4m_5'
names(conc_delayed)[1] <- 'time_m_5'
names(conc_delayed)[2] <- 'time_2m_5'
names(conc_delayed)[3] <- 'time_3m_5'
names(conc_delayed)[4] <- 'time_4m_5'
names(conc_delayed)[5] <- 'conc_m_5'
names(conc_delayed)[6] <- 'conc_2m_5'
names(conc_delayed)[7] <- 'conc_3m_5'
names(conc_delayed)[8] <- 'conc_4m_5'


#create theme
my_theme <- theme_classic() +
  theme(text = element_text(size = 30),legend.title = element_text(size=30),
        legend.key.height= unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))

mPlot1 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = T, y = y, color
                                  = "normal"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5, color
                                   = "missed by m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5, color
                                   = "missed by 2m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5, color
                                   = "missed by 3m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5, color
                                   = "missed by 4m/5"), size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tm, y = ym, color
                                   = "missed completely"), size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ctrough, color
                             = "normal"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times, y = peak, color
                              = "normal"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5, color
                                   = "missed by m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5, y = peak_m_5, color
                                   = "missed by m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5, color
                                   = "missed by 2m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5, y = peak_2m_5, color
                                   = "missed by 2m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5, color
                                   = "missed by 3m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5, y = peak_3m_5, color
                                   = "missed by 3m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5, color
                                   = "missed by 4m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5, y = peak_4m_5, color
                                   = "missed by 4m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT, y = Ctroughm, color
                                   = "missed completely"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm, y = peakm, color
                                   = "missed completely"),shape = 20, size = 3, stroke = 1.5) +
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) +
  labs(title = '500 mg',tag="A")+
  my_theme+ theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold",size = 30))
ggsave(filename = 'missed_dose_500.png', plot = mPlot1, width = 8, height = 6)


###PART B:  plots for 1000 mg

# Load the data that was saved from MATLAB
c_p_values_delayed <- readMat('c_p_values_delayed_1000.mat')
c_p_values_missed <- readMat('c_p_values_missed_1000.mat')
c_p_normal <- readMat('c_p_normal_1000.mat')
conc_normal <- readMat('conc_normal_1000.mat')
conc_delayed <- readMat('conc_delayed_1000.mat')
conc_missed <- readMat('conc_missed_1000.mat')


# Convert from readMat() to a data frame
c_p_values_delayed<- as.data.frame(c_p_values_delayed)
c_p_values_missed<- as.data.frame(c_p_values_missed)
c_p_normal <- as.data.frame(c_p_normal)
conc_normal <- as.data.frame(conc_normal)
conc_delayed<- as.data.frame(conc_delayed)
conc_missed <- as.data.frame(conc_missed)



#rename columns
names(c_p_values_delayed)[1] <- 'ctimes_m_5'
names(c_p_values_delayed)[2] <- 'ctimes_2m_5'
names(c_p_values_delayed)[3] <- 'ctimes_3m_5'
names(c_p_values_delayed)[4] <- 'ctimes_4m_5'
names(c_p_values_delayed)[5] <- 'ctrough_m_5'
names(c_p_values_delayed)[6] <- 'ctrough_2m_5'
names(c_p_values_delayed)[7] <- 'ctrough_3m_5'
names(c_p_values_delayed)[8] <- 'ctrough_4m_5'
names(c_p_values_delayed)[9] <- 'ptimes_m_5'
names(c_p_values_delayed)[10] <- 'ptimes_2m_5'
names(c_p_values_delayed)[11] <- 'ptimes_3m_5'
names(c_p_values_delayed)[12] <- 'ptimes_4m_5'
names(c_p_values_delayed)[13] <- 'peak_m_5'
names(c_p_values_delayed)[14] <- 'peak_2m_5'
names(c_p_values_delayed)[15] <- 'peak_3m_5'
names(c_p_values_delayed)[16] <- 'peak_4m_5'
names(conc_delayed)[1] <- 'time_m_5'
names(conc_delayed)[2] <- 'time_2m_5'
names(conc_delayed)[3] <- 'time_3m_5'
names(conc_delayed)[4] <- 'time_4m_5'
names(conc_delayed)[5] <- 'conc_m_5'
names(conc_delayed)[6] <- 'conc_2m_5'
names(conc_delayed)[7] <- 'conc_3m_5'
names(conc_delayed)[8] <- 'conc_4m_5'



mPlot2 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = T, y = y, color
                                   = "normal"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5, color
                                    = "missed by m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5, color
                                    = "missed by 2m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5, color
                                    = "missed by 3m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5, color
                                    = "missed by 4m/5"), size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tm, y = ym, color
                                   = "missed completely"), size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ctrough, color
                                   = "normal"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times, y = peak, color
                                   = "normal"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5, color
                                           = "missed by m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5, y = peak_m_5, color
                                           = "missed by m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5, color
                                           = "missed by 2m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5, y = peak_2m_5, color
                                           = "missed by 2m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5, color
                                           = "missed by 3m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5, y = peak_3m_5, color
                                           = "missed by 3m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5, color
                                           = "missed by 4m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5, y = peak_4m_5, color
                                           = "missed by 4m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT, y = Ctroughm, color
                                          = "missed completely"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm, y = peakm, color
                                          = "missed completely"),shape = 20, size = 3, stroke = 1.5) +
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) +
  labs(title = '1000 mg',tag="C")+
  my_theme+ theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold",size = 30))
ggsave(filename = 'missed_dose_1000.png', plot = mPlot2, width = 8, height = 6)


###PART C:  plots for 1500 mg

# Load the data that was saved from MATLAB
c_p_values_delayed <- readMat('c_p_values_delayed_1500.mat')
c_p_values_missed <- readMat('c_p_values_missed_1500.mat')
c_p_normal <- readMat('c_p_normal_1500.mat')
conc_normal <- readMat('conc_normal_1500.mat')
conc_delayed <- readMat('conc_delayed_1500.mat')
conc_missed <- readMat('conc_missed_1500.mat')


# Convert from readMat() to a data frame
c_p_values_delayed<- as.data.frame(c_p_values_delayed)
c_p_values_missed<- as.data.frame(c_p_values_missed)
c_p_normal <- as.data.frame(c_p_normal)
conc_normal <- as.data.frame(conc_normal)
conc_delayed<- as.data.frame(conc_delayed)
conc_missed <- as.data.frame(conc_missed)



#rename columns
names(c_p_values_delayed)[1] <- 'ctimes_m_5'
names(c_p_values_delayed)[2] <- 'ctimes_2m_5'
names(c_p_values_delayed)[3] <- 'ctimes_3m_5'
names(c_p_values_delayed)[4] <- 'ctimes_4m_5'
names(c_p_values_delayed)[5] <- 'ctrough_m_5'
names(c_p_values_delayed)[6] <- 'ctrough_2m_5'
names(c_p_values_delayed)[7] <- 'ctrough_3m_5'
names(c_p_values_delayed)[8] <- 'ctrough_4m_5'
names(c_p_values_delayed)[9] <- 'ptimes_m_5'
names(c_p_values_delayed)[10] <- 'ptimes_2m_5'
names(c_p_values_delayed)[11] <- 'ptimes_3m_5'
names(c_p_values_delayed)[12] <- 'ptimes_4m_5'
names(c_p_values_delayed)[13] <- 'peak_m_5'
names(c_p_values_delayed)[14] <- 'peak_2m_5'
names(c_p_values_delayed)[15] <- 'peak_3m_5'
names(c_p_values_delayed)[16] <- 'peak_4m_5'
names(conc_delayed)[1] <- 'time_m_5'
names(conc_delayed)[2] <- 'time_2m_5'
names(conc_delayed)[3] <- 'time_3m_5'
names(conc_delayed)[4] <- 'time_4m_5'
names(conc_delayed)[5] <- 'conc_m_5'
names(conc_delayed)[6] <- 'conc_2m_5'
names(conc_delayed)[7] <- 'conc_3m_5'
names(conc_delayed)[8] <- 'conc_4m_5'


#create theme


mPlot3 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = T, y = y, color
                                   = "normal"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5, color
                                    = "missed by m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5, color
                                    = "missed by 2m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5, color
                                    = "missed by 3m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5, color
                                    = "missed by 4m/5"), size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tm, y = ym, color
                                   = "missed completely"), size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ctrough, color
                                   = "normal"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times, y = peak, color
                                   = "normal"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5, color
                                           = "missed by m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5, y = peak_m_5, color
                                           = "missed by m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5, color
                                           = "missed by 2m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5, y = peak_2m_5, color
                                           = "missed by 2m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5, color
                                           = "missed by 3m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5, y = peak_3m_5, color
                                           = "missed by 3m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5, color
                                           = "missed by 4m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5, y = peak_4m_5, color
                                           = "missed by 4m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT, y = Ctroughm, color
                                          = "missed completely"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm, y = peakm, color
                                          = "missed completely"),shape = 20, size = 3, stroke = 1.5) +
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) + 
  labs(title = '1500 mg',tag="D")+
  my_theme+ theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold",size = 30))
ggsave(filename = 'missed_dose_1500.png', plot = mPlot3, width = 8, height = 6)

###PART D:  plots for 750 mg

# Load the data that was saved from MATLAB
c_p_values_delayed <- readMat('c_p_values_delayed_750.mat')
c_p_values_missed <- readMat('c_p_values_missed_750.mat')
c_p_normal <- readMat('c_p_normal_750.mat')
conc_normal <- readMat('conc_normal_750.mat')
conc_delayed <- readMat('conc_delayed_750.mat')
conc_missed <- readMat('conc_missed_750.mat')


# Convert from readMat() to a data frame
c_p_values_delayed<- as.data.frame(c_p_values_delayed)
c_p_values_missed<- as.data.frame(c_p_values_missed)
c_p_normal <- as.data.frame(c_p_normal)
conc_normal <- as.data.frame(conc_normal)
conc_delayed<- as.data.frame(conc_delayed)
conc_missed <- as.data.frame(conc_missed)



#rename columns
names(c_p_values_delayed)[1] <- 'ctimes_m_5'
names(c_p_values_delayed)[2] <- 'ctimes_2m_5'
names(c_p_values_delayed)[3] <- 'ctimes_3m_5'
names(c_p_values_delayed)[4] <- 'ctimes_4m_5'
names(c_p_values_delayed)[5] <- 'ctrough_m_5'
names(c_p_values_delayed)[6] <- 'ctrough_2m_5'
names(c_p_values_delayed)[7] <- 'ctrough_3m_5'
names(c_p_values_delayed)[8] <- 'ctrough_4m_5'
names(c_p_values_delayed)[9] <- 'ptimes_m_5'
names(c_p_values_delayed)[10] <- 'ptimes_2m_5'
names(c_p_values_delayed)[11] <- 'ptimes_3m_5'
names(c_p_values_delayed)[12] <- 'ptimes_4m_5'
names(c_p_values_delayed)[13] <- 'peak_m_5'
names(c_p_values_delayed)[14] <- 'peak_2m_5'
names(c_p_values_delayed)[15] <- 'peak_3m_5'
names(c_p_values_delayed)[16] <- 'peak_4m_5'
names(conc_delayed)[1] <- 'time_m_5'
names(conc_delayed)[2] <- 'time_2m_5'
names(conc_delayed)[3] <- 'time_3m_5'
names(conc_delayed)[4] <- 'time_4m_5'
names(conc_delayed)[5] <- 'conc_m_5'
names(conc_delayed)[6] <- 'conc_2m_5'
names(conc_delayed)[7] <- 'conc_3m_5'
names(conc_delayed)[8] <- 'conc_4m_5'


#create theme


mPlot4 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = T, y = y, color
                                   = "normal"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5, color
                                    = "missed by m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5, color
                                    = "missed by 2m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5, color
                                    = "missed by 3m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5, color
                                    = "missed by 4m/5"), size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tm, y = ym, color
                                   = "missed completely"), size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ctrough, color
                                   = "normal"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times, y = peak, color
                                   = "normal"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5, color
                                           = "missed by m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5, y = peak_m_5, color
                                           = "missed by m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5, color
                                           = "missed by 2m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5, y = peak_2m_5, color
                                           = "missed by 2m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5, color
                                           = "missed by 3m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5, y = peak_3m_5, color
                                           = "missed by 3m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5, color
                                           = "missed by 4m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5, y = peak_4m_5, color
                                           = "missed by 4m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT, y = Ctroughm, color
                                          = "missed completely"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm, y = peakm, color
                                          = "missed completely"),shape = 20, size = 3, stroke = 1.5) +
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) +
  labs(title = '750 mg',tag="B")+
  my_theme+ theme(legend.position = "none")+
  theme(plot.title = element_text(face="bold",size = 30))
ggsave(filename = 'missed_dose_750.png', plot = mPlot4, width = 8, height = 6)

legend_b <- get_legend(mPlot4 + theme(legend.position = "bottom",legend.text=element_text(size=30)))

#combine all plots
p <- ggarrange(mPlot1, mPlot2,mPlot4, mPlot3, widths = c(7,7), heights = c(8,8),nrow = 2 ,ncol= 2)
title = " Missed Dose Analysis for Concentration in Blood"
p <- ggarrange(p, legend_b, nrow = 2, ncol = 1, heights = c(16,2))
p <- annotate_figure(p,top = text_grob(title,color = "Black",face = "bold",size = 35))
ggsave(filename = 'fig16.png', plot = p, width = 20, height = 20)

