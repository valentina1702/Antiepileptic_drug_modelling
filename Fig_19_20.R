


##For Visualizations of missed dose analysis
#load libraries
library(R.matlab)
library(ggplot2)
library(tidyr)
library(ggpubr)


###PART A:  All 5 patients together

# Load the data that was saved from MATLAB
c_p_values_delayed <- readMat('c_p_values_delayed_pop_old.mat')
c_p_values_missed <- readMat('c_p_values_missed_pop_old.mat')
c_p_normal <- readMat('c_p_normal_pop_old.mat')
conc_normal <- readMat('conc_normal_pop_old.mat')
conc_delayed <- readMat('conc_delayed_pop_old.mat')
conc_missed <- readMat('conc_missed_pop_old.mat')
crcl <- readMat("crcl_values.mat")

# Convert from readMat() to a data frame
c_p_values_delayed<- as.data.frame(c_p_values_delayed)
c_p_values_missed<- as.data.frame(c_p_values_missed)
c_p_normal <- as.data.frame(c_p_normal)
conc_normal <- as.data.frame(conc_normal)
conc_delayed<- as.data.frame(conc_delayed)
conc_missed <- as.data.frame(conc_missed)
crcl <- as.data.frame(crcl)


#rename columns
names(c_p_values_delayed)[1] <- 'ctrough_m_5_1'
names(c_p_values_delayed)[2] <- 'ctrough_m_5_2'
names(c_p_values_delayed)[3] <- 'ctrough_m_5_3'
names(c_p_values_delayed)[4] <- 'ctrough_m_5_4'
names(c_p_values_delayed)[5] <- 'ctrough_m_5_5'
names(c_p_values_delayed)[6] <- 'ctrough_2m_5_1'
names(c_p_values_delayed)[7] <- 'ctrough_2m_5_2'
names(c_p_values_delayed)[8] <- 'ctrough_2m_5_3'
names(c_p_values_delayed)[9] <- 'ctrough_2m_5_4'
names(c_p_values_delayed)[10] <- 'ctrough_2m_5_5'
names(c_p_values_delayed)[11] <- 'ctrough_3m_5_1'
names(c_p_values_delayed)[12] <- 'ctrough_3m_5_2'
names(c_p_values_delayed)[13] <- 'ctrough_3m_5_3'
names(c_p_values_delayed)[14] <- 'ctrough_3m_5_4'
names(c_p_values_delayed)[15] <- 'ctrough_3m_5_5'
names(c_p_values_delayed)[16] <- 'ctrough_4m_5_1'
names(c_p_values_delayed)[17] <- 'ctrough_4m_5_2'
names(c_p_values_delayed)[18] <- 'ctrough_4m_5_3'
names(c_p_values_delayed)[19] <- 'ctrough_4m_5_4'
names(c_p_values_delayed)[20] <- 'ctrough_4m_5_5'


names(c_p_values_delayed)[21] <- 'peak_m_5_1'
names(c_p_values_delayed)[22] <- 'peak_m_5_2'
names(c_p_values_delayed)[23] <- 'peak_m_5_3'
names(c_p_values_delayed)[24] <- 'peak_m_5_4'
names(c_p_values_delayed)[25] <- 'peak_m_5_5'
names(c_p_values_delayed)[26] <- 'peak_2m_5_1'
names(c_p_values_delayed)[27] <- 'peak_2m_5_2'
names(c_p_values_delayed)[28] <- 'peak_2m_5_3'
names(c_p_values_delayed)[29] <- 'peak_2m_5_4'
names(c_p_values_delayed)[30] <- 'peak_2m_5_5'
names(c_p_values_delayed)[31] <- 'peak_3m_5_1'
names(c_p_values_delayed)[32] <- 'peak_3m_5_2'
names(c_p_values_delayed)[33] <- 'peak_3m_5_3'
names(c_p_values_delayed)[34] <- 'peak_3m_5_4'
names(c_p_values_delayed)[35] <- 'peak_3m_5_5'
names(c_p_values_delayed)[36] <- 'peak_4m_5_1'
names(c_p_values_delayed)[37] <- 'peak_4m_5_2'
names(c_p_values_delayed)[38] <- 'peak_4m_5_3'
names(c_p_values_delayed)[39] <- 'peak_4m_5_4'
names(c_p_values_delayed)[40] <- 'peak_4m_5_5'

names(c_p_values_delayed)[41] <- 'ctimes_m_5'
names(c_p_values_delayed)[42] <- 'ctimes_2m_5'
names(c_p_values_delayed)[43] <- 'ctimes_3m_5'
names(c_p_values_delayed)[44] <- 'ctimes_4m_5'
names(c_p_values_delayed)[45] <- 'ptimes_m_5_1'
names(c_p_values_delayed)[46] <- 'ptimes_m_5_2'
names(c_p_values_delayed)[47] <- 'ptimes_m_5_3'
names(c_p_values_delayed)[48] <- 'ptimes_m_5_4'
names(c_p_values_delayed)[49] <- 'ptimes_m_5_5'
names(c_p_values_delayed)[50] <- 'ptimes_2m_5_1'
names(c_p_values_delayed)[51] <- 'ptimes_2m_5_2'
names(c_p_values_delayed)[52] <- 'ptimes_2m_5_3'
names(c_p_values_delayed)[53] <- 'ptimes_2m_5_4'
names(c_p_values_delayed)[54] <- 'ptimes_2m_5_5'
names(c_p_values_delayed)[55] <- 'ptimes_3m_5_1'
names(c_p_values_delayed)[56] <- 'ptimes_3m_5_2'
names(c_p_values_delayed)[57] <- 'ptimes_3m_5_3'
names(c_p_values_delayed)[58] <- 'ptimes_3m_5_4'
names(c_p_values_delayed)[59] <- 'ptimes_3m_5_5'
names(c_p_values_delayed)[60] <- 'ptimes_4m_5_1'
names(c_p_values_delayed)[61] <- 'ptimes_4m_5_2'
names(c_p_values_delayed)[62] <- 'ptimes_4m_5_3'
names(c_p_values_delayed)[63] <- 'ptimes_4m_5_4'
names(c_p_values_delayed)[64] <- 'ptimes_4m_5_5'


names(conc_delayed)[1] <- 'time_m_5'
names(conc_delayed)[2] <- 'time_2m_5'
names(conc_delayed)[3] <- 'time_3m_5'
names(conc_delayed)[4] <- 'time_4m_5'
names(conc_delayed)[5] <- 'conc_m_5_1'
names(conc_delayed)[6] <- 'conc_m_5_2'
names(conc_delayed)[7] <- 'conc_m_5_3'
names(conc_delayed)[8] <- 'conc_m_5_4'
names(conc_delayed)[9] <- 'conc_m_5_5'
names(conc_delayed)[10] <- 'conc_2m_5_1'
names(conc_delayed)[11] <- 'conc_2m_5_2'
names(conc_delayed)[12] <- 'conc_2m_5_3'
names(conc_delayed)[13] <- 'conc_2m_5_4'
names(conc_delayed)[14] <- 'conc_2m_5_5'
names(conc_delayed)[15] <- 'conc_3m_5_1'
names(conc_delayed)[16] <- 'conc_3m_5_2'
names(conc_delayed)[17] <- 'conc_3m_5_3'
names(conc_delayed)[18] <- 'conc_3m_5_4'
names(conc_delayed)[19] <- 'conc_3m_5_5'
names(conc_delayed)[20] <- 'conc_4m_5_1'
names(conc_delayed)[21] <- 'conc_4m_5_2'
names(conc_delayed)[22] <- 'conc_4m_5_3'
names(conc_delayed)[23] <- 'conc_4m_5_4'
names(conc_delayed)[24] <- 'conc_4m_5_5'


#create theme
my_theme <- theme_classic() +
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(color = "#E6AB02"),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))



Plot1 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.1, color
                                   = "normal"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_1, color
                                    = "missed by m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_1, color
                                    = "missed by 2m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_1, color
                                    = "missed by 3m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_1, color
                                    = "missed by 4m/5"), size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.1, color
                                   = "missed completely"), size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.1, color
                                   = "normal"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.1, y = Pk.c.1, color
                                   = "normal"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_1, color
                                           = "missed by m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_1, y = peak_m_5_1, color
                                           = "missed by m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_1, color
                                           = "missed by 2m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_1, y = peak_2m_5_1, color
                                           = "missed by 2m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_1, color
                                           = "missed by 3m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_1, y = peak_3m_5_1, color
                                           = "missed by 3m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_1, color
                                           = "missed by 4m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_1, y = peak_4m_5_1, color
                                           = "missed by 4m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.1, y = Ct.c5.1, color
                                          = "missed completely"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.1, y = Pk.c5.1, color
                                          = "missed completely"),shape = 20, size = 3, stroke = 1.5) +
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.2, color
                                   = "normal"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_2, color
                                    = "missed by m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_2, color
                                    = "missed by 2m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_2, color
                                    = "missed by 3m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_2, color
                                    = "missed by 4m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.2, color
                                   = "missed completely"),alpha=0.8, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.2, color
                                   = "normal"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.2, y = Pk.c.2, color
                                   = "normal"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_2, color
                                           = "missed by m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_2, y = peak_m_5_2, color
                                           = "missed by m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_2, color
                                           = "missed by 2m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_2, y = peak_2m_5_2, color
                                           = "missed by 2m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_2, color
                                           = "missed by 3m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_2, y = peak_3m_5_2, color
                                           = "missed by 3m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_2, color
                                           = "missed by 4m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_2, y = peak_4m_5_2, color
                                           = "missed by 4m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.2, y = Ct.c5.2, color
                                          = "missed completely"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.2, y = Pk.c5.2, color
                                          = "missed completely"),alpha=0.8,shape = 20, size = 3, stroke = 1.5)+
  
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.3, color
                                   = "normal"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_3, color
                                    = "missed by m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_3, color
                                    = "missed by 2m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_3, color
                                    = "missed by 3m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_3, color
                                    = "missed by 4m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.3, color
                                   = "missed completely"),alpha=0.6, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.3, color
                                   = "normal"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.3, y = Pk.c.3, color
                                   = "normal"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_3, color
                                           = "missed by m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_3, y = peak_m_5_3, color
                                           = "missed by m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_3, color
                                           = "missed by 2m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_3, y = peak_2m_5_3, color
                                           = "missed by 2m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_3, color
                                           = "missed by 3m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_3, y = peak_3m_5_3, color
                                           = "missed by 3m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_3, color
                                           = "missed by 4m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_3, y = peak_4m_5_3, color
                                           = "missed by 4m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.3, y = Ct.c5.3, color
                                          = "missed completely"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.3, y = Pk.c5.3, color
                                          = "missed completely"),alpha=0.6,shape = 20, size = 3, stroke = 1.5)+
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.4, color
                                   = "normal"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_4, color
                                    = "missed by m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_4, color
                                    = "missed by 2m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_4, color
                                    = "missed by 3m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_4, color
                                    = "missed by 4m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.4, color
                                   = "missed completely"),alpha=0.4, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.4, color
                                   = "normal"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.4, y = Pk.c.4, color
                                   = "normal"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_4, color
                                           = "missed by m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_4, y = peak_m_5_4, color
                                           = "missed by m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_4, color
                                           = "missed by 2m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_4, y = peak_2m_5_4, color
                                           = "missed by 2m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_4, color
                                           = "missed by 3m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_4, y = peak_3m_5_4, color
                                           = "missed by 3m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_4, color
                                           = "missed by 4m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_4, y = peak_4m_5_4, color
                                           = "missed by 4m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.4, y = Ct.c5.4, color
                                          = "missed completely"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.4, y = Pk.c5.4, color
                                          = "missed completely"),alpha=0.4,shape = 20, size = 3, stroke = 1.5)+
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.5, color
                                   = "normal"),alpha=0.2, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_5, color
                                    = "missed by m/5"),alpha=0.2, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_5, color
                                    = "missed by 2m/5"),alpha=0.2, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_5, color
                                    = "missed by 3m/5"),alpha=0.2, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_5, color
                                    = "missed by 4m/5"),alpha=0.2, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.5, color
                                   = "missed completely"),alpha=0.2, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.5, color
                                   = "normal"),alpha=0.2,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.5, y = Pk.c.5, color
                                   = "normal"),alpha=0.2,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_5, color
                                           = "missed by m/5"),alpha=0.2,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_5, y = peak_m_5_5, color
                                           = "missed by m/5"),alpha=0.2,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_5, color
                                           = "missed by 2m/5"),alpha=0.2,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_5, y = peak_2m_5_5, color
                                           = "missed by 2m/5"),alpha=0.2,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_5, color
                                           = "missed by 3m/5"),alpha=0.2,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_5, y = peak_3m_5_5, color
                                           = "missed by 3m/5"),alpha=0.2,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_5, color
                                           = "missed by 4m/5"),alpha=0.2,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_5, y = peak_4m_5_5, color
                                           = "missed by 4m/5"),alpha=0.2,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.5, y = Ct.c5.5, color
                                          = "missed completely"),alpha=0.2,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.5, y = Pk.c5.5, color
                                          = "missed completely"),alpha=0.2,shape = 20, size = 3, stroke = 1.5)+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) +
  labs(title = 'Missed Dose Analysis for varying population') +
  my_theme +
  theme(legend.position = "bottom")
ggsave(filename = 'missed_dose_pop.png', plot = Plot1, width = 12, height = 10)

###PART B:  First 3 patients together
Plot2 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.1, color
                                   = "normal"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_1, color
                                    = "missed by m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_1, color
                                    = "missed by 2m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_1, color
                                    = "missed by 3m/5"), size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_1, color
                                    = "missed by 4m/5"), size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.1, color
                                   = "missed completely"), size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.1, color
                                   = "normal"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.1, y = Pk.c.1, color
                                   = "normal"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_1, color
                                           = "missed by m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_1, y = peak_m_5_1, color
                                           = "missed by m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_1, color
                                           = "missed by 2m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_1, y = peak_2m_5_1, color
                                           = "missed by 2m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_1, color
                                           = "missed by 3m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_1, y = peak_3m_5_1, color
                                           = "missed by 3m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_1, color
                                           = "missed by 4m/5"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_1, y = peak_4m_5_1, color
                                           = "missed by 4m/5"),shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.1, y = Ct.c5.1, color
                                          = "missed completely"),shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.1, y = Pk.c5.1, color
                                          = "missed completely"),shape = 20, size = 3, stroke = 1.5) +
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.2, color
                                   = "normal"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_2, color
                                    = "missed by m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_2, color
                                    = "missed by 2m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_2, color
                                    = "missed by 3m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_2, color
                                    = "missed by 4m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.2, color
                                   = "missed completely"),alpha=0.6, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.2, color
                                   = "normal"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.2, y = Pk.c.2, color
                                   = "normal"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_2, color
                                           = "missed by m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_2, y = peak_m_5_2, color
                                           = "missed by m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_2, color
                                           = "missed by 2m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_2, y = peak_2m_5_2, color
                                           = "missed by 2m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_2, color
                                           = "missed by 3m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_2, y = peak_3m_5_2, color
                                           = "missed by 3m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_2, color
                                           = "missed by 4m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_2, y = peak_4m_5_2, color
                                           = "missed by 4m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.2, y = Ct.c5.2, color
                                          = "missed completely"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.2, y = Pk.c5.2, color
                                          = "missed completely"),alpha=0.6,shape = 20, size = 3, stroke = 1.5)+
  
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.3, color
                                   = "normal"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_3, color
                                    = "missed by m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_3, color
                                    = "missed by 2m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_3, color
                                    = "missed by 3m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_3, color
                                    = "missed by 4m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.3, color
                                   = "missed completely"),alpha=0.4, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.3, color
                                   = "normal"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.3, y = Pk.c.3, color
                                   = "normal"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_3, color
                                           = "missed by m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_3, y = peak_m_5_3, color
                                           = "missed by m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_3, color
                                           = "missed by 2m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_3, y = peak_2m_5_3, color
                                           = "missed by 2m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_3, color
                                           = "missed by 3m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_3, y = peak_3m_5_3, color
                                           = "missed by 3m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_3, color
                                           = "missed by 4m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_3, y = peak_4m_5_3, color
                                           = "missed by 4m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.3, y = Ct.c5.3, color
                                          = "missed completely"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.3, y = Pk.c5.3, color
                                          = "missed completely"),alpha=0.4,shape = 20, size = 3, stroke = 1.5)+
  
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) +
  labs(title = 'Missed Dose Analysis for varying population') +
  my_theme + annotate("text",label="Creatine clearance for Patient 1 =  82.07301", color="#66A61E", x=40, y=80, alpha=1,fontface=2,size=8)+
  annotate("text",label="Creatine clearance for Patient 2= 162.61268", color="#66A61E", x=40, y=76, alpha=0.6,fontface=2,size=8)+
  annotate("text",label="Creatine clearance for Patient 3= 173.81625", color="#66A61E", x=40, y=72, alpha=0.4,fontface=2,size=8)+
  theme(legend.position = "bottom")
ggsave(filename = 'fig_19.png', plot = Plot2, width = 12, height = 10)

###PART c:  last 3 patients together

Plot3 <- ggplot(data=NULL) +
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.3, color
                                   = "normal"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_3, color
                                    = "missed by m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_3, color
                                    = "missed by 2m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_3, color
                                    = "missed by 3m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_3, color
                                    = "missed by 4m/5"),alpha=0.8, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.3, color
                                   = "missed completely"),alpha=0.8, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.3, color
                                   = "normal"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.3, y = Pk.c.3, color
                                   = "normal"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_3, color
                                           = "missed by m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_3, y = peak_m_5_3, color
                                           = "missed by m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_3, color
                                           = "missed by 2m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_3, y = peak_2m_5_3, color
                                           = "missed by 2m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_3, color
                                           = "missed by 3m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_3, y = peak_3m_5_3, color
                                           = "missed by 3m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_3, color
                                           = "missed by 4m/5"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_3, y = peak_4m_5_3, color
                                           = "missed by 4m/5"),alpha=0.8,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.3, y = Ct.c5.3, color
                                          = "missed completely"),alpha=0.8,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.3, y = Pk.c5.3, color
                                          = "missed completely"),alpha=0.8,shape = 20, size = 3, stroke = 1.5)+
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.4, color
                                   = "normal"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_4, color
                                    = "missed by m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_4, color
                                    = "missed by 2m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_4, color
                                    = "missed by 3m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_4, color
                                    = "missed by 4m/5"),alpha=0.6, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.4, color
                                   = "missed completely"),alpha=0.6, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.4, color
                                   = "normal"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.4, y = Pk.c.4, color
                                   = "normal"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_4, color
                                           = "missed by m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_4, y = peak_m_5_4, color
                                           = "missed by m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_4, color
                                           = "missed by 2m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_4, y = peak_2m_5_4, color
                                           = "missed by 2m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_4, color
                                           = "missed by 3m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_4, y = peak_3m_5_4, color
                                           = "missed by 3m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_4, color
                                           = "missed by 4m/5"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_4, y = peak_4m_5_4, color
                                           = "missed by 4m/5"),alpha=0.6,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.4, y = Ct.c5.4, color
                                          = "missed completely"),alpha=0.6,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.4, y = Pk.c5.4, color
                                          = "missed completely"),alpha=0.6,shape = 20, size = 3, stroke = 1.5)+
  
  geom_line(data = conc_normal,aes(x = Tc, y = Y.c.5, color
                                   = "normal"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_m_5, y = conc_m_5_5, color
                                    = "missed by m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_2m_5, y = conc_2m_5_5, color
                                    = "missed by 2m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_3m_5, y = conc_3m_5_5, color
                                    = "missed by 3m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_delayed,aes(x = time_4m_5, y = conc_4m_5_5, color
                                    = "missed by 4m/5"),alpha=0.4, size = 1.5) +
  geom_line(data = conc_missed,aes(x = Tcm, y = Y.cm.5, color
                                   = "missed completely"),alpha=0.4, size = 1.5) +
  geom_point(data = c_p_normal,aes(x = DoseTimes, y = Ct.c.5, color
                                   = "normal"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_normal,aes(x = times.5, y = Pk.c.5, color
                                   = "normal"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_m_5, y = ctrough_m_5_5, color
                                           = "missed by m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_m_5_5, y = peak_m_5_5, color
                                           = "missed by m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_2m_5, y = ctrough_2m_5_5, color
                                           = "missed by 2m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_2m_5_5, y = peak_2m_5_5, color
                                           = "missed by 2m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_3m_5, y = ctrough_3m_5_5, color
                                           = "missed by 3m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_3m_5_5, y = peak_3m_5_5, color
                                           = "missed by 3m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ctimes_4m_5, y = ctrough_4m_5_5, color
                                           = "missed by 4m/5"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_delayed,aes(x = ptimes_4m_5_5, y = peak_4m_5_5, color
                                           = "missed by 4m/5"),alpha=0.4,shape = 20, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = missed.DoseT.5, y = Ct.c5.5, color
                                          = "missed completely"),alpha=0.4,shape = 1, size = 3, stroke = 1.5) +
  geom_point(data = c_p_values_missed,aes(x = timesm.5, y = Pk.c5.5, color
                                          = "missed completely"),alpha=0.4,shape = 20, size = 3, stroke = 1.5)+
  geom_hline(yintercept=46,linetype='dotted', col = 'red')+
  geom_hline(yintercept=12,linetype='dotted', col = 'red')+
  scale_color_manual(name ='',values = c("#1B9E77" ,"#D95F02" ,"#7570B3" ,"#E7298A" ,"#66A61E","#E6AB02"),
                     limits = c("missed by m/5", "missed by 2m/5", "missed by 3m/5","missed by 4m/5", "missed completely", "normal"))+
  #scale_color_brewer(name = '', palette = 'Dark2') +
  scale_x_continuous(name = 'Time (hr)', breaks = seq(0,1300,10)) +
  scale_y_continuous(name = '[D] (mg/L)', limits = c(0,80), breaks =
                       seq(0,80,10)) +
  labs(title = 'Missed Dose Analysis for varying population') +
  my_theme + annotate("text",label="Creatine clearance for Patient 3= 173.81625", color="#66A61E", x=40, y=80, alpha=0.8,fontface=2,size=8)+
  annotate("text",label="Creatine clearance for Patient 4 =  97.60801", color="#66A61E", x=40, y=76, alpha=0.6,fontface=2,size=8)+
  annotate("text",label="Creatine clearance for Patient 5= 120.70736", color="#66A61E", x=40, y=72, alpha=0.4,fontface=2,size=8)+
  theme(legend.position = "bottom")
ggsave(filename = 'Fig_20.png', plot = Plot3, width = 12, height = 10)

