#######################################################################################################
# PACKAGES

# Load necessary packages
library('R.matlab') 
library('tidyverse') 
library('patchwork')
library('plotly')

# Be sure to change your working directory to the folder where the files are!

#######################################################################################################
# IMPORT DATA FROM MATLAB
heatmap_theme <- theme_minimal() +
  theme(text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.key.size=unit(0.5,'inch'),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.title.x = element_text(size=20, face='bold',margin = margin(10,0,0,0)),
        axis.title.y = element_text(size=20, face='bold',margin = margin(0,10,0,0)),
        axis.text= element_text(size = 20))

# SensTarr_Flag1
case1_SensTarr <- readMat('1_1.mat')
case1_SensTarr <- as.data.frame(case1_SensTarr)
case1_ST_1 <- data.frame(case1_SensTarr$SensTarr.1, case1_SensTarr$t.2)
case1_ST_1$var <- c('ka')
colnames(case1_ST_1) <- c('SensTarr','time','Parameter')

case1_ST_2 <- data.frame(case1_SensTarr$SensTarr.2, case1_SensTarr$t.3)
case1_ST_2$var <- c('kc')
colnames(case1_ST_2) <- c('SensTarr','time','Parameter')

case1_ST_3 <- data.frame(case1_SensTarr$SensTarr.3, case1_SensTarr$t.4)
case1_ST_3$var <- c('V')
colnames(case1_ST_3) <- c('SensTarr','time','Parameter')

case1_ST_total<- rbind(case1_ST_1,case1_ST_2,case1_ST_3)

case1_SensTarr_plot<-ggplot(data=case1_ST_total, aes(x=time, y=SensTarr, fill=Parameter)) +
  geom_line(size=2,aes(color=Parameter))+ylab('Relative sensitivity of concentration(dY/Y)/(dP/P)')+xlab('Time (hrs)')+
  theme(axis.text.x = element_text(size = 20, vjust = 1.5), 
        axis.text.y = element_text(size = 20, vjust = 1.5), 
        axis.title=element_text(size=18,face="bold"), legend.text=element_text(size=18))
png('case1_SensTarr_plot.png')
print(case1_SensTarr_plot)
dev.off()

# Sens_Flag1
case1_Sens <- readMat('1_2.mat')
case1_Sens <- as.data.frame(case1_Sens)
case1_Sens_1 <- data.frame(case1_Sens$Sens.1)
case1_Sens_1$var <- c('ka')
colnames(case1_Sens_1) <- c('Sens','Parameter')

case1_Sens_2 <- data.frame(case1_Sens$Sens.2)
case1_Sens_2$var <- c('kc')
colnames(case1_Sens_2) <- c('Sens','Parameter')

case1_Sens_3 <- data.frame(case1_Sens$Sens.3)
case1_Sens_3$var <- c('V')
colnames(case1_Sens_3) <- c('Sens','Parameter')

case1_Sens_total<- rbind(case1_Sens_1,case1_Sens_2,case1_Sens_3)
case1_Sens_plot<-ggplot(data=case1_Sens_total, aes(x=Parameter, y=Sens, fill=Parameter)) +
  geom_bar(stat="identity")+ xlab('Parameter') + ylab('Relative sensitivity of AUC (dY/Y)/(dP/P)')+
  theme(axis.text.x = element_text(size = 20, vjust = 1.5), 
        axis.text.y = element_text(size = 20, vjust = 1.5), 
        axis.title=element_text(size=20,face="bold"), legend.text=element_text(size=18))


png('case1_Sens_plot.png')
print(case1_Sens_plot)
dev.off()

# SensB_Flag1
case1_SensB <- readMat('1_3.mat')
case1_SensB <- as.data.frame(case1_SensB)
case1_SensB_1 <- data.frame(case1_SensB$SensB.1)
case1_SensB_1$var <- c('ka')
colnames(case1_SensB_1) <- c('SensB','Parameter')

case1_SensB_2 <- data.frame(case1_SensB$SensB.2)
case1_SensB_2$var <- c('kc')
colnames(case1_SensB_2) <- c('SensB','Parameter')

case1_SensB_3 <- data.frame(case1_SensB$SensB.3)
case1_SensB_3$var <- c('V')
colnames(case1_SensB_3) <- c('SensB','Parameter')

case1_SensB_total<- rbind(case1_SensB_1,case1_SensB_2,case1_SensB_3)
case1_SensB_plot<-ggplot(data=case1_SensB_total, aes(x=Parameter, y=SensB, fill=Parameter)) +
  geom_bar(stat="identity") + xlab('Parameters') + ylab('Absolute sensitivity of AUC(dY/dP)')+
  theme(axis.text.x = element_text(size = 20, vjust = 1.5), 
        axis.text.y = element_text(size = 20, vjust = 1.5), 
        axis.title=element_text(size=20,face="bold"),  legend.text=element_text(size=18))

png('case1_Sens_plot.png')
print(case1_SensB_total)
dev.off()

# SensTMax_Flag1
case1_SensTmax <- readMat('1_4.mat')
case1_SensTmax <- as.data.frame(case1_SensTmax)
case1_SensTmax_1 <- data.frame(case1_SensTmax$SensTmax.1)
case1_SensTmax_1$var <- c('ka')
colnames(case1_SensTmax_1) <- c('SensTmax','Parameter')

case1_SensTmax_2 <- data.frame(case1_SensTmax$SensTmax.2)
case1_SensTmax_2$var <- c('kc')
colnames(case1_SensTmax_2) <- c('SensTmax','Parameter')

case1_SensTmax_3 <- data.frame(case1_SensTmax$SensTmax.3)
case1_SensTmax_3$var <- c('V')
colnames(case1_SensTmax_3) <- c('SensTmax','Parameter')

case1_SensTmax_total<- rbind(case1_SensTmax_1,case1_SensTmax_2,case1_SensTmax_3)
case1_SensTmax_plot<-ggplot(data=case1_SensTmax_total, aes(x=Parameter, y=SensTmax, fill=Parameter)) +
  geom_bar(stat="identity")+ xlab('Parameters') + ylab('Time of maximum sensitivity of concentration (hrs)')+
  theme(axis.text.x = element_text(size = 20, vjust = 1.5), 
        axis.text.y = element_text(size = 20, vjust = 1.5), 
        axis.title=element_text(size=18,face="bold"),  legend.text=element_text(size=18))


png('case1_SensTmax_plot.png')
print(case1_SensTmax_plot)
dev.off()

# Sens_Flag2
case2_Sens <- readMat('2_1.mat')
case2_Sens <- as.data.frame(case2_Sens)
parameter <- c('ka', 'kc', 'V')
var <- c('ka', 'kc', 'V')
case2_Sens <- cbind(parameter, case2_Sens)
names(case2_Sens)<-c('Parameter', var)
case2_Sens <- pivot_longer(case2_Sens, !Parameter, names_to = 'var', values_to = 'Sens')
case2_Sens$Parameter <- factor(case2_Sens$Parameter, levels = c('ka', 'kc', 'V'))

case2_Sens_plot<-ggplot(case2_Sens, aes(var, Parameter, fill= Sens)) + 
  geom_tile() + heatmap_theme+ 
  scale_fill_gradient('Sensitivity', limits=c(-1, 0.1), breaks = c(-1, -0.75, -0.5, -0.25, 0),  low = 'lightblue', high = "black")+xlab('Parameter') +
  ylab('Parameter')+ ggtitle("Local sensitivity") +
  theme(axis.text.x = element_text(size = 20, vjust = 1.5), 
        axis.text.y = element_text(size = 20, vjust = 1.5), axis.title=element_text(size=18,face="bold"))

png('case2_Sens_plot.png')
print(case2_Sens_plot)
dev.off()

# SensTmax_Flag2
case2_SensTmax <- readMat('2_2.mat')
case2_SensTmax <- as.data.frame(case2_SensTmax)
parameter <- c('ka', 'kc', 'V')
var <- c('ka', 'kc', 'V')
case2_SensTmax <- cbind(parameter, case2_SensTmax)
names(case2_SensTmax)<-c('Parameter', var)
case2_SensTmax <- pivot_longer(case2_SensTmax, !Parameter, names_to = 'var', values_to = 'Sens')
case2_SensTmax$Parameter <- factor(case2_SensTmax$Parameter, levels = c('ka', 'kc', 'V'))

case2_SensTmax_plot<-ggplot(case2_SensTmax, aes(var, Parameter, fill= Sens)) + 
  geom_tile()  + heatmap_theme+ 
  xlab('Parameter') + ylab('Parameter') + ggtitle("Time of peak sensitivity (hr)") +
  scale_fill_gradient('Sensitivity', limits=c(0, 120), breaks = c(0, 30, 60, 90, 120),  low = 'lightblue', high = "black")+xlab('Parameter') +
  theme(axis.text.x = element_text(size = 20, vjust = 1.5), axis.text.y = element_text(size = 20, vjust = 1.5), axis.title=element_text(size=20,face="bold"))

png('case2_SensTmax_plot.png')
print(case2_SensTmax_plot)
dev.off()

# SensSyn_Flag2
case2_SensSyn <- readMat('2_3.mat')
case2_SensSyn <- as.data.frame(case2_SensSyn)
parameter <- c('ka', 'kc', 'V')
var <- c('ka', 'kc', 'V')
case2_SensSyn <- cbind(parameter, case2_SensSyn)
names(case2_SensSyn)<-c('Parameter', var)
case2_SensSyn <- pivot_longer(case2_SensSyn, !Parameter, names_to = 'var', values_to = 'Sens')
case2_SensSyn$Parameter <- factor(case2_SensSyn$Parameter, levels = c('ka', 'kc', 'V'))

case2_SensSyn_plot<-ggplot(case2_SensSyn, aes(var, Parameter, fill= Sens)) + 
  geom_tile()  + heatmap_theme+ xlab('Parameter') + ylab('Parameter')+ 
  ggtitle("Synergy (Sij/(Si*Sj))") 
png('case2_SensSyn_plot.png')
print(case2_SensSyn_plot)
dev.off()

# Sens_Flag3
case3_Sens <- readMat('3_1.mat')
case3_Sens <- as.data.frame(case3_Sens)
parameter <- c('ka', 'kc', 'V')
doses<- c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500')
doses  <- factor(doses, levels=c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500'))

case3_Sens <- cbind(doses, case3_Sens)
names(case3_Sens)<-c('Parameter', parameter)
print(case3_Sens)
case3_Sens <- pivot_longer(case3_Sens, !Parameter, names_to = 'doses', values_to = 'Sens')

case3_Sens_plot<-ggplot(case3_Sens, aes(doses, Parameter, fill= Sens)) + 
  geom_tile()  + heatmap_theme+ ylab('Dose (mg)') + xlab('Parameter')+ 
  ggtitle("AUC Local (parameters) & Global (dose) sensitivity") +
  scale_fill_gradient('Sensitivity', limits=c(-1.1, 0.1), breaks = c(-1, -0.75, -0.5, -0.25, 0),  low = 'lightblue', high = "black")+xlab('Parameter') +
  
png('case3_Sens_plot.png')
print(case3_Sens_plot)
dev.off()

# SensTmax_Flag3
case3_SensTmax <- readMat('3_2.mat')
case3_SensTmax <- as.data.frame(case3_SensTmax)
parameter <- c('ka', 'kc', 'V')

doses<- c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500')
doses  <- factor(doses, levels=c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500'))

case3_SensTmax <- cbind(doses, case3_SensTmax)
names(case3_SensTmax)<-c('Parameter', parameter)
print(case3_SensTmax)
case3_SensTmax <- pivot_longer(case3_SensTmax, !Parameter, names_to = 'doses', values_to = 'Sens')

case3_SensTmax_plot<-ggplot(case3_SensTmax, aes(doses, Parameter, fill= Sens)) + 
  geom_tile()  + heatmap_theme+ ylab('Dose (mg)') + xlab('Parameter')+ 
  ggtitle("Time of peak sensitivity (hr)") +
  scale_fill_gradient('Sensitivity', limits=c(0, 125), breaks = c(0, 30, 60, 90, 120),  low = 'black', high = "lightblue")+xlab('Parameter') +
  
png('case3_SensTmax_plot.png')
print(case3_SensTmax_plot)
dev.off()

# Sens_Flag4
case4_Sens <- readMat('4_1.mat')
case4_Sens <- as.data.frame(case4_Sens)
parameter <- c('ka', 'kc', 'V')
doses<- c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500')
doses  <- factor(doses, levels=c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500'))
case4_Sens <- cbind(doses, case4_Sens)
names(case4_Sens)<-c('Parameter', parameter)
print(case4_Sens)
case4_Sens <- data.frame(case4_Sens$Parameter, case4_Sens$ka, case4_Sens$kc, case4_Sens$V)
colnames(case4_Sens) <- c('Parameter','ka','kc', 'V')

case4_Sens <- pivot_longer(case4_Sens, !Parameter, names_to = 'doses', values_to = 'Sens')

case4_Sens_plot<-ggplot(case4_Sens, aes(doses, Parameter, fill= Sens)) + 
  geom_tile()  + heatmap_theme+ ylab('Dose (mg)') + xlab('Parameter')+ ggtitle("Local (parameters) & Global (drug concentration) sensitivity") 
png('case4_Sens_plot.png')
print(case4_Sens_plot)
dev.off()
# SensTmax_Flag4
case4_SensTmax <- readMat('4_2.mat')
case4_SensTmax <- as.data.frame(case4_SensTmax)
parameter <- c('ka', 'kc', 'V')
doses<- c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500')
doses  <- factor(doses, levels=c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500'))
case4_SensTmax <- cbind(doses, case4_SensTmax)
names(case4_SensTmax)<-c('Parameter', parameter)
print(case4_SensTmax)
case4_SensTmax <- data.frame(case4_SensTmax$Parameter, case4_SensTmax$ka, case4_SensTmax$kc, case4_SensTmax$V)
colnames(case4_SensTmax) <- c('Parameter','ka','kc', 'V')

case4_SensTmax <- pivot_longer(case4_SensTmax, !Parameter, names_to = 'doses', values_to = 'Sens')

case4_SensTmax_plot<-ggplot(case4_SensTmax, aes(doses, Parameter, fill= Sens)) + 
  geom_tile() + heatmap_theme+ ylab('Drug concentration (mg)') + xlab('Parameter')+ ggtitle("Time of peak sensitivity (hr)") +
  scale_fill_gradient('Sensitivity', limits=c(0, 0.6), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),  low = 'black', high = "lightblue")+xlab('Parameter') +
  
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 1.5))


png('case4_SensTmax_plot.png')
print(case4_SensTmax_plot)
dev.off()

# Ctrough_Flag4
case4_SensTmax <- readMat('4_3.mat')
case4_SensTmax <- as.data.frame(case4_SensTmax)
parameter <- c('ka', 'kc', 'V')
doses<- c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500')
doses  <- factor(doses, levels=c('100', '367', '633', '900', '1167', '1433', '1700', '1967', '2233', '2500'))
case4_SensTmax <- cbind(doses, case4_SensTmax)
names(case4_SensTmax)<-c('Parameter', parameter)
print(case4_SensTmax)
case4_SensTmax <- data.frame(case4_SensTmax$Parameter, case4_SensTmax$ka, case4_SensTmax$kc, case4_SensTmax$V)
colnames(case4_SensTmax) <- c('Parameter','ka','kc', 'V')

case4_SensTmax <- pivot_longer(case4_SensTmax, !Parameter, names_to = 'doses', values_to = 'Ctrough')

case4_SensTmax_plot<-ggplot(case4_SensTmax, aes(doses, Parameter, fill= Ctrough)) + 
  geom_tile() + heatmap_theme+ ylab('Dose (mg)') + xlab('Parameter')+ 
  scale_fill_gradient('Ctrough', limits=c(0, 0.6), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),  low = 'black', high = "lightblue")+xlab('Parameter') +
  
  ggtitle("Time of peak sensitivity (hr)") 
png('case4_Ctrough_plot.png')
print(case4_SensTmax_plot)
dev.off()

# AUC_Flag5
case5_AUC <- readMat('globalAUC.mat')
case5_AUC <- as.data.frame(case5_AUC)
case5_Param <- readMat('globalParam.mat')
case5_Param <- as.data.frame(case5_Param)

Vfinal<- case5_Param$Vfinal
Dosesfinal <- case5_Param$DosesFinal
case5_AUC <- cbind(Vfinal, case5_AUC)
names(case5_AUC) <- c('Vfinal',Dosesfinal)
case5_AUC <- pivot_longer(case5_AUC, !Vfinal, names_to = 'Doses', values_to = 'AUC')
case5_AUC$Doses <- factor(case5_AUC$Doses, levels = unique(case5_AUC$Doses), ordered = TRUE)

globalSensPlot <- ggplot(data = case5_AUC, aes(x = reorder(Doses, Doses), y = Vfinal, fill = AUC)) +
  geom_tile() +
  xlab(label = "Doses(mg)") +
  ylab(label = "Volume of distribution (L)") +
  scale_fill_viridis_c(limits = c(0, 8000)) +
  labs(title = 'AUC Global Sensitivity') +
  heatmap_theme + 
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 1.5))
ggsave('GlobalSensitivityDoseVol.png')
# Ctrough_Flag5
case5_Ctrough <- readMat('globalctrough.mat')
case5_Ctrough <- as.data.frame(case5_Ctrough)
case5_Param <- readMat('globalParam.mat')
case5_Param <- as.data.frame(case5_Param)

Vfinal<- case5_Param$Vfinal
Dosesfinal <- case5_Param$DosesFinal
case5_Ctrough <- cbind(Vfinal, case5_Ctrough)
names(case5_Ctrough) <- c('Vfinal',Dosesfinal)
case5_Ctrough <- pivot_longer(case5_Ctrough, !Vfinal, names_to = 'Doses', values_to = 'Trough')
case5_Ctrough$Doses <- factor(case5_Ctrough$Doses, levels = unique(case5_Ctrough$Doses), ordered = TRUE)

globalTroughPlot <- ggplot(data = case5_Ctrough, aes(x = reorder(Doses, Doses), y = Vfinal, fill = Trough)) +
  geom_tile() +
  xlab(label = "Doses(mg)") +
  ylab(label = "Volume of distribution (L)") +
  scale_fill_viridis_c(limits = c(0, 20)) +
  labs(title = 'AUC Global Sensitivity') +
  heatmap_theme + 
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 1.5))
