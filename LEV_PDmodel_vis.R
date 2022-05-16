
# load packages
library('R.matlab')
library('tidyverse')
library('patchwork')
library('plotly')
library('ggplot2')
library('ggfortify')
library('survival')
library('survminer')

# load data from MATLAB
Tesim <- as.data.frame(readMat('Tesim.mat'))
Tesim_avg <- as.data.frame(readMat('Tesim_avg.mat'))
Tesim_med <- as.data.frame(readMat('Tesim_med.mat'))
survCurvs <- as.data.frame(readMat('survCurvs.mat'))
survCurvs_fit <- as.data.frame(readMat('survCurvs_fit.mat'))

# subset data by dosing regimen
Tesim_1 <- Tesim['Tesim.T.1']
Tesim_2 <- Tesim['Tesim.T.2']
Tesim_3 <- Tesim['Tesim.T.3']
Tesim_4 <- Tesim['Tesim.T.4']
Tesim_5 <- Tesim['Tesim.T.5']
Tesim_6 <- Tesim['Tesim.T.6']

survCurv_1 <- survCurvs['survCurvs.1']
survCurv_2 <- survCurvs['survCurvs.2']
survCurv_3 <- survCurvs['survCurvs.3']
survCurv_4 <- survCurvs['survCurvs.4']
survCurv_5 <- survCurvs['survCurvs.5']
survCurv_6 <- survCurvs['survCurvs.6']

# create dosing regimen labels
l_DR1 <- rep('A: 2000 mg BID', times = nrow(Tesim_1))
l_DR2 <- rep('B: 1500 mg BID', times = nrow(Tesim_2))
l_DR3 <- rep('C: 1000 mg BID', times = nrow(Tesim_3))
l_DR4 <- rep('D: 750 mg BID',  times = nrow(Tesim_4))
l_DR5 <- rep('E: 500 mg BID',  times = nrow(Tesim_5))
l_DR6 <- rep('F: 250 mg BID',  times = nrow(Tesim_6))

# create censoring status labels
l_cens <- rep(1, times = nrow(survCurv_1))

# assign dosing regimen labels
Tesim_1 <- cbind(Tesim_1, l_DR1)
Tesim_2 <- cbind(Tesim_2, l_DR2)
Tesim_3 <- cbind(Tesim_3, l_DR3)
Tesim_4 <- cbind(Tesim_4, l_DR4)
Tesim_5 <- cbind(Tesim_5, l_DR5)
Tesim_6 <- cbind(Tesim_6, l_DR6)

survCurv_1 <- cbind(survCurv_1, l_cens, l_DR1)
survCurv_2 <- cbind(survCurv_2, l_cens, l_DR2)
survCurv_3 <- cbind(survCurv_3, l_cens, l_DR3)
survCurv_4 <- cbind(survCurv_4, l_cens, l_DR4)
survCurv_5 <- cbind(survCurv_5, l_cens, l_DR5)
survCurv_6 <- cbind(survCurv_6, l_cens, l_DR6)

# reassign column names
col_names <- c('Te','DR')
names(Tesim_1) <- col_names
names(Tesim_2) <- col_names
names(Tesim_3) <- col_names
names(Tesim_4) <- col_names
names(Tesim_5) <- col_names
names(Tesim_6) <- col_names

col_names <- c('Te','status','DR')
names(survCurv_1) <- col_names
names(survCurv_2) <- col_names
names(survCurv_3) <- col_names
names(survCurv_4) <- col_names
names(survCurv_5) <- col_names
names(survCurv_6) <- col_names

# recombine data subsets
Tesim <- rbind(Tesim_1, Tesim_2, Tesim_3, Tesim_4, Tesim_5, Tesim_6)
survCurvs <- rbind(survCurv_1, survCurv_2, survCurv_3, survCurv_4, survCurv_5, survCurv_6)

# set up universal plot theme
plot_theme <- theme_classic() +
  theme(text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, margin = margin(10,10,10,10)),
        axis.title.x = element_text(margin = margin(5,0,0,0)),
        axis.title.y = element_text(margin = margin(0,5,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)),
        legend.title = element_blank())

# visualize time to event histograms by dosing regimen (stacked)
Te_histograms <- ggplot(data = Tesim, aes(x = Te, fill = DR)) +
  geom_histogram(binwidth = 1, color = 'white', alpha = 0.6, position = 'identity') +
  coord_cartesian(xlim = c(0,14)) +
  labs(title = 'Time to seizure by dosing regimen',
       x = 'Time (months)',
       y = 'Percent') +
  plot_theme + theme(legend.position = c(0.8, 0.8))

ggsave(filename = 'Te_histograms.png', width = 5, height = 5)

# visualize Kaplan-Meier survival curves by dosing regimen (stacked)
# https://rpkgs.datanovia.com/survminer/reference/ggsurvplot_arguments.html
surv_fit <- survfit(Surv(Te, status) ~ DR, data = survCurvs)
survival_curves <- ggsurvplot(fit = surv_fit,
                              data = survCurvs,
                              conf.int = TRUE,
                              conf.int.style = 'ribbon',
                              conf.int.alpha = 0.15,
                              title = 'Survival curves by dosing regimen',
                              xlab = 'Time (months)',
                              ylab = 'Probability of seizure free survival',
                              xlim = c(0,14),
                              legend = 'right',
                              legend.title = element_blank(),
                              legend.labs = c('2000 mg BID','1500 mg BID','1000 mg BID','750 mg BID','500 mg BID','250 mg BID'),
                              ggtheme = plot_theme)

ggsave(filename = 'survival_curves.png', width = 6.5, height = 5)
