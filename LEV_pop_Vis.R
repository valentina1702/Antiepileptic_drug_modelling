# Final Project - Levitaracetam
# Visualizations for population analysis
# Sophia Nehs

# Don't forget to set working directory!

#######################################################################################################
# PACKAGES

# Load necessary packages
library('R.matlab') 
library('tidyverse') 
library('patchwork')

############################################################################
# POPULATION CONCENTRATION PROFILE VISUALIZATION

# Load the data saved from MATLAB
conc_a <- readMat("conc_a.mat")
pop_a <- readMat("pop_a.mat")
conc_b <- readMat("conc_b.mat")
pop_b <- readMat("pop_b.mat")
conc_c <- readMat("conc_c.mat")
pop_c <- readMat("pop_c.mat")
conc_cIIV <- readMat("conc_cIIV.mat")
pop_cIIV <- readMat("pop_cIIV.mat")

# Convert lists to data frames
conc_a <- as.data.frame(conc_a)
pop_a <- as.data.frame(pop_a)
conc_b <- as.data.frame(conc_b)
pop_b <- as.data.frame(pop_b)
conc_c <- as.data.frame(conc_c)
pop_c <- as.data.frame(pop_c)
conc_cIIV <- as.data.frame(conc_cIIV)
pop_cIIV <- as.data.frame(pop_cIIV)

# Name the columns
names(pop_a) <- c('Subject', 'Weight', 'Height', 'CrCL(mL/min)', 'BMI', 'CrCL', 'Dose', 'AUC','Ct')
CrCL_a <- t(pop_a$`CrCL(mL/min)`)
names(conc_a) <- c('Time', CrCL_a)
names(pop_b) <- c('Subject', 'Weight', 'Height', 'CrCL(mL/min)', 'BMI', 'CrCL', 'Dose', 'AUC','Ct')
Weight_b <- t(pop_b$Weight)
names(conc_b) <- c('Time', Weight_b)
names(pop_c) <- c('Subject', 'Weight', 'Height', 'CrCL(mL/min)', 'BMI', 'CrCL', 'Dose', 'AUC','Ct')
names(pop_cIIV) <- c('Subject', 'Weight', 'Height', 'CrCL(mL/min)', 'BMI', 'CrCL', 'Dose', 'AUC','Ct')
CrCL_c <- t(pop_c$CrCL)
names(conc_c) <- c('Time', CrCL_c)
names(conc_cIIV) <- c('Time', CrCL_c)

# Reformat data so all of the output is in one column
conc_a <- pivot_longer(conc_a, !Time, names_to = 'CrCL (mL/min)', values_to = 'Concentration')
conc_a$`CrCL (mL/min)` <- as.numeric(as.character(conc_a$`CrCL (mL/min)`))

conc_b <- pivot_longer(conc_b, !Time, names_to = 'Weight (kg)', values_to = 'Concentration')
conc_b$`Weight (kg)` <- as.numeric(as.character(conc_b$`Weight (kg)`))

conc_c <- pivot_longer(conc_c, !Time, names_to = 'CrCL (mL/min/1.73m^2)', values_to = 'Concentration')
conc_cIIV <- pivot_longer(conc_cIIV, !Time, names_to = 'CrCL (mL/min/1.73m^2)', values_to = 'Concentration')

conc_c$`CrCL (mL/min/1.73m^2)` <- as.numeric(as.character(conc_c$`CrCL (mL/min/1.73m^2)`))
conc_cIIV$`CrCL (mL/min/1.73m^2)` <- as.numeric(as.character(conc_cIIV$`CrCL (mL/min/1.73m^2)`))

# Set themeing
my_theme <- theme_classic() + 
  theme(text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(10,0,0,0)),
        axis.title.y = element_text(margin = margin(0,10,0,0)),
        axis.text.x = element_text(margin = margin(5,0,0,0)),
        axis.text.y = element_text(margin = margin(0,5,0,0)))


# Create plot
conc_a_plot <- ggplot(conc_a, aes(x = Time, y = Concentration, color = `CrCL (mL/min)`)) +
  geom_line(
    aes(group = `CrCL (mL/min)`),
    size = 1.0,
    alpha = 0.5
  ) +
  scale_x_continuous(name = "Time (hr)") +
  scale_y_continuous(name = "[D] (mg/L)") +
  labs(title = "Population A Concentration Profiles") +
  my_theme +
  theme(legend.title= element_text(size = 12) )


#conc_a_plot
# Save plot
ggsave(filename = 'PopA_conc.png', plot = conc_a_plot, width = 12, height = 8)


conc_b_plot <- ggplot(conc_b, aes(x = Time, y = Concentration, color = `Weight (kg)`)) +
  geom_line(
    aes(group = `Weight (kg)`),
    size = 1.0,
    alpha = 0.5
  ) +
  scale_x_continuous(name = "Time (hr)") +
  scale_y_continuous(name = "[D] (mg/L)") +
  labs(title = "Population B Concentration Profiles") +
  my_theme +
  theme(legend.title= element_text(size = 12) )


#conc_b_plot
# Save plot
ggsave(filename = 'PopB_conc.png', plot = conc_b_plot, width = 12, height = 8)


conc_c_plot <- ggplot(conc_c, aes(x = Time, y = Concentration, color = `CrCL (mL/min/1.73m^2)`)) +
  geom_line(
    aes(group = `CrCL (mL/min/1.73m^2)`),
    size = 1.0,
    alpha = 0.5
  ) +
  scale_x_continuous(name = "Time (hr)") +
  scale_y_continuous(name = "[D] (mg/L)") +
  labs(title = "Population C Concentration Profiles") +
  my_theme +
  theme(legend.title= element_text(size = 12) ) +
  # same Y axis for comparison
  coord_cartesian(ylim= c(0,115))


#conc_c_plot
#Save plot
ggsave(filename = 'PopC_conc.png', plot = conc_c_plot, width = 12, height = 8)

conc_cIIV_plot <- ggplot(conc_cIIV, aes(x = Time, y = Concentration, color = `CrCL (mL/min/1.73m^2)`)) +
  geom_line(
    aes(group = `CrCL (mL/min/1.73m^2)`),
    size = 1.0,
    alpha = 0.5
  ) +
  scale_x_continuous(name = "Time (hr)") +
  scale_y_continuous(name = "[D] (mg/L)") +
  labs(title = "Population C-IIV Concentration Profiles") +
  my_theme +
  theme(legend.title= element_text(size = 12) ) +
  # same Y axis for comparison
  coord_cartesian(ylim= c(0,115))


#conc_cIIV_plot
#Save plot
ggsave(filename = 'PopCIIV_conc.png', plot = conc_cIIV_plot, width = 12, height = 8)

############################################################################
# POPULATION DISTRIBUTIONS OF AUC AND CTROUGH

# Load the data saved from MATLAB
AUC <- readMat("AUCpops.mat")
Ct <- readMat("Ct10pops.mat")


# Convert lists to data frames
AUC <- as.data.frame(AUC)
Ct <- as.data.frame(Ct)


# Name the columns
names(AUC) <- c('Population A', 'Population B', 'Population C', 'Population C-IIV')
names(Ct) <- c('Population A', 'Population B', 'Population C', 'Population C-IIV')


# Reformat data so all of the output is in one column
AUC <- pivot_longer(AUC, cols = everything(), names_to = 'Population', values_to = 'AUC (mg*hr/L)')
Ct <- pivot_longer(Ct, cols = everything(), names_to = 'Population', values_to = 'Ctrough (mg/L)')

# Change parameter column to factor and specify order
AUC$Population <- factor(AUC$Population, levels = c('Population A', 'Population B', 'Population C', 'Population C-IIV'))
Ct$Population <- factor(Ct$Population, levels = c('Population A', 'Population B', 'Population C', 'Population C-IIV'))


# Create Plots


AUC_violin <- ggplot(data = AUC, aes(x = Population, y = `AUC (mg*hr/L)`, fill = Population)) +
  geom_violin(scale = "width") +
  geom_boxplot(fill = "white", width = 0.25) +
  labs(x = 'Population', y = 'AUC (mg*hr/L)', title = 'Population AUC Distributions') +
  my_theme +
  scale_fill_brewer(palette="Set1") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )

AUC_violin
# Save plot
ggsave(filename = 'AUC_vio.png', plot = AUC_violin, width = 12, height = 8)



Ct_violin <- ggplot(data = Ct, aes(x = Population, y = `Ctrough (mg/L)`, fill = Population)) +
  geom_violin(scale = "width") +
  geom_boxplot(fill = "white", width = 0.25) +
  labs(x = 'Population', y = 'Ctrought (mg/L)', title = 'Population Ctrough Distributions') +
  my_theme +
  scale_fill_brewer(palette="Set1") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )

Ct_violin
# Save plot
ggsave(filename = 'Ct_vio.png', plot = Ct_violin, width = 12, height = 8)


keymet_violins <- AUC_violin | Ct_violin +
  plot_layout(guides = 'collect')

keymet_violins
# Save plot
#ggsave(filename = 'keymet_vio.png', plot = keymet_violins, width = 12, height = 8)


###############################################################################
# Separating patients by dosage

# Read in MATLAB file
dose_metrics <- readMat("dose_metrics.mat")

# Convert to data frame
dose_metrics <- as.data.frame(dose_metrics)

# Name columns
names(dose_metrics) <- c('CrCL', 'Dose (mg)', 'AUC', 'Ctrough')

dose_metrics$`Dose (mg)`<- as.factor(dose_metrics$`Dose (mg)`)

doseAUC_violin <- ggplot(data = dose_metrics, aes(x = `Dose (mg)`, y = `AUC`, fill = `Dose (mg)`)) +
  geom_violin(scale = "width") +
  geom_boxplot(fill = "white", width = 0.25) +
  labs(x = 'Dosage (mg)', y = 'AUC (mg*hr/L)', title = 'Dose Based AUC Distributions') +
  my_theme +
  scale_fill_brewer(palette="Set1") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )

doseAUC_violin

doseCt_violin <- ggplot(data = dose_metrics, aes(x = `Dose (mg)`, y = `Ctrough`, fill = `Dose (mg)`)) +
  geom_violin(scale = "width") +
  geom_boxplot(fill = "white", width = 0.25) +
  labs(x = 'Dosage (mg)', y = 'Ctrough (mg/L)', title = 'Dose Based Ct Distributions') +
  my_theme +
  scale_fill_brewer(palette="Set1") +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 16)
  )

doseCt_violin
ggsave(filename = 'dose_ct.png', plot = doseCt_violin, width = 12, height = 8)


# Calculate Correlation between CrCL and Ctrough
# Pearson correlation between 2 variables
p <- cor.test(dose_metrics$CrCL, dose_metrics$Ctrough)
p_val <- p[["p.value"]]
R <- p[["estimate"]][["cor"]]

# Looking at correlation for just the "normal" kidney function
dose_met1500 <- dose_metrics %>% filter(`Dose (mg)` == "1500")
p1500 <- cor.test(dose_met1500$CrCL, dose_met1500$Ctrough)
p_val1500 <- p1500[["p.value"]]
R1500 <- p1500[["estimate"]][["cor"]]

CrCLvCt_plot <- 
  ggplot(data=dose_metrics, aes(
  x = CrCL,
  y = Ctrough,
  color = `Dose (mg)`
)) +
  geom_point(size = 3.5, alpha = 0.7) +
  annotate(geom="text", x=160, y=80, label = "R = -0.615", size=8) +
  #geom_smooth() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "CrCL (mL/min/1.73m^2)", y = "Ctrough (mg/L)", title = "CrCL and Ctrough Relationship") +
  #geom_smooth(color='black') +
  my_theme

CrCLvCt_plot
ggsave(filename = 'dose_ct_scatter.png', plot = CrCLvCt_plot, width = 12, height = 8)


Ctrel_plot <- doseCt_violin | CrCLvCt_plot +
  plot_layout(guides = 'collect')
  
ggsave(filename = 'ct_rel.png', plot = Ctrel_plot, width = 12, height = 8)



CrCLvCt15_plot <- dose_metrics %>% filter(`Dose (mg)` == "1500") %>%
  ggplot(aes(
    x = CrCL,
    y = Ctrough,
    color = `Dose (mg)`
  )) +
  geom_point(size = 3.5, alpha = 0.7, color="#66A61E") +
  annotate(geom="text", x=160, y=80, label = "R = -0.615") +
  #geom_smooth() +
  #scale_color_brewer(palette = "Set1") +
  labs(x = "CrCL (mL/min/1.73m^2)", y = "Ctrough (mg/L)", title = "CrCL and Ctrough Relationship") +
  #geom_smooth(color='black') +
  my_theme

CrCLvCt15_plot
ggsave(filename = 'dose_ct15_scatter.png', plot = CrCLvCt15_plot, width = 12, height = 8)


