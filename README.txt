README.txt ---------------------------------------------------

GROUP 3 - Levetiracetam (CC, VD, CG, SN, SB)


FILES INCLUDED:
MATLAB - 
LEV_eqns.m
LEV_sim.m
LEV_sim_peak.m
LEV_popdriver.m
LEV_sensitivitydriver.m
LEV_PDmodel.m
missed_dose_without_pop.m
missed_dose_pop_5patients.m
missed_dose_diff_set_pop.m

R - 
LEV_pop_Vis.R
LEV_SensitivityPlots.R
LEV_PDmodel_vis.R
Fig_16.R
Fig_17.R
Fig_18.R
Fig_19_20.R
Fig_21.R

Folder - (separate folder for interactive visualization)
lev_misseddose_app


INSTRUCTIONS:

1. Run the following MATLAB codes to save all mat files necessary for visualization (note any figures produced by these files were preliminary, not used in the report)
	1. LEV_popdriver.m - run as is
	2. LEV_PDmodel.m - run as is
	3. LEV_sensitivitydriver.m - run all 5 cases IN ORDER 1-5 (line 5)
	4. missed_dose_without_pop.m - run all 4 cases (line 54), order does not matter
	5. missed_dose_pop_5patients.m - run as is to save all necessary files
	6. missed_dose_diff_set_pop.m - run all 5 cases (line 46), order does not matter

2. Run the following R codes to generate visualizations (note there are some supplementary figures produced that are not included in the report to avoid redundancies)
	1. LEV_pop_Vis.R - run as is (the code is a little slow, but it works!)
	2. LEV_PDmodel_vis.R - run as is (note, you may need to install ggfortify and survminer packages to run)
	3. LEV_sensitivityPlots.R - run as is
	4. Missed dose visualization (ggpubr package is required) - run as is
		1. Fig_16.R
		2. Fig_17.R
		3. Fig_18.R
		4. Fig_19_20.R
		5. Fig_21.R


3. Run the shiny app - navigate to 'lev_misseddose_app' (set new working directory) and Run App. Note, all data used for the interactive visualization is contained in the data folder within the lev_misseddose_app folder.



