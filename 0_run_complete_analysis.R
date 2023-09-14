
# Run the complete analysis for AMR demographics paper

# some global variables to define
sensitivity_0_run <- "OFF" # "ON" or "OFF". If on, runs the analysis including aged 0 
sensitivity_prior_run <- "OFF" # "ON" or "OFF". If on, runs the fits with a regularising prior

##### DATA PREP ####
# clean the data
source("1_1_data_cleaning.R")
# generate initial plots of the data (incl. Fig's X, X)
source("1_2_data_plotting.R")

##### INCIDENCE CALCS ####
# calculates the incidence rates (incl. Figs X, X, X)
source("2_1_incidence.R")

##### MODELLING ####
# creates the model fits. NOTE: manually specify which models to fit at the top of the script
file.edit("3_1_Model_running_main.R")
# generate Model output plots (incl. Figs 3 & 4)
source("3_2_Model_output.R")
# generate comparison plots across models
source("3_3_model_comparisons.R")
# generate further comparison plots across groups
source("3_4_model_raw_comps.R")

######## SIDE ANALYSES #######
# investigate missing data (incl. Figs X)
source("4_1_missing_data.R")
# look at variance components model to understand variation
source("4_2_variance_components_model.R")
# run sensitivity analyses on modelling results
source("4_3_model_out_sensitivities.R")
# subregional analyses
source("4_4_spatial_variation.R")


