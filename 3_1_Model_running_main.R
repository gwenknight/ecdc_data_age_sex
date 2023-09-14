# main script from which to run models. 

# USER input!
pathogen_to_run <-  "staaur"  #"staaur" # pathogen to run (needs to match pathogen" column in data)
resistance_to_run <- "mrsa_R" #"mrsa_R" # resistance to run (needs to match column in data)
models_to_run <-c("D")   #c("A" ,"B", "C", "D") # choose which models you want to run. (Any combo of A, B, C and D). DEFAULT = "D"
how_many_cores_free <- 3 # how many cores to NOT use. I.e. cores used will be total cores available - how_many_cores_free
it_no <- 3000 # DEFAULR 3000 unless it hasn't converged.
sensitivity_0_run_now <- "OFF" # DEFAULT "OFF". Runs the sensitivity analysis including those aged 0
sensitivity_prior_run_now <- "OFF" # DEFAULT "OFF". runs the sensitivity analysis with normal prior (0,1) on the fixed parameters

# checks if already defined globally. if not takes the local value
if(exists("sensitivity_0_run")){
  sensitivity_0 <- sensitivity_0_run
} else {
  sensitivity_0 <- sensitivity_0_run_now
}
# checks if already defined globally. if not takes the local value
if(exists("sensitivity_prior_run")){
  sensitivity_prior <- sensitivity_prior_run
} else {
  sensitivity_prior <- sensitivity_prior_run_now
}


####### SETUP #####
set.seed(500) 
# libraries etc. 
source("3_a_load_needed.R")
#source the functions
source("3_b_modelling_functions.R") # 4 models with default priors
# If prior sensitivity to be run the use other models with normal priors
if(sensitivity_prior == "ON"){
  source("3_c_modelling_functions2.R")
}

# Read in the required data. 
if(sensitivity_0 == "OFF"){
  data<- as.data.table(read_csv("data/data_cleaned.csv")[,-1])
} else if(sensitivity_0 == "ON") {
  data<- as.data.table(read_csv("data/data_cleaned_with0s.csv")[,-1])
}

data_subset <- data[pathogen == pathogen_to_run & name == resistance_to_run]
rm(data)
# format and standardise varaibles
data_subset <- format_standard(data_subset, resistance_to_run)

print(paste0("There are ",nrow(data_subset), " rows in the binomial data"))

# this will run all the specified models and save them. 
run_models(data_binomial = data_subset, models_to_run = models_to_run, 
           it_no = it_no, sensitivity_0 = sensitivity_0)

