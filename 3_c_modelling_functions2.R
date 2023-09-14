# functions for running the bug-drug combo models


format_standard <- function(data_subset, resistance_to_run){
  
  # make all the variables on the same scale
  data_subset[, age_s := age/100]
  data_subset[, age_squared_s := age_s*age_s]
  
  data_subset[year == "2015", year_s := 0]
  data_subset[year == "2016", year_s := 0.25]
  data_subset[year == "2017", year_s := 0.5]
  data_subset[year == "2018", year_s := 0.75]
  data_subset[year == "2019", year_s := 1]
  
  data_subset[, tot_samples := sus+res]
  
  return(data_subset)
}

missingness_function <- function(data_subset){
  # calculate missingness 
  age_missing <- sum(is.na(data_subset$age))
  gender_missing <- sum(!data_subset$gender %in% c("m", "f"))
  missing_both <- sum(data_subset[, is.na(age) & !gender %in% c("m", "f")])
  total_rows <- nrow(data_subset)
  age_percent_missing<- age_missing/total_rows*100
  gender_percent_missing <- gender_missing/total_rows*100
  
  missingness <- c("age_missing" = age_missing, 
                   "gender_missing" = gender_missing, 
                   "both_missing" = missing_both,
                   "total_samples" = total_rows, 
                   "age_percent_missing" = age_percent_missing,
                   "gender_percent" = gender_percent_missing)
  
  saveRDS(missingness,file = paste0("output/missingness_",pathogen_to_run, ".rds"))
}


run_models <- function(data_binomial, models_to_run, it_no, sensitivity_0 = "OFF"){
  
  if("D" %in% models_to_run){
    Model_D <- brm(data = data_binomial, 
                   family = binomial,
                   res | trials(tot_samples) ~ 1 +
                     age_s + 
                     age_squared_s + 
                     gender + 
                     gender*age_s + 
                     year_s + 
                     (1 + age_s | country) + 
                     (1 | laboratorycode),
                   set_prior("normal(0, 1)", class = "b"),
                   file = paste0("fits/D_prior_",pathogen_to_run,"_",resistance_to_run,"_", Sys.Date()),
                   file_refit = getOption("brms.file_refit", "always"),
                   iter = it_no, save_pars = save_pars(all = TRUE))
  }
}