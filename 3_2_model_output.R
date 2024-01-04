# look at model output!

# What bacteria-antibiotic combinations are considered? 
path_drug_combos <- data.frame(read.csv2("path_drug_combos.csv", sep = ",", header = F))

# Which models to consider? 
only_model_D <- T

# Load libraries
source("3_a_load_needed.R")
files_to_search <- list.files("fits")

# finds the latest file with the correct specs
load_model <- function(pathogen_to_run, resistance_to_run, model_c, files_to_search){
  
  file_no <- tail(grep(pattern = paste0(model_c, "_", pathogen_to_run, "_", resistance_to_run), 
                       x = files_to_search),1)
  loaded_model <- brm(file = paste0("fits/",files_to_search[file_no]))
  
  return(loaded_model)
  
}

# storage for parameter values
fixed_effs_storage <- list()

# Read in the required data. 
data<- as.data.table(read_csv("data/data_cleaned.csv")[,-1])
country_reference <- as.data.table(read.csv("countries_reference.csv"))
colnames(country_reference) <- c("country", "country_long","anon_country")
data[country_reference, on = c("country"), country_long := i.country_long]
translate_bug_drugs <-as.data.table(read.csv("data/translate_drug_bugs.csv", header = F))

#### side step - calculate number of bug_drugs by country for supplementary 
data[, total := sus + res]
samples_by_country <- data[, sum(total, na.rm =T), by = c("country", "pathogen", "name")]
country_out <- dcast(samples_by_country, pathogen + name ~ country, value.var = "V1")
country_out[is.na(country_out)] <- 0

### Create storage folder 
dir.create(file.path("output"), showWarnings = FALSE)
write.csv(country_out, file = "output/Country_samples.csv")

# storage for overall parameter table
all_params_out <- data.table()

# and loop through all bacteria-antibiotic combinations
for(i in 1:nrow(path_drug_combos)){
  
  # Combinations
  pathogen_to_run <-  path_drug_combos[i,1]  # choose from: c("acispp", "encfae", "encfai", "esccol", "klepne", "pseaer", "staaur", "strpne")
  resistance_to_run <- path_drug_combos[i,2]
  
  bug_drug <- translate_bug_drugs[V2 == pathogen_to_run & V3 == resistance_to_run][,c("V1", "V6")]
  data_subset <- data[pathogen == pathogen_to_run & name == resistance_to_run]
  
  ############################################## Load models
  if(only_model_D == F){
    Model_A <- load_model(pathogen_to_run = pathogen_to_run, 
                          resistance_to_run = resistance_to_run, 
                          model_c = "A",
                          files_to_search = files_to_search)
    
    Model_B <- load_model(pathogen_to_run = pathogen_to_run, 
                          resistance_to_run = resistance_to_run, 
                          model_c = "B",
                          files_to_search = files_to_search)
    
    Model_C <- load_model(pathogen_to_run = pathogen_to_run, 
                          resistance_to_run = resistance_to_run, 
                          model_c = "C",
                          files_to_search = files_to_search)
  }
  Model_D <- load_model(pathogen_to_run = pathogen_to_run, 
                        resistance_to_run = resistance_to_run, 
                        model_c = "D",
                        files_to_search = files_to_search)
  
  ############################################## Make a summary table of the model output
  if(only_model_D == F){
    
    models_overview <- data.frame(matrix(NA, ncol = 8, nrow = 6))
    rownames(models_overview) <- c("intercept", "year", "age", "age^2", "gender(m)", 
                                   "Age:gender(m)")
    colnames(models_overview) <- c("A estimate", "A se", "B estimate", "B se", 
                                   "C estimate", "C se", "D estimate", "D se")
    models_overview[,"A estimate"] <- c(fixef(Model_A)[c("Intercept", "year_s"),"Estimate"], 
                                        rep(NA, 4))
    models_overview[,"A se"] <- c(fixef(Model_A)[c("Intercept", "year_s"),"Est.Error"], 
                                  rep(NA, 4))
    
    models_overview[,"B estimate"] <- c(fixef(Model_B)[c("Intercept", "year_s",
                                                         "age_s", "age_squared_s"),"Estimate"], 
                                        rep(NA, 2))
    models_overview[,"B se"] <- c(fixef(Model_B)[c("Intercept", "year_s",
                                                   "age_s", "age_squared_s"),"Est.Error"], 
                                  rep(NA, 2))
    
    models_overview[,"C estimate"] <- c(fixef(Model_C)[c("Intercept", "year_s"),"Estimate"], 
                                        rep(NA, 2), fixef(Model_C)[c("genderm"),"Estimate"], 
                                        NA)
    models_overview[,"C se"] <- c(fixef(Model_C)[c("Intercept", "year_s"),"Est.Error"], 
                                  rep(NA, 2), fixef(Model_C)[c("genderm"),"Est.Error"], 
                                  NA)
    
    models_overview[,"D estimate"] <- c(fixef(Model_D)[c("Intercept", "year_s", 
                                                         "age_s", "age_squared_s", 
                                                         "genderm", "age_s:genderm"),"Estimate"])
    models_overview[,"D se"] <- c(fixef(Model_D)[c("Intercept", "year_s", 
                                                   "age_s", "age_squared_s", 
                                                   "genderm", "age_s:genderm"),"Est.Error"])
  }
  
  if(only_model_D == T){
    
    models_overview <- data.frame(matrix(NA, ncol = 5, nrow = 6))
    rownames(models_overview) <- c("intercept", "year", "age", "age^2", "gender(m)", 
                                   "Age:gender(m)")
    colnames(models_overview) <- c("Estimate", "Error", "Lower 95%", "Upper 95%", "95% Credibility Interval")
    
    models_overview[,"Estimate"] <- c(fixef(Model_D)[c("Intercept", "year_s", 
                                                       "age_s", "age_squared_s", 
                                                       "genderm", "age_s:genderm"),"Estimate"])
    models_overview[,"Error"] <- c(fixef(Model_D)[c("Intercept", "year_s", 
                                                    "age_s", "age_squared_s", 
                                                    "genderm", "age_s:genderm"),"Est.Error"])
    
    models_overview <- as.data.table(models_overview)
    models_overview[, `Upper 95%` := Estimate + 1.96*Error]
    models_overview[, `Lower 95%` := Estimate - 1.96*Error]
    models_overview <- as.data.frame(models_overview)
    rownames(models_overview) <- c("intercept", "year", "age", "age^2", "gender(m)", 
                                   "Age:gender(m)")
    models_overview[,1] <- format(round(models_overview[,1],2), nsmall =2)
    models_overview[,2] <- format(round(models_overview[,2],2), nsmall =2)
    models_overview[,3] <- format(round(models_overview[,3],2), nsmall =2)
    models_overview[,4] <- format(round(models_overview[,4],2), nsmall =2)
    
    models_overview[,5] <- paste0("(",models_overview[,3],", ",models_overview[,4],")")
    
    models_table <- as.matrix(models_overview[,c(1,2,5)])
    colnames(models_table) <- c("Estimate","Standard Error","95% Credibility Interval")
    models_table <-tableGrob(models_table)
    title <- textGrob("A",gp=gpar(fontsize=12))
    padding <- unit(8,"mm")
    
    models_table_all <- gtable_add_rows(
      models_table,
      heights = grobHeight(title) + padding,
      pos = 0 )
    models_table_all <- gtable_add_grob(
      models_table_all,
      title,
      1, 1, 1)
  }
  
  ############################################## extract the fixed effects to later compare
  temp <- as.data.frame(fixef(Model_D, summary = F))
  temp$bug_drug <- paste0(pathogen_to_run, "_", resistance_to_run)
  fixed_effs_storage[[i]] <- temp
  
  
  ############################################## decomposition plot. for Model D
  model_to_test <- Model_D
  
  # set up storage
  decomposition <- data.table(age = seq(0,100),
                              age_s= seq(0,100)/100)
  decomposition[, age_squared_s := age_s*age_s]
  
  # extract coefficient values from model
  age_coef <- mean(as.data.frame(model_to_test$fit)$b_age_s)
  age_coef_2 <- mean(as.data.frame(model_to_test$fit)$b_age_squared_s)
  intercept <- mean(as.data.frame(model_to_test$fit)$b_Intercept)
  gender <- mean(as.data.frame(model_to_test$fit)$b_genderm)
  age_gender <- mean(as.data.frame(model_to_test$fit)$`b_age_s:genderm`)
  
  # calculate the decomposed aspects
  decomposition[, age_line := intercept + age_s*age_coef]
  decomposition[, age_squared_line := intercept + age_squared_s*age_coef_2]
  decomposition[, male_age_line := intercept + gender*1  + 1*age_s * age_gender]
  decomposition[, combined_male := intercept + gender*1  + 1*age_s * age_gender + age_squared_s*age_coef_2 + age_s*age_coef ]
  decomposition[, combined_female := intercept + gender*0  + 0*(age_s * age_gender) + age_squared_s*age_coef_2 + age_s*age_coef ]
  
  
  decomposition_m <- melt.data.table(decomposition, id.vars= c("age"),
                                     measure.vars = c("age_line",
                                                      "age_squared_line",
                                                      "male_age_line",
                                                      "combined_male",
                                                      "combined_female"
                                     ))
  
  DECOMPOSITION <- ggplot(decomposition_m, aes(x = age, y = inv.logit(value), colour = variable, linetype=variable)) +
    geom_line() +
    scale_colour_manual(values = c("gold2", "darkorange", "darkorange4", "black", "darkorange")) +
    scale_linetype_manual(values=c(1,1,1,1,2)) +
    labs( x = "Age", y = "predicted value", title = "E")
  
  # make plot folder
  dir.create(file.path("output_figures", "by_bug_drug"), showWarnings = FALSE) # don't want error if already exists
  tiff(paste0("output_figures/by_bug_drug/", pathogen_to_run, "_", resistance_to_run, "_decomposition.tiff"),
       width = 3250, height = 2000, res = 300)
  print(DECOMPOSITION)
  dev.off()
  
  ############################################## extract the random effects for the country age slope
  raneffs_D <- ranef(Model_D, summary = F)
  raneffs_D_country <- as.data.table(raneffs_D$country)
  colnames(raneffs_D_country) <- c("Sample", "country", "variable", "estimate")
  raneffs_D_lab <- as.data.table(raneffs_D$laboratorycode)
  colnames(raneffs_D_lab) <- c("Sample", "lab", "variable", "lab_estimate")
  
  raneffs_D_country_age <- raneffs_D_country[variable == "age_s"]
  raneffs_D_country_intercept <- raneffs_D_country[variable == "Intercept"]
  raneffs_D_country_age[raneffs_D_country_intercept, on = c("Sample", "country"), country_intercept := i.estimate]
  colnames(raneffs_D_country_age) <- c("Sample", "country", "bla", "age_country_slope", "country_intercept")
  raneffs_D_country_age[,bla:=NULL]
  
  ############################################## extract the fixed effects
  fixeds_D <- as.data.table(fixef(Model_D, summary = F))
  fixeds_D$Sample <- 1:nrow(fixeds_D)
  # add the fixed effs to the country level ones. I.e. same table
  calculating_effects <- raneffs_D_country_age[fixeds_D, on = c("Sample"), ]
  
  ############################################## extract subsample for future burden estimates and to compare estimates at age 1 and age 100
  # For future burden
  n <- 100
  calculating_effects_all <- do.call("rbind", replicate(n, calculating_effects, simplify = FALSE))
  calculating_effects_all[, age := rep(1:100, each = nrow(calculating_effects))]
  
  # to compare estimate at age 1 and age 100 
  calculating_effects$age <- 1
  temp <- calculating_effects
  calculating_effects$age <- 100
  calculating_effects <- rbind(temp, calculating_effects)
  
  #### calculate the values
  ## For future burden
  calculating_effects_all[, female :=  Intercept +  country_intercept + age_squared_s*(age/100)*(age/100) + age_s*(age/100) + age_country_slope * (age/100)+ 1*year_s]
  # male int.   # male age slope
  calculating_effects_all[, male := Intercept +  country_intercept + genderm*1  + (age/100) * `age_s:genderm`  + age_squared_s*(age/100)*(age/100) + age_s*(age/100) + age_country_slope * (age/100)+ 1*year_s]
  
  ## For difference age 1 vs 100 
  # age squared term      #age terms.       # country level age effect #2019 effect
  calculating_effects[, female :=  Intercept  +  country_intercept + age_squared_s*(age/100)*(age/100) + age_s*(age/100) + age_country_slope * (age/100) + 1*year_s]
  # male int.   # male age slope
  calculating_effects[, male := Intercept +  country_intercept + genderm*1  + (age/100) * `age_s:genderm`  + age_squared_s*(age/100)*(age/100) + age_s*(age/100) + age_country_slope * (age/100)+ 1*year_s]
  
  
  #### side-point: see which ages have the max age
  which.max(inv.logit(unlist(calculating_effects_all[,mean(female), by = "age"]$V1)))
  which.min(inv.logit(unlist(calculating_effects_all[,mean(female), by = "age"]$V1)))
  which.max(inv.logit(unlist(calculating_effects_all[,mean(male), by = "age"]$V1)))
  which.min(inv.logit(unlist(calculating_effects_all[,mean(male), by = "age"]$V1)))
  
  ## Summmarise and grab IQR + median
  all_summary <- calculating_effects_all[, quantile(female, probs = c(0.025, 0.5, 0.975)), by = c("country", "age")]
  colnames(all_summary)[3] <- "female"
  all_summary$male <- calculating_effects_all[, quantile(male, probs = c(0.025, 0.5, 0.975)),by = c("country", "age")]$V1
  all_summary[,female := inv.logit(female)]
  all_summary[,male := inv.logit(male)]
  all_summary$type <- rep(c("lower", "median", "upper"), length(unique(all_summary$country))*100)
  all_summary_m <- melt.data.table(all_summary, id.vars = c("country", "age", "type"))
  all_summary_mc <- dcast.data.table(all_summary_m, country + age + variable ~ type, value.var = "value")
  
  #reformat so can compare 100 vs 1 year olds
  calculating_effects_100 <- calculating_effects[ age == 100]
  calculating_effects_1 <- calculating_effects[ age == 1]
  calculating_effects_1[calculating_effects_100, on = c("Sample", "country"), female_100 := i.female]
  calculating_effects_1[calculating_effects_100, on = c("Sample", "country"), male_100 := i.male]
  
  calculating_effects_1[, female := inv.logit(female)]
  calculating_effects_1[, male := inv.logit(male)]
  calculating_effects_1[, female_100 := inv.logit(female_100)]
  calculating_effects_1[, male_100 := inv.logit(male_100)]
  
  calculating_effects_1[, female_difference := female_100 - female]
  calculating_effects_1[, male_difference := male_100 - male]
  
  summary_calculations <- calculating_effects_1[, quantile(female_difference, probs = c(0.025, 0.5, 0.975)), by = c("country")]
  colnames(summary_calculations)[2] <- "female_difference"
  summary_calculations$male_difference <- calculating_effects_1[, quantile(male_difference, probs = c(0.025, 0.5, 0.975)), by = c("country")]$V1
  summary_calculations$type <- rep(c("lower", "median", "upper"), length(unique(summary_calculations$country)))
  female <- dcast.data.table(summary_calculations, country ~ type, value.var = "female_difference")
  male <- dcast.data.table(summary_calculations, country ~ type, value.var = "male_difference")
  female$gender <- "female"
  male$gender <- "male"
  
  combo <- rbind(female, male)
  combo <- combo[order(gender, median)]
  
  combo$country <- factor(combo$country, levels = c(combo$country)[1:length(unique(combo$country))])
  combo[country_reference, on = c("country"), country_long := i.country_long]
  combo$country_long <- factor(combo$country_long, levels = c(combo$country_long)[1:length(unique(combo$country_long))])
  combo <- left_join(combo, country_reference)
  combo <- combo[order(gender, median)]
  combo$anon_country <- factor(combo$anon_country, levels = c(combo$anon_country)[1:length(unique(combo$anon_country))])
  
  ACROSS_COUNTRY <- ggplot(combo, aes(x = anon_country, y = median, colour = gender )) +
    # geom_point() +
    geom_hline( yintercept= 0) +
    geom_pointrange(aes(ymin=lower, ymax = upper), position = position_dodge(width = 0.4)) +
    labs( y = "Change in propoortion resistant between ages 1 and 100",
          x = "Country",
          colour = "Gender",
          title = paste0("Bacteria: ", pathogen_to_run, ", antibiotic: ", resistance_to_run)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_colour_manual(values = colours_to_use[1:2])
  
  
  tiff(paste0("output_figures/by_bug_drug/", pathogen_to_run, "_", resistance_to_run, "_countries.tiff"),
       width = 3250, height = 2000, res = 300)
  print(ACROSS_COUNTRY)
  dev.off()
  
  ############################################## save extra info for storing
  combo$bug <- pathogen_to_run
  combo$resistance <- resistance_to_run
  calculating_effects_1$bug <- pathogen_to_run
  calculating_effects_1$resistance <- resistance_to_run
  
  saveRDS(combo, file = paste0("output/comparison_",pathogen_to_run, "_", resistance_to_run, ".rds" ), )
  saveRDS(calculating_effects_1, file = paste0("output/raw_comp_",pathogen_to_run, "_", resistance_to_run, ".rds"))
  saveRDS(all_summary_mc, file = paste0("output/demo_summary_",pathogen_to_run, "_", resistance_to_run, ".rds"))
  
  ############################################## Understand impact of age and gender across models
  if(only_model_D == F){
    raneffs_A <- ranef(Model_A)
    raneffs_A_country <- as.data.frame(raneffs_A$country)
    raneffs_A_country$country <- rownames(raneffs_A_country)
    raneffs_A_country %>% mutate(raneffs_A_country, country = reorder(country, Estimate.Intercept))
    raneffs_A_country$Model <- "A"
    
    raneffs_B <- ranef(Model_B)
    raneffs_B_country <- as.data.frame(raneffs_B$country)
    raneffs_B_country$country <- rownames(raneffs_B_country)
    raneffs_B_country %>% mutate(raneffs_B_country, country = reorder(country, Estimate.Intercept))
    raneffs_B_country$Model <- "B"
    
    raneffs_C <- ranef(Model_C)
    raneffs_C_country <- as.data.frame(raneffs_C$country)
    raneffs_C_country$country <- rownames(raneffs_C_country)
    raneffs_C_country %>% mutate(raneffs_C_country, country = reorder(country, Estimate.Intercept))
    raneffs_C_country$Model <- "C"
    
    raneffs_D <- ranef(Model_D)
    raneffs_D_country <- as.data.frame(raneffs_D$country)
    raneffs_D_country$country <- rownames(raneffs_D_country)
    raneffs_D_country %>% mutate(raneffs_D_country, country = reorder(country, Estimate.Intercept))
    raneffs_D_country$Model <- "D"
    
    raneffs_combo_country <- rbind(raneffs_A_country,
                                   raneffs_B_country[,c("Estimate.Intercept",
                                                        "Est.Error.Intercept",
                                                        "Q2.5.Intercept",
                                                        "Q97.5.Intercept",
                                                        "country",
                                                        "Model")],
                                   raneffs_C_country[,c("Estimate.Intercept",
                                                        "Est.Error.Intercept",
                                                        "Q2.5.Intercept",
                                                        "Q97.5.Intercept",
                                                        "country",
                                                        "Model")],
                                   raneffs_D_country[,c("Estimate.Intercept",
                                                        "Est.Error.Intercept",
                                                        "Q2.5.Intercept",
                                                        "Q97.5.Intercept",
                                                        "country",
                                                        "Model")])
    
    raneffs_combo_country$country <- factor(raneffs_combo_country$country)
    to_plot <- mutate(raneffs_combo_country,country = reorder(country, Estimate.Intercept))
    
    COUNTRY_INTERCEPTS <- ggplot(to_plot, aes(x = country, y = Estimate.Intercept, colour = Model)) +
      geom_pointrange( aes(ymin = Q2.5.Intercept, ymax = Q97.5.Intercept
      ),position = position_dodge(width = 0.4), size = 0.3) +
      scale_color_brewer(type = "qual", palette =2) +
      # scale_color_manual(values = colours_to_use[3:6]) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      labs( x = "Country", title = "A",
            y= "Country random effects")
    
    
    raneffs_A_lab <- as.data.frame(raneffs_A$laboratorycode)
    raneffs_A_lab$lab <- rownames(raneffs_A_lab)
    raneffs_A_lab %>% mutate(raneffs_A_lab, lab = reorder(lab, Estimate.Intercept))
    raneffs_A_lab$Model <- "A"
    raneffs_B_lab <- as.data.frame(raneffs_B$laboratorycode)
    raneffs_B_lab$lab <- rownames(raneffs_B_lab)
    raneffs_B_lab %>% mutate(raneffs_B_lab, lab = reorder(lab, Estimate.Intercept))
    raneffs_B_lab$Model <- "B"
    raneffs_C_lab <- as.data.frame(raneffs_C$laboratorycode)
    raneffs_C_lab$lab <- rownames(raneffs_C_lab)
    raneffs_C_lab %>% mutate(raneffs_C_lab, lab = reorder(lab, Estimate.Intercept))
    raneffs_C_lab$Model <- "C"
    raneffs_D_lab <- as.data.frame(raneffs_D$laboratorycode)
    raneffs_D_lab$lab <- rownames(raneffs_D_lab)
    raneffs_D_lab %>% mutate(raneffs_D_lab, lab = reorder(lab, Estimate.Intercept))
    raneffs_D_lab$Model <- "D"
    
    raneffs_combo_lab <- rbind(raneffs_A_lab,
                               raneffs_B_lab[,c("Estimate.Intercept",
                                                "Est.Error.Intercept",
                                                "Q2.5.Intercept",
                                                "Q97.5.Intercept",
                                                "lab",
                                                "Model")],
                               raneffs_C_lab[,c("Estimate.Intercept",
                                                "Est.Error.Intercept",
                                                "Q2.5.Intercept",
                                                "Q97.5.Intercept",
                                                "lab",
                                                "Model")],
                               raneffs_D_lab[,c("Estimate.Intercept",
                                                "Est.Error.Intercept",
                                                "Q2.5.Intercept",
                                                "Q97.5.Intercept",
                                                "lab",
                                                "Model")])
    
    raneffs_combo_lab$lab <- factor(raneffs_combo_lab$lab)
    to_plot2 <- mutate(raneffs_combo_lab,lab = reorder(lab, Estimate.Intercept))
    
    
    LAB_INTERCEPTS <- ggplot(to_plot2, aes(x = lab, y = Estimate.Intercept, colour = Model)) +
      geom_pointrange( aes(ymin = Q2.5.Intercept, ymax = Q97.5.Intercept
      ),position = position_dodge(width = 0.4)) +
      facet_grid(Model~.) +
      scale_color_manual(values = colours_to_use[3:6]) +
      labs( x = "Laboratory", title = "C",
            y= "Estimated difference in random intercept from all labs mean (with 95% uncertainty)") +
      theme(axis.text.x = element_blank())
    
    raneffs_combo_country_age <- rbind(
      raneffs_B_country[,c("Estimate.age_s",
                           "Est.Error.Intercept",
                           "Q2.5.age_s",
                           "Q97.5.age_s",
                           "country",
                           "Model")],
      
      raneffs_D_country[,c("Estimate.age_s",
                           "Est.Error.Intercept",
                           "Q2.5.age_s",
                           "Q97.5.age_s",
                           "country",
                           "Model")])
    
    raneffs_combo_country_age$country <- factor(raneffs_combo_country_age$country)
    to_plot3 <- mutate(raneffs_combo_country_age,country = reorder(country, Estimate.age_s))
    
    COUNTRY_AGE_INTERCEPTS <- ggplot(to_plot3, aes(x = country, y = Estimate.age_s, colour = Model)) +
      geom_pointrange( aes(ymin = Q2.5.age_s, ymax = Q97.5.age_s
      ),position = position_dodge(width = 0.4), size = 0.4) +
      # facet_grid(Model~.) +
      scale_color_manual(values = colours_to_use[c(4,6)]) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      labs( x = "Country", y= "Country - age random effects",
            title = "B")
    
  }
  
  raneffs_D <- ranef(Model_D)
  raneffs_D_country <- as.data.frame(raneffs_D$country)
  raneffs_D_country$country <- rownames(raneffs_D_country)
  raneffs_D_country %>% mutate(raneffs_D_country, country = reorder(country, Estimate.Intercept))
  raneffs_D_country <- as.data.table(raneffs_D_country)
  
  raneffs_D_lab <- as.data.frame(raneffs_D$laboratorycode)
  raneffs_D_lab$lab <- rownames(raneffs_D_lab)
  #raneffs_D_lab %>% mutate(raneffs_D_lab, lab = reorder(lab, Estimate.Intercept)) # REMOVE? 
  raneffs_D_lab <- as.data.table(raneffs_D_lab)
  
  raneffs_D_lab[, country := substr(lab, 1, 2)]
  
  ############################################## Find country differences
  # country with highest, lowest and median slope
  # middle country to choose (can't use just median as this will give value in between if even numbers)
  selector_mid <- nrow(raneffs_D_country)/2 #take the lower of the two values
  
  max_country <- as.character(raneffs_D_country[which.max(unlist(raneffs_D_country[, "Estimate.age_s"])),]$country)
  min_country <- as.character(raneffs_D_country[which.min(unlist(raneffs_D_country[, "Estimate.age_s"])),]$country)
  mid_country <- as.character(raneffs_D_country[which(unlist(raneffs_D_country[, "Estimate.age_s"]) ==
                                                        sort(unlist(raneffs_D_country[, "Estimate.age_s"]), decreasing=T)[selector_mid]),]$country)
  
  # ages
  min_age <- min(data_subset[country %in% c(max_country, min_country, mid_country),]$age)
  max_age <- max(data_subset[country %in% c(max_country, min_country, mid_country),]$age)
  
  temp = which.max(unlist(raneffs_D_lab[ like(lab, tolower(max_country)), "Estimate.Intercept"]))
  max_lab_max_country = as.character(unlist(raneffs_D_lab[  like(lab, tolower(max_country)), "lab"][temp]))
  temp = which.min(unlist(raneffs_D_lab[ like(lab, tolower(max_country)), "Estimate.Intercept"]))
  min_lab_max_country = as.character(unlist(raneffs_D_lab[ like(lab, tolower(max_country)), "lab"][temp]))
  
  temp = which.max(unlist(raneffs_D_lab[like(lab, tolower(mid_country)), "Estimate.Intercept"]))
  max_lab_mid_country = as.character(unlist(raneffs_D_lab[like(lab, tolower(mid_country)), "lab"][temp]))
  temp = which.min(unlist(raneffs_D_lab[like(lab, tolower(mid_country)), "Estimate.Intercept"]))
  min_lab_mid_country = as.character(unlist(raneffs_D_lab[like(lab, tolower(mid_country)), "lab"][temp]))
  
  temp = which.max(unlist(raneffs_D_lab[ like(lab, tolower(min_country)), "Estimate.Intercept"]))
  max_lab_min_country = as.character(unlist(raneffs_D_lab[ like(lab, tolower(min_country)), "lab"][temp]))
  temp = which.min(unlist(raneffs_D_lab[  like(lab, tolower(min_country)), "Estimate.Intercept"]))
  min_lab_min_country = as.character(unlist(raneffs_D_lab[  like(lab, tolower(min_country)), "lab"][temp]))
  
  plotters <- data_subset[country %in% c(max_country, min_country, mid_country),
                          c("age", "gender", "country", "res", "sus")]
  
  ## Can take this out again? If using all labs and years.
  # group again (so ignorning lab and year)
  plotters_combo <- plotters[, sum(res), by =c("age", "gender", "country")]
  colnames(plotters_combo)[4] <- "res"
  plotters_combo_sus <- plotters[, sum(sus), by =c("age", "gender", "country")]
  colnames(plotters_combo_sus)[4] <- "sus"
  plotters_combo[plotters_combo_sus, on =c("age", "gender", "country"), sus := i.sus]
  plotters_combo[, prop := res/(res+sus)]
  plotters_combo <- left_join(plotters_combo, country_reference)
  
  age_s_total = length(seq(min_age/100, max_age/100, by = 0.01))
  
  
  
  data_test <- data.table(
    country = rep(c(rep(max_country, age_s_total*2),
                    rep(mid_country, age_s_total*2),
                    rep(min_country, age_s_total*2))),
    gender = c(rep("f", age_s_total*3*2),
               rep("m", age_s_total*3*2)),
    laboratorycode = rep(c(rep(min_lab_max_country, age_s_total),
                           rep(max_lab_max_country, age_s_total),
                           rep(min_lab_mid_country, age_s_total),
                           rep(max_lab_mid_country, age_s_total),
                           rep(min_lab_min_country, age_s_total),
                           rep(max_lab_min_country, age_s_total)),),
    age_s = c(rep(seq(min_age/100, max_age/100, by = 0.01),12)),
    year_s = 1,
    tot_samples = 1
  )
  
  
  data_test[,age_squared_s := age_s*age_s]
  
  # Anonymised countries
  data_test <- left_join(data_test, country_reference)
  
  model_to_predict <- Model_D
  
  model_to_predict <- fitted(model_to_predict, newdata = data_test)
  
  data_with_pred <- data.table(model_to_predict)
  data_with_pred <- cbind(data_test, model_to_predict)
  data_with_pred[, age := age_s*100]
  
  data_with_pred[year_s == 0, year := "2015" ]
  data_with_pred[year_s == 0.25, year := "2016" ]
  data_with_pred[year_s == 0.5, year := "2017" ]
  data_with_pred[year_s == 0.75, year := "2018" ]
  data_with_pred[year_s == 1, year := "2019" ]
  data_with_pred$year <- as.character(data_with_pred$year)
  
  data_with_pred[country_reference, on = c("country"), country_long := i.country_long]
  plotters_combo[country_reference, on = c("country"), country_long := i.country_long]
  
  # get countries in right order
  min_country_long <-country_reference[country == min_country]$country_long
  mid_country_long <-country_reference[country == mid_country]$country_long
  max_country_long <-country_reference[country == max_country]$country_long
  
  data_with_pred$country_long <- factor(data_with_pred$country_long,
                                        levels = c(min_country_long, mid_country_long, max_country_long))
  
  plotters_combo$country_long <- factor(plotters_combo$country_long,
                                        levels = c(min_country_long, mid_country_long, max_country_long))
  
  data_with_pred[gender == "f", gender := "female"]
  data_with_pred[gender == "m", gender := "male"]
  
  plotters_combo[gender == "f", gender := "female"]
  plotters_combo[gender == "m", gender := "male"]
  
  ############################################## Plot extreme countries 
  RIBBONS_FITTED2 <- ggplot(data_with_pred , aes(x = age, y = Estimate,
                                                 group= interaction(anon_country, gender, laboratorycode),
                                                 colour = gender)) +
    geom_point(data = plotters_combo, aes ( x = age, y = prop,
                                            colour = gender, group = NULL, linetype = NULL,
                                            size = (res+sus) ), alpha = 0.5)+
    geom_ribbon(aes(ymin =Q2.5, ymax =Q97.5, fill = gender,
                    group = interaction(country_long, gender, laboratorycode)),
                alpha = 0.5, colour = NA) +
    scale_colour_manual(values = colours_to_use[1:2]) +
    facet_grid(gender~anon_country) + 
    guides(colour = "none", fill = "none") +
    geom_line(linewidth=1) +
    labs(y= "Fitted  proportion", x = "age", shape = "country", title = paste0("Bacteria: ", pathogen_to_run, ", antibiotic: ", resistance_to_run),
         size = "Sample size", colour = "sex", fill = "sex") +
    lims(y = c(0,1)) 
  
  tiff(paste0("output_figures/by_bug_drug/", pathogen_to_run, "_", resistance_to_run, "_examples.tiff"),
       width = 3250, height = 2000, res = 300)
  print(RIBBONS_FITTED2)
  dev.off()
  
  
  if(only_model_D == F){
    legend = gtable_filter(ggplot_gtable(ggplot_build(LAB_INTERCEPTS)), "guide-box")
    
    tiff(paste0("output_figures/by_bug_drug/", pathogen_to_run, "_", resistance_to_run, "_combo.tiff"),
         width = 3250, height = 2000, res = 300)
    
    print( grid.arrange(
      COUNTRY_INTERCEPTS + theme(legend.position="none"),
      COUNTRY_AGE_INTERCEPTS + theme(legend.position="none"),
      RIBBONS_FITTED2, legend,
      layout_matrix = rbind(c(1,1,1,1,1,4,3,3,3,3,3,3,3),
                            c(2,2,2,2,2,4,3,3,3,3,3,3,3))
      
    ))
    
    dev.off()
  }
  
  
  legend = gtable_filter(ggplot_gtable(ggplot_build(RIBBONS_FITTED2)), "guide-box")
  
  ############################################## Plot all results combined: figure for paper 
  tiff(paste0("output_figures/by_bug_drug/combined_", pathogen_to_run, "_", resistance_to_run, "_combo.tiff"),
       width = 3850, height = 2250, res = 300)
  print( grid.arrange(
    models_table_all,
    ACROSS_COUNTRY + theme(legend.position="none")+ labs(title="B"),
    RIBBONS_FITTED2 + theme(legend.position="none") + labs(title="C") ,
    legend,
    layout_matrix = rbind(c(1,1,1,1,3,3,3,3,4),
                          c(2,2,2,2,3,3,3,3,4),
                          c(2,2,2,2,3,3,3,3,4)),
    top = textGrob(paste0("Bacteria: ", bug_drug[,1],", Antibiotic: ", bug_drug[,2]) ,gp=gpar(fontsize=20,font=3))
    
  ) )
  
  dev.off()
  
  
  plot_title <- paste0("RESULTS_", pathogen_to_run, "_", resistance_to_run )
  assign(plot_title, grid.arrange(
    models_table_all,
    ACROSS_COUNTRY + theme(legend.position="none")+ labs(title="B"),
    RIBBONS_FITTED2 + theme(legend.position="none") + labs(title="C") ,
    legend,
    layout_matrix = rbind(c(1,1,1,1,3,3,3,3,4),
                          c(2,2,2,2,3,3,3,3,4),
                          c(2,2,2,2,3,3,3,3,4)),
    top = textGrob(paste0("Bacteria: ", bug_drug[,1],", Antibiotic: ", bug_drug[,2]) ,gp=gpar(fontsize=18,font=3))
    
  ) )
  
  ############################################## add info to the big storage table
  models_overview <- as.data.frame(models_overview)
  models_overview$parameter <- rownames(models_overview)
  models_overview$bug <- unlist(bug_drug[,1])
  models_overview$drug <- unlist(bug_drug[,2])
  models_overview <- data.table(models_overview)
  
  all_params_out <- rbind(all_params_out, models_overview)
  
  print(paste0(pathogen_to_run, " ", resistance_to_run, " has been run. Number ", i, " of ",nrow(path_drug_combos)))
  
}
stop()

all_params_out[,Error := NULL]
all_params_out <- all_params_out[, c("bug", "drug", "parameter", "Estimate", "Lower 95%", "Upper 95%")]
all_params_out[, text_out := paste0(
  Estimate,
  " (",
  `Lower 95%`,
  " - ",
  `Upper 95%`,
  ")")]

all_params_out_c <- dcast.data.table(all_params_out, bug + drug ~ parameter, value.var = "text_out")


write.csv(all_params_out_c, file = "output/all_params_out.csv")


fixed_effs_storage <- rbindlist(fixed_effs_storage)
saveRDS(fixed_effs_storage, file = "output/fixed_effs_storage.RDS")
