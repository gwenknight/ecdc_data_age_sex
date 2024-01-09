# make the country level plots
# note this takes a while to run (not written very efficiently)

# What bacteria-antibiotic combinations are considered? 
path_drug_combos <- data.frame(read.csv2("path_drug_combos.csv", sep = ",", header = F))
dir.create(file.path("output_figures", "country_specific"), showWarnings = FALSE)
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


# Read in the required data. 
data<- as.data.table(read_csv("data/data_cleaned.csv")[,-1])
country_reference <- as.data.table(read.csv("countries_reference.csv"))
colnames(country_reference) <- c("country", "country_long")
data[country_reference, on = c("country"), country_long := i.country_long]
translate_bug_drugs <-as.data.table(read.csv("data/translate_drug_bugs.csv", header = F))



# and loop through all bacteria-antibiotic combinations
for(j in 1:nrow(country_reference)){
  #counter only used for storage data frames
  counter <- 0
  for(i in 1:nrow(path_drug_combos)){
    
    target_country <- country_reference[j, "country_long"]
    target_country_short <- country_reference[j, "country"]
    # Combinations
    pathogen_to_run <-  path_drug_combos[i,1]  # choose from: c("acispp", "encfae", "encfai", "esccol", "klepne", "pseaer", "staaur", "strpne")
    resistance_to_run <- path_drug_combos[i,2]
    
    bug_drug <- translate_bug_drugs[V2 == pathogen_to_run & V3 == resistance_to_run][,c("V1", "V6")]
    data_subset <- data[pathogen == pathogen_to_run & name == resistance_to_run]
    
    # don't make if there's no data for that subset
    if(nrow(data_subset[country_long == target_country]) >0){
      #counter only used for storage data frames
      counter <- counter+1
      
      Model_D <- load_model(pathogen_to_run = pathogen_to_run, 
                            resistance_to_run = resistance_to_run, 
                            model_c = "D",
                            files_to_search = files_to_search)
      
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
      
      # ages
      min_age <- min(data_subset[country_long==target_country,]$age)
      max_age <- max(data_subset[country_long == target_country,]$age)
      
      temp = which.max(unlist(raneffs_D_lab[ like(lab, tolower(target_country_short)), "Estimate.Intercept"]))
      max_lab_country = as.character(unlist(raneffs_D_lab[  like(lab, tolower(target_country_short)), "lab"][temp]))
      temp = which.min(unlist(raneffs_D_lab[ like(lab, tolower(target_country_short)), "Estimate.Intercept"]))
      min_lab_country = as.character(unlist(raneffs_D_lab[ like(lab, tolower(target_country_short)), "lab"][temp]))
      
      plotters <- data_subset[country == target_country_short,
                              c("age", "gender", "country", "res", "sus")]
      
      # group again (so ignorning lab and year)
      plotters_combo <- plotters[, sum(res), by =c("age", "gender", "country")]
      colnames(plotters_combo)[4] <- "res"
      plotters_combo_sus <- plotters[, sum(sus), by =c("age", "gender", "country")]
      colnames(plotters_combo_sus)[4] <- "sus"
      plotters_combo[plotters_combo_sus, on =c("age", "gender", "country"), sus := i.sus]
      plotters_combo[, prop := res/(res+sus)]
      
      age_s_total = length(seq(min_age/100, max_age/100, by = 0.01))
      
      data_test <- data.table(
        country = rep(c(rep(unlist(target_country_short), age_s_total*2))),
        gender = c(rep("f", age_s_total*2),
                   rep("m", age_s_total*2)),
        laboratorycode = rep(c(rep(min_lab_country, age_s_total),
                               rep(max_lab_country, age_s_total)),),
        age_s = c(rep(seq(min_age/100, max_age/100, by = 0.01),4)),
        year_s = 1,
        tot_samples = 1
      )
      
      data_test[,age_squared_s := age_s*age_s]
      
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
      
      data_with_pred[, bug := bug_drug[,1]]
      data_with_pred[, drug := bug_drug[,2]]
      
      plotters_combo[, bug := bug_drug[,1]]
      plotters_combo[, drug := bug_drug[,2]]
      
      if(counter == 1){
        data_storage_pred <- data_with_pred
        data_storage_data <- plotters_combo
      } else{
        data_storage_pred <- rbind(data_storage_pred, data_with_pred)
        data_storage_data <- rbind(data_storage_data, plotters_combo)
      }
    }
  }
  
  data_storage_pred[country_reference, on = c("country"), country_long := i.country_long]
  data_storage_data[country_reference, on = c("country"), country_long := i.country_long]
  
  data_storage_pred[gender == "f", gender := "female"]
  data_storage_pred[gender == "m", gender := "male"]
  
  data_storage_data[gender == "f", gender := "female"]
  data_storage_data[gender == "m", gender := "male"]
  
  ############################################## Plot extreme countries 
  RIBBONS_F <- ggplot(data_storage_pred[gender == "female"] , aes(x = age, y = Estimate,
                                                                  group= interaction(bug, drug, gender, laboratorycode),
                                                                  colour = gender)) +
    geom_point(data = data_storage_data[gender=="female"], aes ( x = age, y = prop,
                                                                 colour = gender, group = NULL, linetype = NULL,
                                                                 size = (res+sus) ), alpha = 0.5)+
    geom_ribbon(aes(ymin =Q2.5, ymax =Q97.5, fill = gender,
                    group = interaction(bug, drug, country_long, gender, laboratorycode)),
                alpha = 0.5, colour = NA) +
    scale_colour_manual(values = colours_to_use[1]) +
    scale_fill_manual(values = colours_to_use[1]) +
    facet_wrap(.~bug +drug, ncol =9) +
    geom_line(linewidth=1) +
    labs(y= "Fitted  proportion", x = "age", shape = "country", title = paste0("Country: ", target_country, ", sex: female"),
         size = "sample size", colour = "sex", fill = "sex") +
    lims(y = c(0,1))+ 
    theme(strip.text.x = element_text(size = 13), 
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 15), 
          legend.text = element_text(size =15), 
          legend.title = element_text(size = 15), 
          title = element_text(size = 18))
  
  RIBBONS_M <- ggplot(data_storage_pred[gender == "male"] , aes(x = age, y = Estimate,
                                                                group= interaction(bug, drug, gender, laboratorycode),
                                                                colour = gender)) +
    geom_point(data = data_storage_data[gender=="male"], aes ( x = age, y = prop,
                                                               colour = gender, group = NULL, linetype = NULL,
                                                               size = (res+sus) ), alpha = 0.5)+
    geom_ribbon(aes(ymin =Q2.5, ymax =Q97.5, fill = gender,
                    group = interaction(bug, drug, country_long, gender, laboratorycode)),
                alpha = 0.5, colour = NA) +
    scale_colour_manual(values = colours_to_use[2]) +
    scale_fill_manual(values = colours_to_use[2]) +
    facet_wrap(.~bug +drug, ncol =9) +
    geom_line(linewidth=1) +
    labs(y= "Fitted  proportion", x = "age", shape = "country", title = paste0("Country: ", target_country, ", sex: male"),
         size = "sample size", colour = "sex", fill = "sex") +
    lims(y = c(0,1)) + 
    theme(strip.text.x = element_text(size = 13), 
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 15), 
          legend.text = element_text(size =15), 
          legend.title = element_text(size = 15), 
          title = element_text(size = 18))
  
  
  tiff(paste0("output_figures/country_specific/", target_country_short, ".tiff"),
       width = 7000, height = 9000, res = 300)
  grid.arrange(RIBBONS_F, RIBBONS_M, ncol=1)
  dev.off()
  
  print(paste0("country ", target_country, " is complete and figure saved."))
  
}
