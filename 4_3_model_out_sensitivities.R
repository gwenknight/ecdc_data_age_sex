# model output for sensitivity analyses

# run loos? (may need to do this for the first time using, if model fit loo haven't already been run)
run_loo_calcs <- "YES"

# libraries etc. 
source("3_a_load_needed.R")
# read in info for all sensitivities

# Read in the required data. 
data<- as.data.table(read_csv("data/data_cleaned.csv")[,-1])
# lookup for country names
country_reference <- as.data.table(read.csv("countries_reference.csv"))
colnames(country_reference) <- c("country", "country_long")
data[country_reference, on = c("country"), country_long := i.country_long]
# look up for bug and drug names
translate_bug_drugs <-as.data.table(read.csv("data/translate_drug_bugs.csv", header = F))

pathogen_to_run <- "staaur"  
resistance_to_run <- "mrsa_R"

# look up and subset
bug_drug <- translate_bug_drugs[V2 == pathogen_to_run & V3 == resistance_to_run][,c("V1", "V6")]
data_subset <- data[pathogen == pathogen_to_run & name == resistance_to_run]




####### PRIOR SENSITIVITY ANALYSIS ######

# find models
files_to_search <- list.files("fits")

file_no <- tail(grep(pattern = paste0("D_", pathogen_to_run, "_", resistance_to_run), 
                     x = files_to_search),1)

Model_default <- brm(file = paste0("fits/",files_to_search[file_no]))

file_no <- tail(grep(pattern = paste0("D_prior_", pathogen_to_run, "_", resistance_to_run), 
                     x = files_to_search),1)

Model_prior <- brm(file = paste0("fits/",files_to_search[file_no]))


# 
models_overview <- data.table(matrix(nrow = 12, ncol = 6))
colnames(models_overview) <- c("Parameter", "Estimate", "Error", "Lower", "Upper", "Model")
models_overview$Parameter <- rep(c("intercept", "year", "age", "age^2", "gender(m)", 
                               "Age:gender(m)"),2)
models_overview$Model <- rep(c("Default", "Regularising prior"), each = 6)

models_overview$Estimate <- c(c(fixef(Model_default)[c("Intercept", "year_s", 
                                                     "age_s", "age_squared_s", 
                                                     "genderm", "age_s:genderm"),"Est.Error"]), 
                              c(fixef(Model_prior)[c("Intercept", "year_s", 
                                                       "age_s", "age_squared_s", 
                                                       "genderm", "age_s:genderm"),"Est.Error"]))
models_overview$Error <- c(c(fixef(Model_default)[c("Intercept", "year_s", 
                                                "age_s", "age_squared_s", 
                                                "genderm", "age_s:genderm"),"Est.Error"]), 
                           c(fixef(Model_prior)[c("Intercept", "year_s", 
                                                  "age_s", "age_squared_s", 
                                                  "genderm", "age_s:genderm"),"Est.Error"]))
models_overview[, Lower := Estimate - 1.96*Error]
models_overview[, Upper := Estimate + 1.96*Error]

SENSITIVITY_PRIOR <- ggplot(models_overview, aes(x = Estimate, y = Parameter, colour = Model)) + 
  geom_pointrange(aes(xmin = Lower, xmax = Upper), position = position_dodge(width = 0.3))

tiff(paste0("output_figures/sensitivity_prior_mrsa.tiff"), 
     width = 3250, height = 2000, res = 300)
print(SENSITIVITY_PRIOR)
dev.off()

###### AGED 0 SENSITIVITY #######

files_to_search <- list.files("fits")

file_no <- tail(grep(pattern = paste0("D_sensitivity0_staaur_mrsa"), 
                     x = files_to_search),1)

Model_aged0 <- brm(file = paste0("fits/",files_to_search[file_no]))

models_overview <- data.table(matrix(nrow = 12, ncol = 6))
colnames(models_overview) <- c("Parameter", "Estimate", "Error", "Lower", "Upper", "Model")
models_overview$Parameter <- rep(c("intercept", "year", "age", "age^2", "gender(m)", 
                                   "Age:gender(m)"),2)
models_overview$Model <- rep(c("Default", "Including age 0"), each = 6)

models_overview$Estimate <- c(c(fixef(Model_default)[c("Intercept", "year_s", 
                                                       "age_s", "age_squared_s", 
                                                       "genderm", "age_s:genderm"),"Est.Error"]), 
                              c(fixef(Model_aged0)[c("Intercept", "year_s", 
                                                     "age_s", "age_squared_s", 
                                                     "genderm", "age_s:genderm"),"Est.Error"]))
models_overview$Error <- c(c(fixef(Model_default)[c("Intercept", "year_s", 
                                                    "age_s", "age_squared_s", 
                                                    "genderm", "age_s:genderm"),"Est.Error"]), 
                           c(fixef(Model_aged0)[c("Intercept", "year_s", 
                                                  "age_s", "age_squared_s", 
                                                  "genderm", "age_s:genderm"),"Est.Error"]))
models_overview[, Lower := Estimate - 1.96*Error]
models_overview[, Upper := Estimate + 1.96*Error]

SENSITIVITY_AGED0 <- ggplot(models_overview, aes(x = Estimate, y = Parameter, colour = Model)) + 
  geom_pointrange(aes(xmin = Lower, xmax = Upper), position = position_dodge(width = 0.3)) +
  labs(title="A") +
  theme(legend.position = "none")


ranefs_aged0 <- as.data.frame(ranef(Model_aged0)$country)[,c("Estimate.age_s", "Q2.5.age_s", 
                                                             "Q97.5.age_s")]
ranefs_aged0$country <- rownames(ranefs_aged0)
ranefs_aged0$model <- "Including age 0"

ranefs_default <- as.data.frame(ranef(Model_default)$country)[,c("Estimate.age_s", "Q2.5.age_s", 
                                                             "Q97.5.age_s")]
ranefs_default$country <- rownames(ranefs_default)
ranefs_default$model <- "Default"

combined_ranefs <- rbind(ranefs_aged0, ranefs_default)
combined_ranefs <- as.data.table(combined_ranefs)
temp <- combined_ranefs[model == "Default"]
new_order <- temp[order(Estimate.age_s, country)]$country
combined_ranefs$country <- factor(combined_ranefs$country, levels = new_order)
#to_plot_ranefs <-mutate(combined_ranefs, country = reorder(country, Estimate.age_s))

RANDEFS <- ggplot(combined_ranefs, aes(x = country, y = Estimate.age_s, colour = model)) +
  geom_pointrange(aes(ymin = Q2.5.age_s, ymax = Q97.5.age_s),
                  position = position_dodge(width = 0.7)) +
  labs( x = "Country", y = "Country-level age slope variation", colour = "Model", 
        title = "B")


tiff(paste0("output_figures/sensitivity_aged0_mrsa.tiff"), 
     width = 3250, height = 2000, res = 300)
grid.arrange(SENSITIVITY_AGED0, RANDEFS, layout_matrix = rbind(c(1,2,2)))
dev.off()

###### PHILOSOPHICAL SENSITIVITY #####


# find models
files_to_search <- list.files("fits")

# load all 4 models
file_no <- tail(grep(pattern = paste0("A_", pathogen_to_run, "_", resistance_to_run), 
                     x = files_to_search),1)
Model_A <- brm(file = paste0("fits/",files_to_search[file_no]))

file_no <- tail(grep(pattern = paste0("B_", pathogen_to_run, "_", resistance_to_run), 
                     x = files_to_search),1)
Model_B <- brm(file = paste0("fits/",files_to_search[file_no]))

file_no <- tail(grep(pattern = paste0("C_", pathogen_to_run, "_", resistance_to_run), 
                     x = files_to_search),1)
Model_C <- brm(file = paste0("fits/",files_to_search[file_no]))

file_no <- tail(grep(pattern = paste0("D_", pathogen_to_run, "_", resistance_to_run), 
                     x = files_to_search),1)
Model_D <- brm(file = paste0("fits/",files_to_search[file_no]))

# Create table of output
#### Make a summary table of the model output
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
temp_table <- signif(models_overview, digits = 2)
temp_table[is.na(temp_table)] <- "-"
models_table <-tableGrob(temp_table)

title <- textGrob("A",
                  gp=gpar(fontsize=15))
padding <- unit(8,"mm")

models_table_all <- gtable_add_rows(
  models_table, 
  heights = grobHeight(title) + padding,
  pos = 0 )
models_table_all <- gtable_add_grob(
  models_table_all, 
  title, 
  1, 1, 1)

# only run the loos for the first time!
if(run_loo_calcs == "YES"){
  Model_A <- add_criterion(Model_A,
                criterion = c("loo"),
                cores = 1,
                seed = 12345,
                ndraws =3000,
                overwrite = TRUE
  )

  Model_B <- add_criterion(Model_B,
                criterion = c("loo"),
                cores = 1,
                seed = 12345,
                ndraws = 3000,
                overwrite = TRUE
  )

  Model_C <- add_criterion(Model_C,
                criterion = c("loo"),
                cores = 1,
                seed = 12345,
                ndraws = 3000,
                overwrite = TRUE
  )

  Model_D <- add_criterion(Model_D,
                criterion = c("loo"),
                cores = 1,
                seed = 12345,
                ndraws = 3000,
                overwrite = TRUE
  )

}

loo_compare(Model_A, Model_B, Model_C, Model_D, criterion = "loo")

model_weights(Model_A, Model_B, Model_C, Model_D, weights = "loo")
  
ll_loo <- loo_compare(Model_A, Model_B, Model_C, Model_D, criterion = "loo")
temp <- as.data.table(ll_loo) %>% mutate(model = rownames(ll_loo))

LOOS_AGE <- ggplot(temp[model %in% c("Model_A", "Model_B", "Model_D")], 
               aes(x = model, y = looic)) + geom_point() + 
  geom_errorbar(aes(ymin=looic - se_looic, ymax = looic + se_looic)) + coord_flip() + 
  labs(y = "LOOIC", x = "Model", title = "B")

LOOS_GENDER <- ggplot(temp[model %in% c("Model_A", "Model_C", "Model_D")], 
                   aes(x = model, y = looic)) + geom_point() + 
  geom_errorbar(aes(ymin=looic - se_looic, ymax = looic + se_looic)) + coord_flip() + 
  labs(y = "LOOIC", x = "Model", title = "C")
  
tiff(paste0("output_figures/sensitivity_philosophy.tiff"), 
     width = 3250, height = 2000, res = 300)
grid.arrange(models_table_all, LOOS_AGE, LOOS_GENDER, 
             layout_matrix= rbind(c(1,1), 
                                  c(2,3)))
dev.off()


