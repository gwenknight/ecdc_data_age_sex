# missing data

library(data.table)
library(tidyverse)
library(gridExtra)

run_subset <- NA #"staaur" #NA #  # or NA for whole

################## MISSING AGE AND GENDER############################
#### Read in data
## data_orig <- read_csv("data/ecdc_ears_net_2002.csv")[,-1] # Apply to TESSY for access
dim(data_orig) 
# Make sure numeric values correct
data <- data_orig
data$age <- as.numeric(data$age)
data <- as.data.table(data)
data <- data[DateUsedForStatisticsYear >= 2015]

name <- "all"

if(!is.na(run_subset)){
data <- data[pathogen == run_subset]
name <- run_subset
}

data_missing_gender <- data[!(gender %in% c("f","m"))]
# 156133 missing gender
data_missing_age <- data[is.na(age)]
# 70120 missing age
data_missing_both <- data[is.na(age) & !(gender %in% c("f","m"))]
# 41265 missing both

# numbers of missing
# missing age only
nrow(data_missing_age) - nrow(data_missing_both)
# missing gender only
nrow(data_missing_gender) - nrow(data_missing_both)
# missing both
nrow(data_missing_both)


data[!gender %in% c("f","m"), missing_gender := "Sex missing"]
data[gender %in% c("f","m"), missing_gender := "Sex included"]
data[, gender := as.factor(gender)]
# of those missing gender, is the distribution of age the same as those not missing gender
MISSING_GENDER <- ggplot(data, aes(x = age, colour = missing_gender, group = missing_gender)) + 
  geom_density(linewidth = 0.8) + 
  labs(x = "Age", y="Density", colour = "", title = paste0("A: ", name)) + 
  theme_linedraw() +
  theme(legend.position = "bottom") 
 

no_age <- prop.test(nrow(data[is.na(age) & gender == "f"]), 
                    nrow(data[is.na(age) & gender %in% c("f","m")]))
with_age <- prop.test(nrow(data[(!is.na(age)) & gender == "f"]), 
                      nrow(data[(!is.na(age)) & gender %in% c("f","m")]))

age_test <- data.frame( type = c("Age missing", "Age included"),
                        estimate = c(no_age[[4]], with_age[[4]]), 
                        min = c(no_age[[6]][[1]], with_age[[6]][[1]]),
                        max = c(no_age[[6]][[2]], with_age[[6]][[2]] ))

MISSING_AGE <- ggplot(age_test, aes( x = type, y = estimate, ymin = min, ymax = max)) + 
  geom_errorbar() + geom_point() +
  labs(x = "", y = "Proportion female", title =paste0("B: ", name)) + 
  lims(y = c(0, 1))+ 
  theme_linedraw()

tiff(paste0("plots/missing_data_",name,".tiff"), height = 2000, width = 3000, res = 300)
grid.arrange(MISSING_GENDER, MISSING_AGE, layout_matrix = rbind(c(1,1,1,2,2)))
dev.off()

# try to run a regression with the missing data as a variable outcome
data_dt <- as.data.table(data)

# choose one bug to look at for now. # mrsa_resistance. 
data_dt_mrsa <- data_dt[,c("DateUsedForStatisticsYear", "reportingcountry", "gender",
                      "hospitalid", "hospitalunittype", "laboratorycode", "pathogen", "mrsa_R", 
                      "age")]


##### Run a regression on the missingness (MRSA) ####
data_dt_mrsa <- data_dt_mrsa[pathogen == "staaur" & 
                          DateUsedForStatisticsYear %in% c("2019", "2018", 
                                                           "2017", "2016", 
                                                           "2015"),]

data_dt_mrsa[!(gender %in% c("m", "f")), gender := "unk"]
data_dt_mrsa[, age_s := age/100]
data_dt_mrsa[, age_squared_s := age_s *age_s]

data_dt_mrsa[DateUsedForStatisticsYear == "2015", year_s := 0]
data_dt_mrsa[DateUsedForStatisticsYear == "2016", year_s := 0.25]
data_dt_mrsa[DateUsedForStatisticsYear == "2017", year_s := 0.5]
data_dt_mrsa[DateUsedForStatisticsYear == "2018", year_s := 0.75]
data_dt_mrsa[DateUsedForStatisticsYear == "2019", year_s := 1]

data_dt_mrsa[, gender_missing := 0]
data_dt_mrsa[!(gender %in% c("f", "m")), gender_missing := 1]

data_dt_mrsa[, age_missing := 0]
data_dt_mrsa[is.na(age), age_missing := 1]

gender_missing <- glm(gender_missing ~
                        1  + 
                        age_s + 
                        gender ,
                        data = data_dt_mrsa, family = binomial(link = "logit"))
age_missing <- glm(age_missing ~
                        1  + 
                        age_s + 
                        gender ,
                      data = data_dt_mrsa, family = binomial(link = "logit"))

###### LOOK AT MISSINGNESS OF SUSCEPTIBILITY TESTING ####


data_dt_mrsa[is.na(mrsa_R), susceptibility_missing := "missing"]
data_dt_mrsa[!is.na(mrsa_R), susceptibility_missing := "not missing"]
data_dt_mrsa$susceptibility_missing <- as.factor(data_dt_mrsa$susceptibility_missing)

#
ggplot(data_dt_mrsa, aes(x = age, colour = susceptibility_missing, group = susceptibility_missing)) + 
  geom_density() + 
  labs(x = "Age", y = "Density", colour = "Susceptibility")


mrsa_missing <- prop.test(nrow(data_dt_mrsa[susceptibility_missing == "missing" & gender == "f"]), 
                    nrow(data_dt_mrsa[susceptibility_missing == "missing" & gender %in% c("f","m")]))
mrsa_not_missing <- prop.test(nrow(data_dt_mrsa[susceptibility_missing == "not missing" & gender == "f"]), 
                          nrow(data_dt_mrsa[susceptibility_missing == "not missing" & gender %in% c("f","m")]))

sus_missing <- data.frame( type = c("Susceptibility missing", "Susceptibility included"),
                        estimate = c(mrsa_missing[[4]], mrsa_not_missing[[4]]), 
                        min = c(mrsa_missing[[6]][[1]], mrsa_not_missing[[6]][[1]]),
                        max = c(mrsa_missing[[6]][[2]], mrsa_not_missing[[6]][[2]] ))

SUS_MISSNG_GENDER <- ggplot(sus_missing, aes(x = type, y = estimate, ymin = min, ymax = max)) + 
  geom_errorbar() + geom_point() +
  labs(x = "", y = "Proportion female", title =paste0("B")) + 
  lims(y = c(0, 1))

SUS_MISSING_AGE <- ggplot(data_dt_mrsa, aes(x = age, colour = susceptibility_missing, group = susceptibility_missing)) + 
  geom_density() + 
  labs(x = "Age", y = "Density", colour = "Susceptibility", title = "A") +
  theme(legend.position = "bottom")

tiff(paste0("plots/missing_susceptibility_mrsa.tiff"), height = 2000, width = 3000, res = 300)
grid.arrange(SUS_MISSING_AGE, SUS_MISSNG_GENDER, layout_matrix = rbind(c(1,1,1,2,2)))
dev.off()

######## Calculate the proportion of missing data by bug_drug combo. ######
translate_bug_drugs <-as.data.table(read.csv("data/translate_drug_bugs.csv", header = F))

# melt for subsetting
data_for_subsetting <- melt.data.table(data_dt, id.vars = c("DateUsedForStatisticsYear" ,"DateUsedForStatisticsQuarter","DateUsedForStatisticsMonth",  "DateUsedForStatisticsWeek"  , 
"DateUsedForStatisticsDay"    , "reportingcountry"           , "age"                      ,   "DateOfHospitalisationYear" ,  
"DateOfHospitalisationQuarter", "DateOfHospitalisationMonth" , "DateOfHospitalisationWeek",   "DateOfHospitalisationDay"  ,  
"dateofhospitalisation"       , "diskload"                   , "esbl"                     ,   "gender"                    ,  
"hospitalid"                  , "hospitalunittype"           , "laboratorycode"           ,   "pathogen"                  ,  
"patientcounter"              , "patienttype"                , "referenceguidelinessir"   ,   "resultcarbapenemases"      ,  
"serotype"                    , "specimen"                   , "dateusedforstatistics"    ,   "year"                      ,  
"bd"), measure.vars = c( "fq_ent_R"               ,     "fq_pseudo_R"          ,       "fq_strpne_R"          ,       
"fq_staaur_R"            ,     "cefIII_entero_R"      ,       "cefIII_strpne_R"      ,        "ceftaz_R"         ,           
"aminogl_R"              ,     "amika_R"              ,       "ureidopen_R"          ,        "carbapen_R"       ,           
"ert_R"                  ,     "mrsa_R"               ,       "penic_RI"             ,        "macrol_R"         ,           
"aminopen_R"             ,     "genta_high"           ,       "vanco_R"              ,        "rifamp_R"         ,           
"esccol_multi_R"         ,     "klepne_multi_R"       ,       "pseaer_multi_R"       ,        "acispp_multi_R"   ,           
"strpne_multi_R"))
  
  
#path_drug_combos <- data.frame(read.csv2("path_drug_combos.csv", sep = ",", header = F))
storage_missing <- data.frame(matrix(ncol = 3, nrow = nrow(path_drug_combos)))

for(i in 1:nrow(translate_bug_drugs)){
  
  pathogen_to_run <-  unlist(translate_bug_drugs[i,2])
  resistance_to_run <- unlist(translate_bug_drugs[i,3])

  bug_drug <- translate_bug_drugs[V2 == pathogen_to_run & V3 == resistance_to_run][,c("V1", "V6")]
  data_subset <- data_for_subsetting[pathogen == pathogen_to_run & variable == resistance_to_run]
  # calculate missingness 
  Susceptibility_missing <- nrow(data_subset[is.na(value)])
  total_samples <-nrow(data_subset)
  percent_not_tested <- (Susceptibility_missing/total_samples)*100
  
  temp <- c(pathogen_to_run, resistance_to_run, percent_not_tested)
  storage_missing[i,] <- temp
  }

colnames(storage_missing) <- c("pathogen", "resistance", "Percent_NA")
storage_missing$Percent_NA <- as.numeric(storage_missing$Percent_NA)
storage_missing$inverse_percent <- 100 - storage_missing$Percent_NA
storage_missing$Percent_NA <- round(storage_missing$Percent_NA, 1)
storage_missing$inverse_percent <- round(storage_missing$inverse_percent, 1)

  
