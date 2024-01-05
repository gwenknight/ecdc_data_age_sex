### ********** ECDC data exploration and cleaning  ****************###############
##### The data was explored and cleaned to remove outliers and multiple antibiotic groupings. 
#### Refer to the supplementary for information on the steps / justification 

#### Load in libraries and set plotting rules
library(tidyverse)
theme_set(theme_bw(base_size = 11))
sensitivity_0_run_now <- "OFF" # This creates the cleaned data including 0 age if ON. Default should be OFF!

# checks if already defined globally. if not takes the local value
if(exists("sensitivity_0_run")){
  sensitivity_0 <- sensitivity_0_run
} else {
  sensitivity_0 <- sensitivity_0_run_now
}


##################********************* DATA ***************************#####################################
#### Read in and explore raw data
##################**********************************************************#####################################
#### Read in data
## Use original => no need to screen this as looking at totals
# data_orig <- read_csv("data/ecdc_ears_net_2002.csv")[,-1] #### CONTACT ECDC TO APPLY FOR THIS 
dim(data_orig) # Original = 3,549,617 isolates
str(data_orig) # Explore

# Make sure numeric values correct
data <- data_orig
data$age <- as.numeric(data$age)
data <- data %>% filter(!is.na(age)) %>% filter(gender %in% c("f","m"), age > 0)


### Plot raw data
data_long <- data %>% pivot_longer(cols = fq_ent_R:strpne_multi_R) %>% select(DateUsedForStatisticsYear, age, gender, pathogen, reportingcountry, name, value,laboratorycode) %>% 
  group_by(DateUsedForStatisticsYear, age, gender, pathogen, name, value) %>% summarise(n = n()) %>%
  pivot_wider(names_from = value, values_from = n) %>% mutate(proportion = `1`/(`1` + `0`)) %>% 
  filter(!is.na(proportion))

ggplot(data_long, aes(x = age, y = proportion, group = interaction(gender,DateUsedForStatisticsYear))) + geom_line(aes(col = DateUsedForStatisticsYear, linetype = gender)) +
  geom_point(aes(col = DateUsedForStatisticsYear),size = 0.1) +
  facet_grid(pathogen ~ name)
ggsave("plots/rawdata_all.png")

# Pivot data 
data_long_cntry <- data %>% pivot_longer(cols = fq_ent_R:strpne_multi_R) %>% select(year, age, gender, pathogen, reportingcountry, name, value) %>% 
  group_by(year, age, gender, pathogen, name, value, reportingcountry) %>% summarise(n = n()) %>%
  pivot_wider(names_from = value, values_from = n) %>% mutate(proportion = `1`/(`1` + `0`)) %>% 
  filter(!is.na(proportion))

##################********************* SUMMARY ***************************#####################################
#### What is the distribution of key characteristics?
##################**********************************************************#####################################

### What percentage of the isolates are not bloodstream but CSF? 
100 * table(data_orig$specimen)["csf"] / sum(table(data_orig$specimen)) # < 1% are CSF

# Laboratory code exploration 
which(is.na(data$laboratorycode)) # all had a laboratory code
lc <- data %>% group_by(reportingcountry, laboratorycode) %>% summarise(n = n())
summary(lc$n)
length(unique(lc$laboratorycode))

# Bacteria - antibiotic combinations
bug_drug <- unique(data_long_cntry[,c("pathogen","name")])
dim(bug_drug)

# Number of groupings
sum(table(data_long_cntry$gender))

# Number of results (0 or 1)
n_results <- data_long_cntry %>% summarise(totals = sum(`0` + `1`)) %>% ungroup() %>% summarise(sum(totals)) 
n_results
n_results / dim(data_orig)[1] # => 2.7 tests per sample 

# Number of m / f
table(data$gender)
100 * round(table(data$gender)[1]/sum(table(data$gender)),2)

# Summary of data
summ_data <- data_long_cntry %>% group_by(year, name, reportingcountry) %>% 
  summarise(total = sum(`0` + `1`))

# Add in label for bug-drug and bug-drug-yr combinations
bug_drug$combo <- paste0(bug_drug$pathogen, bug_drug$name)
data_long_cntry$combo <- paste0(data_long_cntry$pathogen, data_long_cntry$name)
data_long_cntry$comboyr <- paste0(data_long_cntry$pathogen, data_long_cntry$name, data_long_cntry$year)

##################********************* EXPLORATION ***************************#####################################
#### What is the distribution in data across each category? 
##################**********************************************************#####################################
#### (1) How many in total for each RESISTANCE
length(unique(data_long_cntry$name))
total_bugdrug <- data_long_cntry %>% group_by(name) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) %>% mutate(name=factor(name, levels=name)) %>% print(n=Inf) 
ggplot(total_bugdrug, aes(x=name, y = total)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# => bottom 3: acispp_multi_R // fq_strpne_R // cefIII_strpne_R

#### (2) Remove multi-resistance classifications
## Labelled as "combined" resistances in EARS-NET reports (e.g. Table 3 in 2019 report: https://www.ecdc.europa.eu/en/publications-data/surveillance-antimicrobial-resistance-europe-2019)
data %>% filter(pathogen == "strpne") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(strpne_multi_R == 1) # From EARS-NET: multi_R = combined R for strepne = Pen + macrolide R 
data %>% filter(pathogen == "strpne") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(strpne_multi_R == 1, penic_RI == 0, macrol_R == 0) # size = 0 => multi_R is a combination column so can be safely removed without affecting prevalence

data %>% filter(pathogen == "klepne") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(klepne_multi_R == 1) # From EARS-NET: multi_R = combined R for klepne = Fq + cefIII + aminogl R 
data %>% filter(pathogen == "klepne") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(klepne_multi_R == 1, fq_ent_R == 0, cefIII_entero_R == 0, aminogl_R == 0) # size = 0 => multi_R is a combination column so can be safely removed without affecting prevalence

# multi_R for esccol: Combined resistance to third-generation cephalosporins, fluoroquinolones, and aminoglycosides
data %>% filter(pathogen == "esccol") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(klepne_multi_R == 1) # size = 0 => multi_R is a bug specific thing (well one example to check)
data %>% filter(pathogen == "klepne") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(acispp_multi_R == 1) # size = 0 => multi_R is a bug specific thing (well one example to check)
data %>% filter(pathogen == "esccol") %>% select(pathogen, fq_ent_R:strpne_multi_R) %>% filter(esccol_multi_R == 1, fq_ent_R == 0, cefIII_entero_R == 0, aminogl_R == 0) # size = 0 => multi_R is a combination column so can be safely removed without affecting prevalence

# Other two: 
# multi_R for psaer: Combined resistance to >3 antimicrobial groups (among piperacillin + tazobactam, ceftazidime, carbapenems, fluoroquinolones and aminoglycosides)
# multi_R for acispp: Combined resistance to carbapenems, fluoroquinolones and aminoglycosides

multi_rs <- c("esccol_multi_R","strpne_multi_R","klepne_multi_R","pseaer_multi_R","acispp_multi_R")
bug_drug <- bug_drug %>% filter(!name %in% multi_rs)

# what was removed? 
data_multi_out <- data_long_cntry %>% filter(name %in% multi_rs)
100 * sum(data_multi_out$`0` + data_multi_out$`1`) / n_results # 12% of "results" removed by this ("results" as this was a column that was a summary of other columns)

#####********* FILTERING ******* #######
data_long_cntry1 <- data_long_cntry %>% filter(combo %in% bug_drug$combo)

#### (3) How many in total for each BACTERIA
total_bugdrug <- data_long_cntry1 %>% group_by(pathogen) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) %>% mutate(pathogen=factor(pathogen, levels=pathogen)) %>% print(n=Inf) 
ggplot(total_bugdrug, aes(x=pathogen, y = total)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
## Keeping acispp in for now but low numbers

#### (4) How many in total for each RESISTANCE + BACTERIA
total_bugdrug <- data_long_cntry1 %>% group_by(pathogen, name) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) %>% print(n=Inf)
ggplot(total_bugdrug, aes(x=name, y = total)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(~pathogen, scales="free")
ggplot(total_bugdrug, aes(x=reorder(interaction(pathogen,name),total), y = total, group = interaction(pathogen, name))) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
## All above 14K

#### (5) How many in total for each COUNTRY
total_bugdrug <- data_long_cntry1 %>% group_by(reportingcountry) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) %>% 
  mutate(reportingcountry=factor(reportingcountry, levels=reportingcountry)) %>% print(n=Inf) 
ggplot(total_bugdrug, aes(x=reportingcountry, y = total)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Explore
hist(total_bugdrug$total, breaks = seq(0,2000000,10000))

# => 10 isolates x 20 age groups x 19 years x 3 pathogens = 11400 threshold for exclusion? not for now
bot_3000cntry <- as.character(unlist(total_bugdrug %>% filter(total < 11400) %>% select(reportingcountry))) 
countries <- unique(as.character(unlist(data %>% dplyr::select(reportingcountry)))) # list of final countries

#### (6) How many in total for each YEAR
total_bugdrug <- data_long_cntry1 %>% group_by(year) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) %>% mutate(year=factor(year, levels=year)) %>% print(n=Inf) 
ggplot(total_bugdrug, aes(x=year, y = total)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# Don't remove any years - need to also split by bug+drug

## remove early years
data_long_cntry2 <- data_long_cntry1 %>% filter(year > 2015, year < 2020)

# what does this filtering do to the number of results in the data? 
n_results1 <- data_long_cntry1 %>% summarise(totals = sum(`0` + `1`)) %>% ungroup() %>% summarise(sum(totals)) 
n_results2 <- data_long_cntry2 %>% summarise(totals = sum(`0` + `1`)) %>% ungroup() %>% summarise(sum(totals)) 
n_results
n_results1
n_results2
100 * (n_results1 - n_results2)/n_results1
100 * n_results2/n_results


#### (7) How many in total for each YEAR + BUG-DRUG
total_bugdrug <- data_long_cntry2 %>% group_by(combo, year) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) 
tail(total_bugdrug)
ggplot(total_bugdrug, aes(x=year, y = total)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~combo, scales = "free")

ggplot(total_bugdrug, aes(x=year, y = total)) + geom_bar(stat = "identity", position = "dodge",aes(fill = combo)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~year, scales = "free")

ggplot(total_bugdrug, aes(x=total)) + geom_histogram(binwidth = 500) + scale_x_continuous("Number in a pathogen, antibiotic and year grouping")
ggsave("plots/appendix_data_cleaning1.pdf")

# Explore distribution 
h <- hist(total_bugdrug$total, breaks = seq(0,160000,500)) #break at < 3000 => two peaks 
ggplot(total_bugdrug, aes(x=total)) + geom_histogram(binwidth = 500) + scale_x_continuous("Number in a pathogen, antibiotic and year grouping", lim = c(0,30000))
h$counts
h$breaks
# Add column and use to explore if < 200 isolates (10 per 20 age groupings)
total_bugdrug$comboyr <- paste0(total_bugdrug$combo, total_bugdrug$year)
lowest_bdy <- as.character(unlist(total_bugdrug %>% ungroup() %>% filter(total < 200) %>% select(comboyr)))
length(lowest_bdy)

# output those with low numbers 
pp <- data_long_cntry2 %>% filter(comboyr %in% lowest_bdy) %>% ungroup() %>% select(pathogen, name, year)
u <- unique(pp[,c("pathogen","name","year")])
write.csv(u, "data/removed_path_drug_year.csv")

#### (8) How many in total for each YEAR + BUG-DRUG + AGE 
total_bugdrug <- data_long_cntry2 %>% filter(combo %in% bug_drug$combo) %>% group_by(combo, year, age) %>% summarise(total = sum(`0` + `1`)) %>% arrange(desc(total)) 
tail(total_bugdrug)
ggplot(total_bugdrug, aes(x=total)) + geom_histogram(binwidth = 2) + scale_x_continuous("Number in a pathogen, antibiotic, year and age grouping", lim = c(0,100))
# Many have < 25 

##################********************* CLEANING ***************************#####################################
#### Based on the above - remove / clean data 
##################**********************************************************#####################################
# Original with
data <- data_orig
# (a) filtered for data with age and gender values 
# (b) with all aged older than 0
data$age <- as.numeric(data$age)
if(sensitivity_0 == "OFF"){
  data <- data %>% filter(!is.na(age)) %>% filter(gender %in% c("f","m"), age > 0)
} else if(sensitivity_0 == "ON"){ # Include those aged 0 in cleaned data 
  data <- data %>% filter(!is.na(age)) %>% filter(gender %in% c("f","m"))
}

# (c) pivoted into long format and multi resistances removed 
# Multi resistances are: 
multi_rs <- c("esccol_multi_R","strpne_multi_R","klepne_multi_R","pseaer_multi_R","acispp_multi_R")
data_long_final0 <- data %>% 
  pivot_longer(cols = fq_ent_R:rifamp_R) %>%# just don't take last 5 columns = multi resistance
  select(DateUsedForStatisticsYear, age, gender, year, pathogen, reportingcountry, name, value,laboratorycode, patientcounter)

# (d) just 2015-2019
data_long_final1 <- data_long_final0 %>% 
  filter(year > 2014, year < 2020) ## NEW: 75958584 => 34048728   
rm(data_long_final0)

### No more filtering in the below: so this is the final number of patients
data_long_final1 %>% summarise(length(unique(patientcounter))) # 944520

data_long_final1 %>% group_by(gender) %>% summarise(length(unique(patientcounter))) # f = 444778, m = 538723

# (e) with proportions and combos
data_long_final2 <- data_long_final1 %>% 
  select(DateUsedForStatisticsYear, age, gender, pathogen, reportingcountry, name, value,laboratorycode) %>% # don't need patientcounter now
  group_by(DateUsedForStatisticsYear, age, gender, pathogen, name, value, laboratorycode, reportingcountry) %>% 
  summarise(n = n()) %>% # sum number of susceptibility results by each of these groupings
  pivot_wider(names_from = value, values_from = n)
rm(data_long_final1)

data_long_final2$`1`[is.na(data_long_final2$`1`)] <- 0
data_long_final2$`0`[is.na(data_long_final2$`0`)] <- 0

data_long_final <- data_long_final2 %>% mutate(proportion = `1`/(`1` + `0`)) %>% # calculate proportions
  filter(!is.na(proportion))

colnames(data_long_final) <- c("year", "age","gender","pathogen", "name", "laboratorycode","country", "sus", "res", "NA", "proportion")
sum(data_long_final$sus + data_long_final$res) # 6862577

# (f) add in bug drug names 
data_long_finalnames <- data_long_final %>% 
  mutate(pathogen = recode(pathogen, "encfae" = "Enterococcus faecalis", "encfai" = "Enterococcus faecium", 
                           "esccol" = "Escherichia coli", "klepne" = "Klebsiella pneumoniae", 
                           "staaur" = "Staphylococcus aureus", "strpne" = "Streptococcus pneumoniae", 
                           "acispp" = "Acinetobacter spp", "pseaer"="Pseudomonas aeruginosa"),
         drug = recode(name, "amika_R" = "amikacin","aminogl_R" = "aminoglycoside","carbapen_R" = "carbapenem",
                       "fq_pseudo_R" = "fluoroquinolone", "fq_staaur_R" = "fluoroquinolone", "fq_strpne_R" = "fluoroquinolone","fq_ent_R" = "fluoroquinolone",
                       "aminopen_R" = "aminopenicillins","genta_high" = "high_level_aminoglycoside","vanco_R" = "vancomycin",
                       "cefIII_entero_R" = "Third_generation_cephalosporins", "cefIII_strpne_R" = "Third_generation_cephalosporins",
                       "ert_R" = "ertapenem","ureidopen_R" = "piperacillin_tazobactam", "ceftaz_R" = "ceftazidime",
                       "mrsa_R" = "methicillin","rifamp_R"="rifampicin","macrol_R" = "macrolides","penic_RI" = "penicillins")) # check out table S3 in appendix

# SAVE
if(sensitivity_0 == "ON"){
  write.csv(data_long_finalnames, "data/data_cleaned_with0s.csv")
} else if(sensitivity_0 == "OFF"){
  write.csv(data_long_finalnames, "data/data_cleaned.csv")
}

#####********* Basic analysis of final data ******* #######
length(unique(data_long_finalnames$country))
length(unique(data_long_finalnames$pathogen))
length(unique(data_long_finalnames$name))
mean(data_long_finalnames$proportion)
sd(data_long_finalnames$proportion)
sum(data_long_finalnames$sus + data_long_final1$res)
n_results2

f <- data_long_finalnames %>% ungroup() %>% group_by(pathogen, name) %>% summarise(total = sus + res) %>%
  group_by(pathogen, name) %>% summarise(totals = sum(total)) 
write.csv(f, "data/summary_final_cleaned.csv")

f %>% group_by(pathogen) %>% summarise(sum = sum(totals)) %>% arrange(desc(sum))
f %>% mutate(bugdrug = paste0(pathogen,name)) %>% arrange(desc(totals)) %>% print(n=Inf)

