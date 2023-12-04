#### Incidence curves

### Use the EARS-NET data to calculate the incidence of infection with certain bacteria over time by age and sex
## This is the incidence of BSI (and CSF infection (<10%)).
## Will need to be inflated by the coverage values in the ECDC reports to give national values.  

## Aim: to explore the trends in incidence of infection by age per country over time 

## Libraries
library(tidyverse)
library(patchwork)
library(here)
library(RColorBrewer)
library(lme4)
library(broom)
theme_set(theme_bw(base_size = 16))

########################################################***** DATA *********#################################################################################
#########****************** Read in data on coverage of surveillance + infection numbers + population size  ******************###################################
##################################################################################################################################################################

#########****************** Data on number of infections: from TESSY ******************###################
#### Read in data
## Use original => no need to screen this as looking at totals
# data_orig <- read_csv("data/ecdc_ears_net_2002.csv")[,-1] #### CONTACT ECDC TO APPLY FOR THIS 

dim(data_orig) # Original = 3,549,617 isolates x 54 columns
# Make sure numeric values correct
data_inf <- data_orig %>% filter(!is.na(age)) %>% rename(sex = gender) %>% filter(sex %in% c("f","m"), age > 0)

data_inf$age <- as.numeric(data_inf$age)
data_inf$year <- as.numeric(data_inf$year)
data_inf <- rename(data_inf, country = reportingcountry)

# How much of the original data used? 
100 * dim(data_inf)[1] / dim(data_orig)[1]
100 * dim(data_orig %>% filter(is.na(age)))/ dim(data_orig)[1]
100 * dim(data_orig %>% filter(!gender %in% c("f","m")))/ dim(data_orig)[1]
dim(data_inf)[1]

# Summarise data = number of isolates per year / country in each age / sex for each pathogen: NOT just number tested for resistance 
data_infs <- data_inf %>% dplyr::select(year, age, country, sex, pathogen, patienttype,hospitalunittype) %>% 
  group_by(year, country, age, sex, pathogen, patienttype, hospitalunittype) %>% 
  summarise(n = n())

#########****************** Population coverage estimates ******************###################
######### Use ECDC estimates of population coverage to push these numbers up?
## Country specific estimates if they exist for 2018,2019,2020
data_popcov <- read_csv("data/pop_rep_ecdc.csv") 
ggplot(data_popcov , aes(x=country, y = as.numeric(est_pop_cov_perc))) + geom_point(aes(col = factor(year))) + 
  geom_line() + coord_flip() + 
  scale_color_discrete("Year") + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Estimated population coverage")
ggsave("plots/pop_coverage_estimates_1820.pdf")

popcovsd <- data_popcov %>% group_by(country) %>% summarise(sd_yr = sd(est_pop_cov_perc, na.rm = TRUE), diff_yr = max(est_pop_cov_perc, na.rm = TRUE) - min(est_pop_cov_perc, na.rm = TRUE))
ggplot(popcovsd, aes(x=country, y = as.numeric(sd_yr))) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  ggtitle("Standard deviation between years")
ggplot(popcovsd, aes(x=country, y = as.numeric(diff_yr))) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  ggtitle("Biggest difference between years")

###### Cassini data 
cass_data <- read_csv("data/est_pop_cov_cassini_tables3.csv") 

# Explore 
cass_data_long <- cass_data %>% pivot_longer(cols = "S. pneumoniae":"Acinotobacter spp", names_to = "bacteria",values_to = "est_cov")

ggplot(cass_data_long, aes(x=Country, y = est_cov)) + 
  geom_bar(stat = "identity", position = "dodge",aes(fill = bacteria)) + 
  geom_line() + coord_flip() + 
  scale_fill_brewer(palette = "Accent","Bacteria", labels = expression(italic("Acinotobacter spp"),italic("E coli"),italic("Enterococci"),italic("K pnuemoniae"),italic("P aeruginosa"),italic("S. aureus"),italic("S. pneumoniae"))) + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Estimated population coverage")
ggsave("plots/pop_coverage_estimates_cassini_SUPPFIG.pdf")

# How likely between bug variation? 
bug_var <- cass_data %>% select(-Country) %>% #remove country so can just compare coverage
  rowwise %>% # by row
  mutate(n_diff = n_distinct(unlist(cur_data()))) # number different
plot(table(bug_var$n_diff))
round(100 * table(bug_var$n_diff) / dim(cass_data)[1],0) # 70% have same coverage for all bacteria

max_var <- cass_data %>% select(-Country) %>% #remove country so can just compare coverage
  rowwise %>% # by row
  mutate(max_diff = max(unlist(cur_data())) - min(unlist(cur_data()))) %>%
  arrange(desc(max_diff)) %>% 
  left_join(cass_data) # add back in country
max_var # only a few have big differences

### Combine all pop coverage data 
cass_data_long <- rename(cass_data_long, country = Country)
cass_data_long_adj <- cass_data_long %>% group_by(country) %>% summarise(est_pop_cov_perc = mean(est_cov)) %>% mutate(year = 2015) 
combined <- rbind(data_popcov[,c("country", "est_pop_cov_perc","year")], cass_data_long_adj)

ggplot(combined , aes(x=country, y = as.numeric(est_pop_cov_perc))) + geom_point(aes(col = factor(year)), size = 3) + 
  geom_line() + coord_flip() + 
  scale_color_discrete("Year") + 
  geom_point(data = combined %>% filter(year == 2015), pch = "x", size = 5) + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Estimated population coverage")
ggsave("plots/pop_coverage_estimates_1520_SUPPFIG.pdf")

### Is there a trend to increasing coverage with time? 
n_data = combined %>% ungroup() %>% group_by(country) %>% summarise(n = n_distinct(year))
combined_n <- left_join(combined,n_data) 
fit <- lmList(est_pop_cov_perc ~ year | country, data=combined_n %>% filter(!is.na(est_pop_cov_perc)))
r2 <- sapply(fit,function(x) summary(x)$r.squared)
s <- sapply(fit,function(x) tidy(x))
u <- unique(combined$country)
table_store<- c()

for(i in u){
  if(i %in% colnames(s)){
    w_col = which(colnames(s) == i)
    table_store <- rbind(table_store, 
                         c(i, round(s[,w_col][[2]],4),signif(s[,w_col][[5]][2],3),round(r2[which(names(r2)==i)],4)))
  }
}

colnames(table_store) <- c("country","intercept","value","pval_value","Rsquared")
table_store <- as.data.frame(table_store) %>% filter(!is.na(c(value))) %>% filter(!value == 0, !is.nan(pval_value)) %>% 
  filter(Rsquared > 0.6, Rsquared < 1, pval_value < 0.01)

table_store # only in Italy! 

## Export summary data 
final_pop_covs <- combined %>% filter(!is.na(est_pop_cov_perc))
write.csv(final_pop_covs,"data/est_pop_cov_final.csv")
# Use year based coverage where exists
# Use last year for projections

#########****************** Inflate number of infections to reflect population coverage ******************###################
key <- read.csv("data/country_key.csv")
final_pop_covs <- final_pop_covs %>% mutate(inflation_factor = 100/est_pop_cov_perc)
final_pop_covs[which(final_pop_covs$country == "UK"),"country"] <- "United Kingdom"
pop_covs_key <- left_join(final_pop_covs, key) %>% select(-c(country))
pop_covs_key <- rename(pop_covs_key, country = code)

data_infs_inflate <- left_join(data_infs, pop_covs_key) %>% filter(!is.na(est_pop_cov_perc)) %>%
  mutate(n_inflat = round(n * inflation_factor,0))


###########################****************** Data on population sizes: from World Bank ******************###################
# key = projection,unit,sex,age,geo\time
data_popsize_percs <- read.csv("data/world_bank_data_pop.csv") %>% 
  separate(Series.Code, into = c("projection", "unit","agerange","sex","agegrouping")) %>% 
  mutate(lowage = substr(agerange,1,2), upage1 = substr(agerange,3,4),
         upage = ifelse(upage1 == "UP",120,upage1)) %>% 
  pivot_longer(cols =c(X2000:X2021), names_to = "xyear",values_to = "perc_pop") %>% 
  mutate(year = as.numeric(sub('.', '', xyear)),
         sex = tolower(substr(sex,1,1))) %>% 
  select(-c(upage1, xyear,agerange,agegrouping,projection,unit,Series.Name,Country.Code))

data_popsize_ns_orig <- read.csv("data/world_bank_data_pop2.csv")
data_popsize_ns <- data_popsize_ns_orig %>% filter(!sex == "t") %>% 
  pivot_longer(cols =c(X2000:X2021), names_to = "xyear",values_to = "popsize") %>% 
  mutate(year = as.numeric(sub('.', '', xyear))) %>% 
  select(-c(xyear,Series.Name,Country.Code))

data_popsize <- left_join(data_popsize_ns, data_popsize_percs) %>% 
  mutate(n_inage_gdr = perc_pop * popsize / 100) %>% 
  filter(!is.na(n_inage_gdr)) 
data_popsize$lowage <- as.numeric(data_popsize$lowage)
data_popsize$upage <- as.numeric(data_popsize$upage)

ggplot(data_popsize, aes(x=lowage, y = n_inage_gdr, group = year)) + geom_line(aes(colour = factor(year))) +
  facet_grid(Country.Name~sex, scales = "free") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Number in age group") + 
  scale_color_discrete("Year") +
  theme(strip.text.y = element_text(angle = 0))
ggsave("plots/population_by_age_sex_2000_2021_SUPPFIG.pdf", width = 14, height = 20)

mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(length(unique(data_popsize$lowage)))
ggplot(data_popsize, aes(x=year, y = n_inage_gdr, group = year)) + 
  geom_bar(stat = "identity",aes(fill = factor(lowage)), position = "fill") + 
  scale_fill_manual(values = mycolors, "Age") + 
  facet_grid(Country.Name~sex, scales = "free") + 
  scale_x_continuous("Year") + 
  scale_y_continuous("Proportion in each age group") + 
  theme(strip.text.y = element_text(angle = 0))
ggsave("plots/population_by_age_sex_bar_2000_2021_SUPPFIG.pdf", width = 14, height = 20)

# Do they match? check for random year. Yep for 2005
data_popsize %>% group_by(Country.Name, year) %>% summarise(total = sum(n_inage_gdr)) %>% filter(year == 2005)
data_popsize_ns_orig %>% filter(sex == "t") %>% select(c("Country.Name","X2005"))

# Rename columns 
data_popsize <- rename(data_popsize, country = Country.Code.EARS)
data_popsize$lowage <- as.numeric(data_popsize$lowage)
data_popsize$upage <- as.numeric(data_popsize$upage)
write.csv(data_popsize, "data/data_popsize_clean.csv")

########################################################***** COMBINE *********#################################################################################
#########****************** COMBINE coverage + infection numbers + population size to estimate incidence ******************###################################
##################################################################################################################################################################

#########*# Label pathogens
data_infs_inflate <- data_infs_inflate %>% mutate(pathogen = recode(pathogen, "encfae" = "Enterococcus faecalis", "encfai" = "Enterococcus faecium", 
                                               "esccol" = "Escherichia coli", "klepne" = "Klebsiella pneumoniae", 
                                               "staaur" = "Staphylococcus aureus", "strpne" = "Streptococcus pneumoniae", 
                                                           "acispp" = "Acinetobacter spp", "pseaer"="Pseudomonas aeruginosa"))


#########*
## to match age - put in age groupings
data_infs_inflate <- data_infs_inflate %>% mutate(age_group = cut(age, breaks = seq(-1,120,5)))
data_infs_inflate$age_group <- as.character(data_infs_inflate$age_group)
data_infs_inflate[which(data_infs_inflate$age > 80),"age_group"] <- "(80,120]"

data_popsize <- data_popsize %>% mutate(age_f = 0.5 * (lowage + upage), age_group = cut(age_f, breaks = seq(-1,120,5)))
data_popsize$age_group <- as.character(data_popsize$age_group)
data_popsize[which(data_popsize$lowage == 80),"age_group"] <- "(80,120]"

## Which countries? 
pc <- unique(data_popsize$country)
pi <- unique(data_infs_inflate$country)
length(intersect(pc,pi))
setdiff(pc,pi)
dput(unique(data_popsize$Country.Name))
data_infs_inflate %>% group_by(year) %>% summarise(n = n_distinct(country)) # different number of countries in each 

## Sensitivity: put min coverage and use this 
data_infs_inflate <- data_infs_inflate %>% group_by(country) %>% mutate(sens_min = min(est_pop_cov_perc), inflation_factor_max = 100/sens_min, n_sens = inflation_factor_max * n)

# Plot subset of data to check min extracted 
#ggplot(data_infs_inflate %>% filter(age < 10), aes(x=country, y = est_pop_cov_perc)) + geom_point(aes(col = year)) + geom_point(data = data_infs_inflate %>% filter(age < 10), aes(y = sens_min), pch = "x", size = 3)


## Store inflated data
write.csv(data_infs_inflate, "data/data_infs_inflated.csv")


############ Combine - with just age and sex split and bacteria: MAIN ANALYSIS AND SPLIT at EUROPEAN LEVEL 
# Colour palette for bacteria
brewer_col = c("#d73027","#f46d43","#fdae61","darkgoldenrod","#74add1","#4575b4","purple","violet")

data_popsize_eu <- data_popsize %>% group_by(year, sex, lowage, age_group) %>% summarise(total_inage_gdr = sum(n_inage_gdr)) # sum populations over country
data_infs_inflate_eu <- data_infs_inflate %>% filter(age > 0) %>% group_by(year,sex,pathogen,age_group) %>% summarise(total_inf = sum(n_inflat))
combine_eu = left_join(data_infs_inflate_eu,data_popsize_eu, by = c("year","sex","age_group")) %>% 
  mutate(incidence = (total_inf / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

# Store 
write.csv(combine_eu, "output/pathogen_age_incidence_estimates.csv")

# check
combine_eu %>% filter(sex == "m", lowage == 15, pathogen == "Escherichia coli")
# relabel 
combine_eu[which(combine_eu$sex == "f"),"sex"] <- "Female" 
combine_eu[which(combine_eu$sex == "m"),"sex"] <- "Male" 
combine_eu$pathogen <- factor(combine_eu$pathogen, levels = c("Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae", "Enterococcus faecium", 
                                                              "Pseudomonas aeruginosa", "Enterococcus faecalis",  
                                                              "Streptococcus pneumoniae","Acinetobacter spp"))

## Interesting younger age patterns 
## Interesting younger age patterns 
g2mid<- ggplot(combine_eu %>% filter(lowage < 55, year == 2019), 
               aes(x = lowage, y = incidence, group = sex,col = factor(sex), 
                   fill = factor(sex))) + 
  geom_point() + 
  geom_smooth(aes(linetype = sex), alpha = 0.2) + 
  facet_wrap(~pathogen, scales = "free", ncol = 2) + 
  scale_fill_manual("Sex", labels = c("Female","Male"), values  = c("#E41A1C", "#377EB8")) +  
  scale_color_manual("Sex", labels = c("Female","Male"), values = c("#E41A1C", "#377EB8")) +  
  scale_linetype_manual("Sex", labels = c("Female","Male"), values = c(1,2), guide = guide_legend(override.aes = list(linewidth = 1))) +  
  theme(strip.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_log10("Incidence (infections / 100,000 population per year)") + 
  theme(strip.background = element_rect(fill= "white", size=1.5))
ggsave("plots/incidence_sex_sep_combine_eu_mid_years.pdf")

## All incidence
g_group <- ggplot(combine_eu %>% filter(year == 2019), aes(x = lowage, y = incidence, group = pathogen)) + 
  geom_point(aes(col = factor(pathogen))) + 
  geom_smooth(aes(fill = factor(pathogen),col = factor(pathogen)), alpha = 0.2) + 
  facet_wrap(~sex,  ncol = 2) + 
  scale_fill_manual("Bacteria", values = brewer_col) + 
  scale_color_manual("Bacteria",values = brewer_col) + 
  theme(legend.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_log10("Incidence (infections / 100,000 population per year)") + 
  theme(strip.background = element_rect(fill= "white", size=1.5))
ggsave("plots/incidence_sex_group_combine_eu_all.pdf")

g2mid + g_group + plot_layout(widths=c(1,1.9)) + plot_annotation(tag_levels = "A")
ggsave("plots/Fig2.tiff", dpi = 300, width = 19, height = 11)


########################################################***** SENSITIVITY *********#################################################################################
#########****************** USE MIN coverage over 2015-2020 to inflate incidence           ******************###################################
##################################################################################################################################################################

data_popsize_eu <- data_popsize %>% group_by(year, sex, lowage, age_group) %>% summarise(total_inage_gdr = sum(n_inage_gdr)) # sum populations over country
## HERE sum over n_sens: number using the smallest coverage value
data_infs_inflate_eu <- data_infs_inflate %>% filter(age > 0) %>% group_by(year,sex,pathogen,age_group) %>% summarise(total_inf = sum(n_sens))
combine_eu_sens = left_join(data_infs_inflate_eu,data_popsize_eu, by = c("year","sex","age_group")) %>% 
  mutate(incidence = (total_inf / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

# relabel 
combine_eu_sens[which(combine_eu_sens$sex == "f"),"sex"] <- "Female" 
combine_eu_sens[which(combine_eu_sens$sex == "m"),"sex"] <- "Male" 
combine_eu_sens$pathogen <- factor(combine_eu_sens$pathogen, levels = c("Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae", "Enterococcus faecium", 
                                                              "Pseudomonas aeruginosa", "Enterococcus faecalis",  
                                                              "Streptococcus pneumoniae","Acinetobacter spp"))

# compare
combine_eu %>% filter(sex == "Male", lowage == 15, pathogen == "Escherichia coli")
combine_eu_sens %>% filter(sex == "Male", lowage == 15, pathogen == "Escherichia coli")


## Interesting younger age patterns 
g2mids<- ggplot(combine_eu_sens %>% filter(lowage < 55, year == 2019), aes(x = lowage, y = incidence, group = sex)) + 
  geom_point(aes(col = factor(sex))) + 
  geom_smooth(aes(fill = factor(sex),col = factor(sex)), alpha = 0.2) + 
  facet_wrap(~pathogen, scales = "free", ncol = 2) + 
  scale_fill_discrete("Sex") + 
  scale_color_discrete("Sex") + 
  theme(strip.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")
ggsave("plots/incidence_sex_sep_combine_eu_sens_mid_years.pdf")

## All incidence
g_groups <- ggplot(combine_eu_sens %>% filter(year == 2019), aes(x = lowage, y = incidence, group = pathogen)) + 
  geom_point(aes(col = factor(pathogen))) + 
  geom_smooth(aes(fill = factor(pathogen),col = factor(pathogen)), alpha = 0.2) + 
  facet_wrap(~sex,  ncol = 2) + 
  scale_fill_manual("Bacteria", values = brewer_col) + 
  scale_color_manual("Bacteria",values = brewer_col) + 
  theme(legend.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")
ggsave("plots/incidence_sex_group_combine_eu_sens_all.pdf")

g2mids + g_groups + plot_layout(widths=c(1,1.9)) + plot_annotation(tag_levels = "A")
ggsave("plots/incidence_SUPPFIGURE.pdf", width = 19, height = 11)

# ((g2mid + g_group + plot_layout(widths=c(1,1.9)) ) / (g2mids + g_groups + plot_layout(widths=c(1,1.9))) + plot_annotation(tag_levels = "A"))
# ggsave("plots/incidence_SUPPFIGURE.pdf", width = 19, height = 11)

# Compare
combine_eu_sens <- combine_eu_sens %>% rename(incidence_s = incidence)
left_join(combine_eu[,c("year","sex","pathogen","age_group","incidence")], combine_eu_sens[,c("year","sex","pathogen","age_group","incidence_s")]) %>% 
  mutate(perc_diff = 100 * (incidence_s - incidence) / incidence) %>% 
  group_by(sex, pathogen) %>% 
  summarise(mean(perc_diff))

left_join(combine_eu[,c("year","sex","pathogen","age_group","incidence")], combine_eu_sens[,c("year","sex","pathogen","age_group","incidence_s")]) %>% 
  mutate(perc_diff = 100 * (incidence_s - incidence) / incidence) %>% 
  group_by(sex, pathogen) %>% 
  summarise(mean = mean(perc_diff)) %>% summarise(mean(mean))

########################################################***** NUMBER OF R infections *********#################################################################################
#########****************** Combine with resistance prevalence to get number of R infections by age and sex and country ******************###################################
##################################################################################################################################################################
######### Combine - with country, age and sex split
incidence_est_data_cntryageglevel = left_join(data_infs_inflate %>% filter(age > 0),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, country, sex, pathogen, Country.Name, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  mutate(incidence = (total / n_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence)) 

write.csv(incidence_est_data_cntryageglevel, "output/country_pathogen_age_incidence_estimates.csv")

##### Resistance data 
amr_data_long_orig <- read_csv("data/data_cleaned.csv")

# relabel pathogens and add age groupings
amr_data_long <- amr_data_long_orig %>% rename(sex = gender) %>%
  mutate(pathogen = recode(pathogen, "encfae" = "Enterococcus faecalis", "encfai" = "Enterococcus faecium", 
                           "esccol" = "Escherichia coli", "klepne" = "Klebsiella pneumoniae", 
                           "staaur" = "Staphylococcus aureus", "strpne" = "Streptococcus pneumoniae", 
                           "acispp" = "Acinetobacter spp", "pseaer"="Pseudomonas aeruginosa")) %>%
  mutate(age_group = cut(age, breaks = seq(-1,120,5))) %>% # Add age groupings 
  filter(age > 0) 
amr_data_long$age_group <- as.character(amr_data_long$age_group)
amr_data_long[which(amr_data_long$age > 80),"age_group"] <- "(80,120]"
amr_data_long_sumcag <- amr_data_long %>% group_by(sex, age_group, combo, pathogen, name, country, year) %>% 
  summarise(sust = sum(sus), rest = sum(res), proportion = rest / (sust + rest))

################ Combine ###############
combine_inc_r <- left_join(amr_data_long_sumcag, incidence_est_data_cntryageglevel) %>% 
  mutate(total_r_cases = proportion * n_inage_gdr * incidence / 100000, 
         incidence_r = (total_r_cases / n_inage_gdr) * 100000) #%>% 
  #group_by(combo, country, sex, year, a) %>% 
  #filter(year == 2019) %>% 
  #summarise(tots = sum(total_r_cases, na.rm = TRUE))
write.csv(combine_inc_r, "output/combine_incidence_R_country_level.csv")

ggplot(combine_inc_r %>% filter(year == 2019), aes(x = lowage, y = incidence_r, group = interaction(sex,pathogen))) + 
  geom_point(aes(col = country)) + 
  geom_smooth(aes(fill = country,col = country), alpha = 0.2) + 
  facet_grid(pathogen~sex,  scales = "free") + 
  #scale_fill_manual("Bacteria", values = brewer_col) + 
  #scale_color_manual("Bacteria",values = brewer_col) + 
  theme(legend.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")

########################################################***** NUMBER OF R infections *********#################################################################################
#########****************** Combine with resistance prevalence to get number of R infections by age and sex only ******************###################################
##################################################################################################################################################################
## R data not by country
amr_data_long_sumag <- amr_data_long %>% group_by(sex, age_group, combo, pathogen, name, year) %>% 
  summarise(sust = sum(sus), rest = sum(res), proportion = rest / (sust + rest))

################ Combine ###############
amr_data_long_sumag[which(amr_data_long_sumag$sex == "f"),"sex"] <- "Female" 
amr_data_long_sumag[which(amr_data_long_sumag$sex == "m"),"sex"] <- "Male" 
combine_inc_r <- left_join(amr_data_long_sumag, combine_eu) %>% 
  mutate(total_r_cases = proportion * total_inage_gdr * incidence / 100000, 
         incidence_r = (total_r_cases / total_inage_gdr) * 100000) 
write.csv(combine_inc_r, "output/combine_incidence_R_level.csv")

ggplot(combine_inc_r %>% filter(year == 2019), aes(x = lowage, y = incidence_r, 
                                                   group = interaction(combo,pathogen))) + 
  geom_point(aes(col = factor(name))) + 
  geom_smooth(aes(fill = factor(name),col = factor(name)), alpha = 0.2) + 
  facet_grid(pathogen~sex,  scales = "free") + 
  theme(legend.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")

ggplot(combine_inc_r %>% filter(year == 2019), aes(x = lowage, y = total_r_cases, group = interaction(combo,pathogen))) + 
  geom_point(aes(col = factor(name))) + 
  geom_smooth(aes(fill = factor(name),col = factor(name)), alpha = 0.2) + 
  facet_grid(pathogen~sex,  scales = "free") + 
  #scale_fill_manual("Bacteria", values = brewer_col) + 
  #scale_color_manual("Bacteria",values = brewer_col) + 
  theme(legend.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Total number of R infections")

## Like Figure 1: 
combine_inc_r <- combine_inc_r %>% mutate(drug = recode(name, "amika_R" = "amikacin","aminogl_R" = "aminoglycoside","carbapen_R" = "carbapenem",
                                                        "fq_pseudo_R" = "fluoroquinolone", "fq_staaur_R" = "fluoroquinolone", "fq_strpne_R" = "fluoroquinolone","fq_ent_R" = "fluoroquinolone",
                                                        "aminopen_R" = "aminopenicillins","genta_high" = "gentamicin","vanco_R" = "vancomycin",
                                                        "cefIII_entero_R" = "3G cephalosporins", "cefIII_strpne_R" = "3G cephalosporins",
                                                        "ert_R" = "ertapenem","ureidopen_R" = "ureidopenicillin", "ceftaz_R" = "ceftazidime",
                                                        "mrsa_R" = "methicillin","rifamp_R"="rifampicin","macrol_R" = "macrolides","penic_RI" = "penicillins"))

g1a <- ggplot(combine_inc_r %>% filter(year == 2019) %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=lowage, y = total_r_cases, group = sex)) + 
  geom_point(aes(colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Number of resistant infections") + 
  scale_color_manual("Sex", values = c("#E41A1C", "#377EB8")) + 
  scale_fill_manual("Sex", values = c("#E41A1C", "#377EB8")) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ drug, scales = "free")
  

g2a <- ggplot(combine_inc_r %>% filter(year == 2019) %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=lowage, y = total_r_cases, group = sex)) + 
  geom_point(aes( colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Number of resistant infections") + 
  scale_color_manual("Sex", values = c("#E41A1C", "#377EB8")) + 
  scale_fill_manual("Sex", values = c("#E41A1C", "#377EB8")) +
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ drug, scales = "free")

g3a <- ggplot(combine_inc_r %>% filter(year == 2019) %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=lowage, y = total_r_cases, group = sex)) + 
  geom_point(aes(colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Number of resistant infections") + 
  scale_color_manual("Sex", values = c("#E41A1C", "#377EB8")) + 
  scale_fill_manual("Sex", values = c("#E41A1C", "#377EB8")) +
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ drug, scales = "free")

g4a <- ggplot(combine_inc_r %>% filter(year == 2019) %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=lowage, y = total_r_cases, group = sex)) + 
  geom_point(aes(colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Number of resistant infections") + 
  scale_color_manual("Sex", values = c("#E41A1C", "#377EB8")) + 
  scale_fill_manual("Sex", values = c("#E41A1C", "#377EB8")) +
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ drug, scales = "free") 

g1a + g2a + g3a + g4a + plot_layout(guides = "collect") 
ggsave("plots/total_r_cases_SUPPFIG.png", width = 24, height = 15)




########################################################***** INPATIENT / OUTPATIENT *********#################################################################################
#########****************** Explore patient type variable                             ******************###################################
##################################################################################################################################################################

######### Combine - PATIENT TYPE with country, age and sex split
data_popsize <- data_popsize %>% rename(country = Country.Code.EARS)
combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(patienttype), !patienttype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, country, sex, pathogen, patienttype, Country.Name, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  mutate(incidence = (total / n_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$patienttype <- factor(combine$patienttype, levels = c("INPAT","OUTPAT","O"))

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(Country.Name~sex + patienttype, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_patienttype_country_lines.pdf", height = 15, width = 10)

# how many per country? 
# Sum over age groups
manypercn <- combine %>% group_by(year, country, sex, patienttype, Country.Name) %>% 
  summarise(totals = sum(total)) %>% 
  group_by(year, sex, Country.Name) %>% mutate(perc = totals / sum(totals))
inpat_vals <- manypercn %>% filter(patienttype == "INPAT") %>% rename(perc_inpat = perc)
outpat_vals <- manypercn %>% filter(patienttype == "OUTPAT") %>% rename(perc_outpat = perc)
manypercn <- manypercn %>% left_join(inpat_vals) %>% fill(perc_inpat) %>% left_join(outpat_vals) %>% group_by(year, country, sex) %>%
  arrange((totals), .by_group = TRUE) %>% fill(perc_outpat)
ggplot(manypercn %>% filter(year == 2019), 
       aes(x=fct_reorder(Country.Name, perc_inpat, .desc = FALSE), y = totals, group = patienttype)) + 
  geom_bar(stat = "identity",position = "fill", aes(fill = patienttype)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_discrete("Country, sex") + 
  scale_y_continuous("Total number of isolates") + 
  scale_fill_discrete("Patient type") + 
  facet_wrap(~sex)
ggsave("plots/inpatient_dist_by_country_SUPPFIG.pdf")

# which countries don't have this data? 
all_countries_with <- data_infs_inflate %>% filter(!patienttype=="UNK") %>% ungroup() %>% distinct(country) %>% unlist()
all_countries <- data_infs_inflate %>% ungroup() %>% distinct(country) %>% unlist()
setdiff(all_countries,all_countries_with) # only UK excluded 


## Remove country split but keep patient type
combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$patienttype <- factor(combine$patienttype, levels = c("INPAT","OUTPAT","O"))

combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(patienttype), !patienttype == "UNK"),
                    data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, sex, pathogen, patienttype, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  group_by(age_group, year, sex, pathogen, patienttype, lowage) %>% 
  summarise(totals = sum(total)) %>% # Can't do this: doesn't give same denominator for each grouping: nage_grp_overcountries = sum(n_inage_gdr)) %>% 
  left_join(., data_popsize_eu_as) %>% # Add in denominator = number of people of that age and sex in Europe
  mutate(incidence = (totals / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))


ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(patienttype ~ sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_patienttype_totals_lines_SUPPFIG.pdf", width = 8, height = 8)

########################################################***** HOSPITAL UNIT TYPE *********#################################################################################
#########****************** Explore hospital unit type variable                             ******************###################################
##################################################################################################################################################################

######### Combine - HOSPITAL UNIT TYPE with country, age and sex split
round(100 * table(data_orig$hospitalunittype) / length(data_orig$hospitalunittype),2)

data_orig %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% filter(reportingcountry == "UK") %>% mutate(percsamp = 100 * nsamp / sum(nsamp))
data_orig %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% filter(reportingcountry == "DE") %>% mutate(percsamp = 100 * nsamp / sum(nsamp))


data_orig %>% filter(year == 2019) %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% 
  ggplot(aes(x=reportingcountry, y = nsamp, group = hospitalunittype)) + geom_bar(position = "dodge", stat = "identity", aes(fill = hospitalunittype)) + scale_x_discrete("Country") + scale_y_continuous("Number of isolates") + scale_fill_discrete("Hospital unit type") 

data_orig %>% filter(year == 2019) %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% 
  ggplot(aes(x=reportingcountry, y = nsamp, group = hospitalunittype)) + geom_bar(position = "stack", stat = "identity", aes(fill = hospitalunittype)) + scale_x_discrete("Country") + scale_y_continuous("Number of isolates") + scale_fill_discrete("Hospital unit type") 

data_orig %>% filter(year == 2019) %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% 
  ggplot(aes(x=reportingcountry, y = nsamp, group = hospitalunittype)) + geom_bar(position = "stack", stat = "identity", aes(fill = hospitalunittype)) + scale_x_discrete("Country") + scale_y_continuous("Number of isolates") + scale_fill_discrete("Hospital unit type") 
ggsave("plots/hospitalunittype_distribution_SUPPFIGURE.pdf", width = 20, height = 10)

data_orig %>% filter(year == 2019) %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% 
  ggplot(aes(x=reportingcountry, y = nsamp, group = hospitalunittype)) + geom_bar(position = "fill", stat = "identity", aes(fill = hospitalunittype)) + scale_x_discrete("Country") + scale_y_continuous("Proportion of isolates") + scale_fill_discrete("Hospital unit type") 
ggsave("plots/hospitalunittype_distribution_prop_SUPPFIGURE.pdf", width = 20, height = 10)


combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(hospitalunittype), !hospitalunittype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, country, sex, pathogen, hospitalunittype, Country.Name, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  mutate(incidence = (total / n_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$hospitalunittype <- factor(combine$hospitalunittype, levels = c("INTMED", "UNK", "ED",  "ICU", "URO", "SURG", "ONCOL", "OBGYN", "O", 
                                                                        "INFECT", "PEDS", "PEDSICU", NA))

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(Country.Name~sex + hospitalunittype, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_hospitalunittype_country_lines_SUPPFIG.pdf", height = 15, width = 20)

### Add incidences? think ok to do as all standardised? 
eleven_cols = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence)) + 
  geom_bar(position = "fill", stat = "identity", aes(fill = factor(hospitalunittype))) + 
  facet_grid(sex ~ country, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of total incidence (infections / 100,000 population per year)") + 
  scale_fill_manual("Hospital unit type", values = eleven_cols) + 
  geom_hline(yintercept = c(0.10, 0.30, 0.5, 0.70), lty = "dashed", alpha = 0.5)
ggsave("plots/incidence_by_hospitalunittype_totals_cntry_bars.pdf", width = 20, height = 8)




## Remove country split but keep hospital type
combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(hospitalunittype), !hospitalunittype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, sex, pathogen, n_inage_gdr,hospitalunittype, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  group_by(age_group, year, sex, pathogen, hospitalunittype, lowage) %>% 
  summarise(totals = sum(total)) %>% #, nage_grp_overcountries = sum(n_inage_gdr)) %>% 
  left_join(., data_popsize_eu_as) %>% # Add in denominator = number of people of that age and sex in Europe
  mutate(incidence = (totals / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine %>% filter(year == 2019, lowage == 20, sex == "f") # checking denominator ok 

combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$hospitalunittype <- factor(combine$hospitalunittype, levels = c("INTMED", "UNK", "ED",  "ICU", "URO", "SURG", "ONCOL", "OBGYN", "O", 
                                                                        "INFECT", "PEDS", "PEDSICU", NA))



ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence, group = interaction(hospitalunittype, pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  facet_grid(hospitalunittype ~ sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_hospitalunittype_totals_lines.pdf", width = 8, height = 8)


### Add incidences? think ok to do as all standardised? 
eleven_cols = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence)) + 
  geom_bar(position = "fill", stat = "identity", aes(fill = factor(hospitalunittype))) + 
  facet_grid( ~ sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_fill_manual("Hospital unit type", values = eleven_cols) + 
  geom_hline(yintercept = c(0.10, 0.30, 0.5, 0.70), lty = "dashed", alpha = 0.5)
ggsave("plots/incidence_by_hospitalunittype_totals_bars.pdf", width = 12, height = 8)
















########################################################***** EXPLORE *********#################################################################################
#########****************** Exploration of other combinations / summaries of incidence as checks / data exploration ******************###################################
##################################################################################################################################################################

# sex split over time
ggplot(combine_eu , aes(x = lowage, y = incidence, group = year)) + 
  geom_point(aes(col = factor(year))) + 
  geom_smooth(aes(fill = factor(year),col = factor(year)), alpha = 0.2) + 
  facet_grid(pathogen~sex, scales = "free") + 
  scale_fill_discrete("Year") + 
  scale_color_discrete("Year") + 
  theme(strip.text.y = element_text(angle = 0,face="italic")) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")
ggsave("plots/incidence_sex_sep_combine_eu.pdf")

# Bug level sex split
g2all<- ggplot(combine_eu %>% filter(year == 2019), aes(x = lowage, y = incidence, group = sex)) + 
  geom_point(aes(col = factor(sex))) + 
  geom_smooth(aes(fill = factor(sex),col = factor(sex)), alpha = 0.2) + 
  facet_wrap(~pathogen, scales = "free", ncol = 2) + 
  scale_fill_discrete("Sex") + 
  scale_color_discrete("Sex") + 
  theme(strip.text = element_text(angle = 0,face="italic"), legend.position = "bottom") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")
ggsave("plots/incidence_sex_sep_combine_eu_all.pdf")


#########*## Combine - with no BACTERIA split to compare to Europe wide BSI population estimates in literature ****#############################################
data_popsize_eu_all <- data_popsize %>% group_by(year) %>% summarise(total_inage_gdr = sum(n_inage_gdr))
data_infs_inflate_eu_all <- data_infs_inflate %>% filter(age > 0) %>% group_by(year) %>% summarise(total_inf = sum(n_inflat))
combine_eu_all = left_join(data_infs_inflate_eu_all,data_popsize_eu_all, by = c("year")) %>% 
  mutate(incidence = (total_inf / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

ggplot(combine_eu_all, aes(x = year, y = incidence)) + 
  geom_bar(stat = "identity",position = "dodge") +
  scale_x_continuous("Year") + 
  scale_y_continuous("Incidence (infections per 100,000 people per year)") + 
  ggtitle("Europe")
ggsave("plots/incidence_total_europe.pdf")


#########*## Combine - with just country split to compare to country BSI population estimates in literature ****#############################################
data_popsize_country_only <- data_popsize %>% group_by(year, country) %>% summarise(total_inage_gdr = sum(n_inage_gdr))
data_infs_inflate_country_only <- data_infs_inflate %>% filter(age > 0) %>% group_by(country,year) %>% summarise(total_inf = sum(n_inflat))
combine_country_only = left_join(data_infs_inflate_country_only,data_popsize_country_only, by = c("year","country")) %>% 
  mutate(incidence = (total_inf / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

ggplot(combine_country_only, aes(x = year, y = incidence)) + 
  geom_bar(stat = "identity",position = "dodge", aes()) + 
  facet_wrap(~country)
write.csv(combine_country_only, "output/incidence_country_level.csv")

g1 <- ggplot(combine_country_only, aes(x = year, y = incidence, group = country)) + 
  geom_point(aes(col = country)) + geom_line(aes(col = country)) + 
  scale_x_continuous("Year") + 
  scale_y_continuous("Incidence (infections per 100,000 people per year)") + 
  scale_color_discrete("Country")

ggplot(combine_country_only, aes(x = year, y = incidence)) + 
  geom_point(aes(col = country)) + 
  geom_smooth() + 
  scale_x_continuous("Year") + 
  scale_y_continuous("Incidence (infections per 100,000 people per year)") + 
  scale_color_discrete("Country")

g2 <- ggplot(combine_country_only %>% filter(!year == 2020), aes(x = year, y = incidence)) + 
  geom_point(aes(col = country)) + 
  geom_smooth() + 
  geom_point(data = combine_country_only) + 
  scale_x_continuous("Year") + 
  scale_y_continuous("Incidence (infections per 100,000 people per year)") + 
  scale_color_discrete("Country") + theme(legend.position = "none") 

g1 + g2 
ggsave("plots/incidence_country_trend_SUPPFIG.pdf",width = 10, height = 5)
 

combine_country_only %>% filter(country == 'FI') # Goto review => 168 in 2007 => similar ball park
combine_country_only %>% filter(country == 'DK') # Goto review => 166 in 2006 => similar ball park but increasing


#########*## Combine - with bug + country split to compare to bacteria and country population estimates in literature ****#############################################
data_popsize_country <- data_popsize %>% group_by(year, country) %>% summarise(total_inage_gdr = sum(n_inage_gdr))
data_infs_inflate_country <- data_infs_inflate %>% filter(age > 0) %>% group_by(country,year,pathogen) %>% summarise(total_inf = sum(n_inflat))
combine_country = left_join(data_infs_inflate_country,data_popsize_country, by = c("year","country")) %>% 
  mutate(incidence = (total_inf / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

ggplot(combine_country, aes(x = year, y = incidence)) + 
  geom_bar(stat = "identity",position = "dodge", aes(fill = factor(pathogen))) + 
  facet_grid(pathogen~country, scales = "free") + 
  scale_x_continuous("Year") + 
  scale_y_continuous("Incidence (infections per 100,000 people per year)") + 
  scale_fill_discrete("Bacteria") +
  theme(legend.text = element_text(angle = 0,face="italic")) + 
  theme(strip.text.y = element_text(face="italic")) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("plots/incidence_country_bacteria_trend_SUPPFIG.pdf",width = 20, height = 20)

combine_country %>% filter(country == 'UK', pathogen == "Escherichia coli") # MacKinnon review => 60.4 in 2012-13 for England for Esccol => exactly the same number but for UK in 2015
combine_country %>% filter(country == 'DK', pathogen == "Escherichia coli") # MacKinnon review => 56-70 <2010 for region in Denmark for Esccol => similar ball park but higher now
combine_country %>% filter(country == 'SE', pathogen == "Escherichia coli") # MacKinnon review => 67 in 2011 for region in Sweden for Esccol => similar ball park but higher now
combine_country %>% filter(country == 'FI', pathogen == "Escherichia coli") # MacKinnon review => 67 in 2011 for region in Sweden for Esccol => similar ball park but higher now

combine_country %>% filter(country == 'IS', pathogen == "Staphylococcus aureus") # Laupland 22-29 up to 2008 for Iceland => same ball park 
combine_country %>% filter(country == 'DK', pathogen == "Staphylococcus aureus") # Laupland 26 up to 2008 for Denmark => higher now...
combine_country %>% filter(country == 'FI', pathogen == "Staphylococcus aureus") # Laupland 26 up to 2008 for Finland => higher now... 

# For appendix 
write_csv(combine_country,"output/country_pathogen_incidence_estimates.csv")


############ Combine - with just age and sex split
data_popsize_eu_as <- data_popsize %>% group_by(year, sex, lowage, age_group) %>% summarise(total_inage_gdr = sum(n_inage_gdr))
data_infs_inflate_eu_as <- data_infs_inflate %>% filter(age > 0) %>% group_by(year,sex,age_group) %>% summarise(total_inf = sum(n_inflat))
combine_eu_as = left_join(data_infs_inflate_eu_as,data_popsize_eu_as, by = c("year","sex","age_group")) %>% 
  mutate(incidence = (total_inf / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine_eu_as[which(combine_eu_as$sex == "f"),"sex"] <- "Female" 
combine_eu_as[which(combine_eu_as$sex == "m"),"sex"] <- "Male" 
ggplot(combine_eu_as , aes(x = lowage, y = incidence, group = year)) + 
  geom_point(aes(col = factor(year))) + 
  geom_smooth(aes(fill = factor(year),col = factor(year)), alpha = 0.2) + 
  facet_grid(~sex, scales = "free") + 
  scale_fill_discrete("Year") + 
  scale_color_discrete("Year") + 
  theme(strip.text.y = element_text(angle = 0,face="italic")) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")
ggsave("plots/incidence_sex_years_combine_eu.pdf")

ggplot(combine_eu_as %>% filter(year == 2019), aes(x = lowage, y = incidence, group = sex)) + 
  geom_point(aes(col = factor(sex))) + 
  geom_smooth(aes(fill = factor(sex),col = factor(sex)), alpha = 0.2) + 
  scale_fill_discrete("Sex") + 
  scale_color_discrete("Sex") + 
  theme(strip.text.y = element_text(angle = 0,face="italic")) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year) for 2019")
ggsave("plots/incidence_sex_2019_combine_eu.pdf")


g3 <- ggplot(combine_eu , aes(x = lowage, y = incidence, group = interaction(sex,year))) + 
  geom_point(aes(col = factor(year))) + 
  geom_smooth(aes(fill = factor(year),col = factor(year), linetype = sex), alpha = 0.2) + 
  facet_wrap(~pathogen, scales = "free") + 
  scale_fill_discrete("Year") + 
  scale_color_discrete("Year") + 
  theme(strip.text.x = element_text(angle = 0,face="italic")) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_linetype_discrete("Sex")
ggsave("plots/incidence_sex_tog_combine_eu_SUPPFIG.pdf")

ggplot(combine_eu, aes(x = year, y = incidence, group = lowage)) + 
  geom_point(aes(col = factor(lowage))) + 
  geom_smooth(aes(fill = factor(lowage),col = factor(lowage)), alpha = 0.2) + 
  facet_grid(pathogen~sex, scales = "free") + 
  scale_fill_discrete("Age") + 
  scale_color_discrete("Age") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Year") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)")
ggsave("plots/incidence_year_sex_sep_combine_eu.pdf")

combine_eu %>% filter(sex == "m", lowage == 15, pathogen == "esccol")

######### Combine - with country, age and sex split - used above for R calcs too but 
combine <- read_csv("output/country_pathogen_age_incidence_estimates.csv")
combine <- rename(combine, sex = gender)

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(sex~Country.Name, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0))
ggsave("plots/incidence_by_country_lines.pdf")

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_smooth(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(Country.Name~sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0))
ggsave("plots/incidence_by_country_smooth.pdf")

## Is men > women in all countries? yes
ggplot(combine %>% filter(year == 2019), 
       aes(x=lowage, y = incidence, group = interaction(sex,Country.Name,year)), col = sex) + 
  geom_line(aes(col = sex)) + 
  facet_grid(Country.Name~pathogen, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0))

## women > men for E coli in 25-45? yes
ggplot(combine %>% filter(year == 2019, lowage < 50, pathogen == "Escherichia coli"), 
       aes(x=lowage, y = incidence, group = interaction(sex,Country.Name,year)), col = sex) + 
  geom_line(aes(col = sex)) + 
  facet_grid(Country.Name~pathogen, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0))


## Does the value vary by year? 
ggplot(combine %>% filter(sex == "m"), aes(x=lowage, y = incidence, group = interaction(pathogen,year, Country.Name))) + 
  geom_smooth(aes(col = factor(pathogen)),alpha = 0.2) + 
  facet_wrap(~Country.Name, scales = "free") + 
  ggtitle("Men")

ggplot(combine %>% filter(sex == "f"), aes(x=lowage, y = incidence, group = interaction(pathogen,year, Country.Name))) + 
  geom_smooth(aes(col = factor(pathogen)),alpha = 0.2) + 
  facet_wrap(~Country.Name, scales = "free") + 
  ggtitle("Women")

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_smooth(aes(col = factor(year)),alpha = 0.05) + 
  facet_grid(Country.Name~sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0))

ggplot(combine %>% filter(Country.Name == "Finland"), aes(x=lowage, y = incidence, group = interaction(pathogen,year, Country.Name))) + 
  geom_smooth(aes(col = factor(pathogen)),alpha = 0.2) + 
  facet_grid(~sex, scales = "free")

ggplot(combine, aes(x=year, y = incidence, group = interaction(pathogen,lowage))) + 
  geom_line(aes(col = factor(lowage))) + 
  facet_grid(pathogen~sex + country, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0))




######### Combine - PATIENT TYPE with country, age and sex split
combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(patienttype), !patienttype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, country, sex, pathogen, patienttype, Country.Name, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  mutate(incidence = (total / n_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$patienttype <- factor(combine$patienttype, levels = c("INPAT","OUTPAT","O"))

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(Country.Name~sex + patienttype, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_patienttype_country_lines.pdf", height = 15, width = 10)

# how many per country? 
manypercn <- combine %>% group_by(year, country, sex, patienttype, Country.Name) %>% 
  summarise(totals = sum(total)) %>% 
  group_by(year, sex, Country.Name) %>% mutate(perc = totals / sum(totals))
inpat_vals <- manypercn %>% filter(patienttype == "INPAT") %>% rename(perc_inpat = perc)
outpat_vals <- manypercn %>% filter(patienttype == "OUTPAT") %>% rename(perc_outpat = perc)
manypercn <- manypercn %>% left_join(inpat_vals) %>% fill(perc_inpat) %>% left_join(outpat_vals) %>% group_by(year, country, sex) %>%
  arrange((totals), .by_group = TRUE) %>% fill(perc_outpat)
ggplot(manypercn %>% filter(year == 2019), 
       aes(x=fct_reorder(Country.Name, perc_inpat, .desc = FALSE), y = totals, group = patienttype)) + 
  geom_bar(stat = "identity",position = "fill", aes(fill = patienttype)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_discrete("Country, sex") + 
  scale_y_continuous("Total number of isolates") + 
  scale_fill_discrete("Patient type") + 
  facet_wrap(~sex)
ggsave("plots/inpatient_dist_by_country.pdf")

# which countries don't have this data? 
all_countries_with <- data_infs_inflate %>% filter(!patienttype=="UNK") %>% ungroup() %>% distinct(country) %>% unlist()
all_countries <- data_infs_inflate %>% ungroup() %>% distinct(country) %>% unlist()
setdiff(all_countries,all_countries_with) # only UK excluded 


## Remove country split but keep patient type
combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$patienttype <- factor(combine$patienttype, levels = c("INPAT","OUTPAT","O"))

combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(patienttype), !patienttype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, sex, pathogen, patienttype, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  group_by(age_group, year, sex, pathogen, patienttype, lowage) %>% 
  summarise(totals = sum(total)) %>% # Can't do this: doesn't give same denominator for each grouping: nage_grp_overcountries = sum(n_inage_gdr)) %>% 
  left_join(., data_popsize_eu_as) %>% # Add in denominator = number of people of that age and sex in Europe
  mutate(incidence = (totals / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))


ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(patienttype ~ sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_patienttype_totals_lines.pdf", width = 8, height = 8)

######### Combine - HOSPITAL UNIT TYPE with country, age and sex split
round(100 * table(data_orig$hospitalunittype) / length(data_orig$hospitalunittype),2)

data_orig %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% 
  ggplot(aes(x=reportingcountry, y = nsamp, group = hospitalunittype)) + geom_bar(position = "dodge", stat = "identity", aes(fill = hospitalunittype)) + scale_x_discrete("Country") + scale_y_continuous("Number of isolates") + scale_fill_discrete("Hospital unit type") 
ggsave("plots/hospitalunittype_distribution.pdf")

data_orig %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% filter(reportingcountry == "UK")
data_orig %>% group_by(reportingcountry, hospitalunittype) %>% select(reportingcountry, hospitalunittype) %>% summarise(nsamp = n()) %>% filter(reportingcountry == "DE")

combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(hospitalunittype), !hospitalunittype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, country, sex, pathogen, hospitalunittype, Country.Name, n_inage_gdr, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  mutate(incidence = (total / n_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$hospitalunittype <- factor(combine$hospitalunittype, levels = c("INTMED", "UNK", "ED",  "ICU", "URO", "SURG", "ONCOL", "OBGYN", "O", 
                                                                        "INFECT", "PEDS", "PEDSICU", NA))

ggplot(combine, aes(x=lowage, y = incidence, group = interaction(pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  geom_point(aes(col = factor(pathogen)),alpha = 0.05) + 
  facet_grid(Country.Name~sex + hospitalunittype, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_hospitalunittype_country_lines.pdf", height = 15, width = 20)

### Add incidences? think ok to do as all standardised? 
eleven_cols = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence)) + 
  geom_bar(position = "fill", stat = "identity", aes(fill = factor(hospitalunittype))) + 
  facet_grid(sex ~ country, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of total incidence (infections / 100,000 population per year)") + 
  scale_fill_manual("Hospital unit type", values = eleven_cols) + 
  geom_hline(yintercept = c(0.10, 0.30, 0.5, 0.70), lty = "dashed", alpha = 0.5)
ggsave("plots/incidence_by_hospitalunittype_totals_cntry_bars.pdf", width = 20, height = 8)


## Remove country split but keep hospital type
combine = left_join(data_infs_inflate %>% filter(age > 0, !is.na(hospitalunittype), !hospitalunittype == "UNK"),data_popsize, by = c("year","country","sex","age_group")) %>% 
  group_by(age_group, year, sex, pathogen, n_inage_gdr,hospitalunittype, lowage) %>% 
  summarise(total = sum(n_inflat)) %>% 
  group_by(age_group, year, sex, pathogen, hospitalunittype, lowage) %>% 
  summarise(totals = sum(total)) %>% #, nage_grp_overcountries = sum(n_inage_gdr)) %>% 
  left_join(., data_popsize_eu_as) %>% # Add in denominator = number of people of that age and sex in Europe
  mutate(incidence = (totals / total_inage_gdr) * 100000) %>% 
  filter(!is.na(incidence))

combine %>% filter(year == 2019, lowage == 20, sex == "f") # checking denominator ok 

combine[which(combine$sex == "f"),"sex"] <- "Female" 
combine[which(combine$sex == "m"),"sex"] <- "Male" 
combine$hospitalunittype <- factor(combine$hospitalunittype, levels = c("INTMED", "UNK", "ED",  "ICU", "URO", "SURG", "ONCOL", "OBGYN", "O", 
                                                                        "INFECT", "PEDS", "PEDSICU", NA))



ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence, group = interaction(hospitalunittype, pathogen,year))) + 
  geom_line(aes(col = factor(pathogen))) + 
  facet_grid(hospitalunittype ~ sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_color_discrete("Bacteria") + 
  theme(legend.text = element_text(angle = 0,face="italic")) 
ggsave("plots/incidence_by_hospitalunittype_totals_lines.pdf", width = 8, height = 8)


### Add incidences? think ok to do as all standardised? 
eleven_cols = c("#a50026","#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4","#313695")
ggplot(combine %>% filter(year == 2019), aes(x=lowage, y = incidence)) + 
  geom_bar(position = "fill", stat = "identity", aes(fill = factor(hospitalunittype))) + 
  facet_grid( ~ sex, scales = "free") + 
  theme(strip.text.y = element_text(angle = 0)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Incidence (infections / 100,000 population per year)") + 
  scale_fill_manual("Hospital unit type", values = eleven_cols) + 
  geom_hline(yintercept = c(0.10, 0.30, 0.5, 0.70), lty = "dashed", alpha = 0.5)
ggsave("plots/incidence_by_hospitalunittype_totals_bars.pdf", width = 12, height = 8)





