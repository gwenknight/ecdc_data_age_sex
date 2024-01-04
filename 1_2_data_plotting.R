#### Basic data plotting 
## Explore the trends in AMR for the data

########*********  Setup  ******* #######
## Libraries and plot setup
library(tidyverse)
library(gganimate)
library(gifski)
library(ggformula) # for geom_spline
library(patchwork)
theme_set(theme_bw(base_size = 11))
colours_to_use <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

## Data
data <- read_csv("data/data_cleaned.csv")

######*********** Summary stats ******* #######
# Number of susceptibility tests
data %>% ungroup() %>% summarise(total = sus + res) %>% summarise(sum(total))

# By drug / bug
data %>% ungroup() %>% group_by(pathogen, name) %>% summarise(total = sus + res) %>%
  group_by(pathogen, name) %>% summarise(totals = sum(total))

# Demographics
range(data$age)
data %>% group_by(gender) %>% summarise(mean(age), range(age), median(age), sd(age))

d <- data %>% group_by(gender) %>% summarise(total = sus+res) %>% summarise(sum(total))
d[1,2] / (d[1,2] + d[2,2])*100 # 47% are from females


########*********  Usual plots  ******* #######
#(1) Resistance proportions over time across all of Europe 
europe <- data %>% ungroup() %>% group_by(year, pathogen, name) %>% dplyr::summarise(totals = sum(sus), totalr = sum(res)) %>% 
  mutate(total= totals+totalr, proportion = totalr / total)

ggplot(europe, aes(x=year, y = proportion, size = total, colour = name)) + geom_point(alpha = 0.7, aes()) + 
  facet_wrap(~pathogen, ncol = 4) + 
  scale_size("Number of isolates") + 
  scale_color_discrete("Type of resistance")
ggsave("plots/basic_europe.pdf")

#(2) Resistance proportions over time in each country
country_time <- data %>% ungroup() %>% group_by(year, pathogen, name, country) %>% dplyr::summarise(totals = sum(sus), totalr = sum(res)) %>%
  mutate(total= totals+totalr, proportion = totalr / total)

ggplot(country_time, aes(x=year, y = proportion,colour = country)) + geom_line() + 
  facet_wrap(name~pathogen, ncol = 4) + 
  scale_color_discrete("Country")
ggsave("plots/basic_country_time_all.pdf")

ggplot(country_time %>% filter(pathogen %in% c("esccol", "staaur"), name %in% c("mrsa_R", "fq_ent_R")), aes(x=year, y = proportion,colour = country)) + geom_line() + 
  facet_wrap(name~pathogen, ncol = 4) + 
  scale_color_discrete("Country")
ggsave("plots/basic_country_time_egs.pdf")

ggplot(country_time, aes(x=reorder(country, proportion), y = proportion)) + geom_boxplot() + 
  scale_x_discrete("Country") + 
  scale_y_continuous("Proportion of isolates resistant across data")
ggsave("plots/basic_country_time_boxpot.pdf")

ggplot(country_time, aes(x=reorder(country, proportion), y = proportion)) + geom_boxplot() + 
  scale_x_discrete("Country") + 
  facet_wrap(name~pathogen, ncol = 4) + 
  scale_y_continuous("Proportion of isolates resistant across data")
ggsave("plots/basic_country_time_boxpot_all.pdf")

### (3) Animate
p1 <- ggplot(country_time, aes(x=country, y = proportion, size = totals)) + geom_point(alpha=0.6, aes(colour = country), show.legend = FALSE) +  
  facet_wrap(~pathogen + name) +
  theme(axis.text.x = element_text(angle = 90)) + 
  labs(title = 'Year: {round(frame_time,0)}', y = "Proportion resistant", x = "Country") +
  transition_time(year) +
  ease_aes('linear')

animate(p1)
anim_save("plots/animate_per_cntry_yer.gif", width = 15, height = 10)

### Exploring Sex and age => need to account for population size underlying
# For incidence
data_popsize <- read_csv("data/data_popsize_clean.csv")[,-1]
data_popsize_eu <- data_popsize %>% filter(year == 2019) %>% group_by(year, sex, lowage) %>% summarise(total_inage_gdr = sum(n_inage_gdr)) %>% 
  mutate(age_group = cut(lowage, breaks = seq(-5,115,5))) 

# MRSA 
data <- data %>% mutate(sex = gender)
data_mrsa <- data %>% filter(pathogen == "staaur", year == 2019, name == "mrsa_R") %>% mutate(age_group = cut(age, breaks = seq(-5,130,5))) %>% 
  group_by(age_group, sex) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(age_group, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t))

data_inci <- left_join(data_mrsa, data_popsize_eu) %>% mutate(incidence_infection = 100000 * total / total_inage_gdr, incidence_R_infection = 100000 * res_t / total_inage_gdr, )
g1 <- ggplot(data_inci, aes(x=age_group, y = proportion, group = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Proportion of isolates tested\n that are resistant") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
g2 <- ggplot(data_inci, aes(x=age_group, y = total, group = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Total number of isolates \n ~ number of BSI") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
g3 <- ggplot(data_inci, aes(x=age_group, y = res_t, group = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Total number of resistant isolates") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
g4 <- ggplot(data_inci, aes(x=age_group, y = incidence_infection, group = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Incidence of infection per 100,000 inhabitants") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
g5 <- ggplot(data_inci, aes(x=age_group, y = incidence_R_infection, group = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Incidence of R infection per 100,000 inhabitants") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))

g1 + g2 + g4 + guide_area() + g3 + g5 + plot_layout(guides = "collect") + plot_annotation(title = "MRSA")

# What does this say re: sex?
# says that despite women getting more antibiotics, resistance proportion same - but incidence much higher in men of infections with R bacteria

# run through all combinations 
dir.create(file.path("plots", "age_sex_summaries"), showWarnings = FALSE)
setwd(file.path("plots", "age_sex_summaries"))
bugs <- unique(data$pathogen)
data <- data %>% filter(age > 0) %>% mutate(age_group = cut(age, breaks = seq(0,120,5)))

for(i in 1:length(bugs)){
  data_bug <- data %>% filter(pathogen == bugs[i])
  drugs <- unique(data_bug$name)
  for(j in 1:length(drugs)){
    print(paste0(bugs[i], ", ", drugs[j]))
    data_drug_bug <- data_bug %>% filter(year == 2019, name == drugs[j]) %>% 
      group_by(age_group, sex) %>% 
      summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
      group_by(age_group, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t))
    
    data_inci <- left_join(data_drug_bug, data_popsize_eu) %>% mutate(incidence_infection = 100000 * total / total_inage_gdr, incidence_R_infection = 100000 * res_t / total_inage_gdr, )
    g1 <- ggplot(data_inci, aes(x=age_group, y = proportion, group = sex)) + geom_point(aes(col = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Proportion of isolates tested\n that are resistant") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
    g2 <- ggplot(data_inci, aes(x=age_group, y = total, group = sex))  + geom_point(aes(col = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Total number of isolates \n ~ number of BSI") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
    g3 <- ggplot(data_inci, aes(x=age_group, y = res_t, group = sex))  + geom_point(aes(col = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Total number of resistant isolates") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
    g4 <- ggplot(data_inci, aes(x=age_group, y = incidence_infection, group = sex))  + geom_point(aes(col = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Incidence of infection per 100,000 inhabitants") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
    g5 <- ggplot(data_inci, aes(x=age_group, y = incidence_R_infection, group = sex))  + geom_point(aes(col = sex)) + geom_smooth(aes(col = sex, fill = sex)) + scale_y_continuous("Incidence of R infection per 100,000 inhabitants") + scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + theme(axis.text.x = element_text(angle = 90))
    
    g <- g1 + g2 + g4 + guide_area() + g3 + g5 + plot_layout(guides = "collect") + plot_annotation(title = paste0(bugs[i], ", ", drugs[j]))
    ggsave(paste0(bugs[i], "_", drugs[j],"summary.pdf"))
  }
}
setwd("../..")
# What does this say re: sex?
# Proportion R in women >> men: encfae + encfai vanc_R + escol, carb_R at ages < 50
# No bug where incidence in women >> men: all men higher or similar 

##### All in one 
data_long <- data %>% filter(age > 0) %>% 
  mutate(pathogen = recode(pathogen, "encfae" = "Enterococcus faecalis", "encfai" = "Enterococcus faecium", 
                           "esccol" = "Escherichia coli", "klepne" = "Klebsiella pneumoniae", 
                           "staaur" = "Staphylococcus aureus", "strpne" = "Streptococcus pneumoniae", 
                           "acispp" = "Acinetobacter spp", "pseaer"="Pseudomonas aeruginosa"),
         drug = recode(name, "amika_R" = "amikacin","aminogl_R" = "aminoglycoside","carbapen_R" = "carbapenem",
                       "fq_pseudo_R" = "fluoroquinolone", "fq_staaur_R" = "fluoroquinolone", "fq_strpne_R" = "fluoroquinolone","fq_ent_R" = "fluoroquinolone",
                       "aminopen_R" = "aminopenicillins","genta_high" = "gentamicin","vanco_R" = "vancomycin",
                       "cefIII_entero_R" = "3G cephalosporins", "cefIII_strpne_R" = "3G cephalosporins",
                       "ert_R" = "ertapenem","ureidopen_R" = "ureidopenicillin", "ceftaz_R" = "ceftazidime",
                       "mrsa_R" = "methicillin","rifamp_R"="rifampicin","macrol_R" = "macrolides","penic_RI" = "penicillins"),
         family = ifelse(drug %in% c("penicillins","aminopenicillins","ureidopenicillin","methicillin","carbapenem","ertapenem","ceftazidime","3G cephalosporins"),"beta-lactam",
                         ifelse(drug %in% c("aminoglycoside","amikacin","gentamicin"),"aminoglycoside",
                                ifelse(drug == "fluoroquinolone","fluoroquinolone",
                                       ifelse(drug == "vancomycin","glycopeptide",
                                              ifelse(drug == "macrolides","macrolides","rifamycin"))))),
         aware = ifelse(drug %in% c("penicillins", "aminoglycoside"), "Access + Watch",
                        ifelse(drug %in% c("3G cephalosporins"), "Watch + Reserve",
                               ifelse(drug %in% c("amikacin","aminopenicillins","gentamicin", "methicillin"), "Access","Watch")))) # check out table S3 in appendix


data_bybugdrug <- data_long %>% filter(year < 2020, year > 2014) %>% 
  mutate(age_group = cut(age, breaks = seq(-5,130,5))) %>% 
  group_by(pathogen, drug, age_group, sex) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, drug, age_group, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10)


g1a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=age_group, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ drug)

g2a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age_group, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ drug)

g3a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age_group, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ drug)

g4a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age_group, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ drug) 

g1a + g2a + g3a + g4a + plot_layout(guides = "collect") + ggtitle("European level (2015-2019)")
ggsave("plots/European_level_byage_bug_drug_5yr.png", width = 15, height = 10)

##########################################################################################################################################################
## individual age
##########################################################################################################################################################
data_bybugdrug <- data_long %>% filter(year < 2020, year > 2014) %>% 
  group_by(pathogen, drug, family, age, sex) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, drug, family, age, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10)


g1a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_size(breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ family + drug)

g2a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_size(breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ family + drug)

g3a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_size(breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  facet_grid(pathogen ~ family + drug)

g4a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_size(breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug) 

(g1a + g2a + g3a + g4a + plot_layout(guides = "collect")) + plot_annotation("European level (2015-2019)")
ggsave("plots/European_level_byage_bug_drug_1yr.png", width = 15, height = 10)

##########################################################################################################################################################
### colour by family
##########################################################################################################################################################
data_bybugdrug$family <- factor(data_bybugdrug$family)
familydrugpal <- c('#c51b7d','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd')
g1a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g2a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g3a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) +  
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g4a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", "rifamycin", "macrolides"), 
                     labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug) 

(g1a + g2a + g3a + g4a + plot_layout(guides = "collect")) + plot_layout(guides = "collect") + plot_annotation("European level (2015-2019)")
ggsave("plots/European_level_byage_bug_drug_1yr_drugfamily.png", width = 15, height = 10)




## individual age
data_bybugdrug <- data_long %>% filter(year < 2020, year > 2014) %>% 
  group_by(pathogen, family, age, sex) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, family, age, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10)

ggplot(data_bybugdrug , aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = sex)) + 
  geom_smooth(aes(col = sex, fill = sex)) + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male")) + 
  scale_size(breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family) 

##########################################################################################################################################################
#### Subregion
##########################################################################################################################################################
## ## Add in subregion
subregion <- read.csv("data/subregions.csv")
#cols_subr <- c('#1a9641','#0066FF','#9933FF','#d7191c','#a6d96a')
cols_subr <- c('#abd9e9','#2c7bb6','#d7191c','#fdae61') # 5 diverging colourblind safe from colorbrewer2: removed Western Asia so only need 4 now (removed yellow in middle): Cyprus in Southern Europe

data_bybugdrugsubr <- data_long %>% filter(year < 2020, year > 2014) %>% 
  left_join(., subregion) %>% 
  group_by(pathogen, drug, family, age, sex, subregion) %>%
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, drug, family, age, sex,subregion) %>% 
  mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10) 
data_bybugdrugsubr$subregion <- factor(data_bybugdrugsubr$subregion, levels = c("Northern Europe", "Western Europe", "Southern Europe", "Eastern Europe"))

g1a <- ggplot(data_bybugdrugsubr %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)+ 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 


g2a <- ggplot(data_bybugdrugsubr %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion", values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion", values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)+ 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 

g3a <- ggplot(data_bybugdrugsubr %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion", values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion", values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)+ 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 

g4a <-  ggplot(data_bybugdrugsubr %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion", values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion", values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 

(g1a + g2a + g3a + g4a + plot_layout(guides = "collect", width = c(1,1.2))) + plot_annotation("Subregion level (2015-2019)")
ggsave("plots/European_level_byage_bug_drug_1yr_drugfamily_SUBREGION.png", width = 22, height = 10)

ggplot(data_bybugdrugsubr %>% filter(pathogen %in% c("Acinetobacter spp")), aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)+ 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 
ggsave("plots/European_level_byage_bug_drug_1yr_drugfamily_SUBREGION_acinetobacter.png", width = 22, height = 10)

data_eg <- rbind(data_bybugdrugsubr %>% filter(pathogen %in% c("Acinetobacter spp")) %>% filter(drug %in% c("carbapenem")),
                 data_bybugdrugsubr %>% filter(pathogen %in% c("Escherichia coli","Staphylococcus aureus")) %>%
                   filter(drug %in% c("methicillin","aminopenicillins")),
                 data_bybugdrugsubr %>% filter(pathogen %in% c("Klebsiella pneumoniae")) %>% filter(drug %in% c("3G cephalosporins")))
g1eg <- ggplot(data_eg %>% filter(pathogen == "Acinetobacter spp"),aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  ggtitle("Acinetobacter spp, carbapenem") + 
  theme(plot.title = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 
g4eg <- ggplot(data_eg %>% filter(pathogen == "Klebsiella pneumoniae"),aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  ggtitle("Klebsiella pneumoniae, 3G cephalosporins") + 
  theme(plot.title = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 
g2eg <- ggplot(data_eg %>% filter(pathogen == "Escherichia coli"),aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  ggtitle("Escherichia coli, aminopenicillins") + 
  theme(plot.title = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 
g3eg <- ggplot(data_eg %>% filter(pathogen == "Staphylococcus aureus"),aes(x=age, y = proportion, group = interaction(subregion, sex))) + 
  geom_point(aes(size = total, alpha = 0.5, colour = subregion)) + 
  geom_smooth(aes(fill = subregion, col = subregion, lty = sex)) + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant", lim = c(0,1)) + 
  scale_color_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_fill_manual("Subregion",  values = cols_subr, drop = FALSE) +
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  ggtitle("Staphylococcus aureus, methicillin") + 
  theme(plot.title = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none") 

g1eg + g4eg + g2eg + g3eg + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'a')
ggsave("plots/European_level_byage_bug_drug_1yr_drugfamily_SUBREGION_egs.png", width = 22, height = 18)

##########################################################################################################################################################
############ Check certain countries don't drive patterns
##########################################################################################################################################################
#### Importance of countries in this data plot 
# top ones are UK / DE / IT / FR 
cntn = data %>% mutate(total = sus + res) %>% 
  group_by(country) %>% summarise(n = sum(total)) %>% arrange(desc(n))
ggplot(cntn, aes(x=reorder(country, -n), y = n)) + geom_bar(stat = "identity")
# UK much higher than others 
### colour by family - use code from above on data without UK 
data_bybugdrug <- data_long %>% filter(!country == "UK", year < 2020, year > 2014) %>% 
  group_by(pathogen, drug, family, age, sex) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, drug, family, age, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10)
data_bybugdrug$family <- factor(data_bybugdrug$family)

g1a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"),type = c("#E41A1C", "#377EB8")) + 
  
  g2a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) +  
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g3a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g4a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", "rifamycin", "macrolides"), 
                     labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug) 

(g1a + g2a + g3a + g4a + plot_layout(guides = "collect")) + plot_layout(guides = "collect") + 
  plot_annotation("European level (2015-2019) without UK")
ggsave("plots/European_level_byage_bug_drug_1yr_drugfamily_withoutUK.png")

#### Remove top four 
data_bybugdrug <- data_long %>% filter(!country %in% c("UK","DE","IT","FR"), year < 2020, year > 2014) %>% 
  group_by(pathogen, drug, family, age, sex) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, drug, family, age, sex) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10)
data_bybugdrug$family <- factor(data_bybugdrug$family)

g1a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g2a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g3a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                               "rifamycin", "macrolides"), labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", 
                                                                                      "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"),type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug)

g4a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age, y = proportion, group = sex)) + 
  geom_point(aes(size = total, colour = family)) + 
  geom_smooth(aes(fill = sex), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family", breaks = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", "rifamycin", "macrolides"), 
                     labels = c("aminoglycoside", "beta-lactam", "fluoroquinolone", "glycopeptide", "rifamycin", "macrolides"),
                     values = familydrugpal, drop = FALSE) + 
  scale_fill_discrete("Sex", breaks = c("f","m"), labels = c("female","male"), type = c("#E41A1C", "#377EB8")) + 
  scale_size("# samples", breaks = seq(2500,20000,2500)) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  facet_grid(pathogen ~ family + drug) 

(g1a + g2a + g3a + g4a + plot_layout(guides = "collect")) + plot_layout(guides = "collect") + 
  plot_annotation("European level (2015-2019) without UK, IT, FR, DE")
ggsave("plots/European_level_byage_bug_drug_1yr_drugfamily_withoutUKITFRDE.png")

##########################################################################################################################################################
###### AWARE grouping
##########################################################################################################################################################
data_long %>% filter(family == "beta-lactam") %>% summarise(unique(aware))
data_long <- data_long %>% 
  mutate(family_aware = ifelse(family == "beta-lactam", paste0(paste(family, aware,sep="\n("),")"), family))
# familydrugawarepal <- c('#c51b7d','#e6f598',
#                         '#FF0000','#fc8d59',"#FFCC00",'#FFF000',
#                         '#000CCC','#99d594','#3288bd',"#654CFF")
familydrugawarepal <- c('#4575b4','#ffffbf','#d73027','#f46d43','#fdae61','#fee090','#74add1','#000CCC',"#654CFF")

### colour by aware within beta-lactam
data_bybugdrug <- data_long %>% 
  group_by(pathogen, drug, family, family_aware, age, gender) %>% 
  summarise(sus_t = sum(sus), res_t = sum(res)) %>% 
  group_by(pathogen, drug, family, family_aware, age, gender) %>% mutate(proportion = res_t / (sus_t + res_t), total = sum(sus_t+res_t)) %>%
  filter(total > 10)
data_bybugdrug$family <- factor(data_bybugdrug$family)
data_bybugdrug$family_aware <- factor(data_bybugdrug$family_aware)

g1a <- 
  ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), 
         aes(x=age, y = proportion, group = gender)) + 
  geom_point(aes(size = total, shape = gender),alpha = 0.2, colour = "black") + #geom_point(aes(size = total, alpha = 0.5, colour = family_aware)) + 
  geom_smooth(aes(fill = family_aware, linetype = gender), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_linetype_manual(name="Sex", values = c(1, 2), labels = c("Female","Male"),
                        guide = guide_legend(override.aes = list(size = 10))) +
  scale_fill_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_size("# samples", breaks = seq(2500,20000,2500), limits = c(0,12000), guide = guide_legend(override.aes = list(alpha = 0.3) )) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none", fill = "none") + 
  facet_grid(pathogen ~ drug)+ 
  theme(strip.background = element_rect(fill= "white", size=1.5))

g2a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Klebsiella pneumoniae","Escherichia coli")), aes(x=age, y = proportion, group = gender)) + 
  geom_point(aes(size = total, shape = gender),alpha = 0.05, colour = "black") + #geom_point(aes(size = total, alpha = 0.5, colour = family_aware)) + 
  geom_smooth(aes(fill = family_aware, linetype = gender), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_linetype_manual(name="Sex", values = c(1, 2), labels = c("Female","Male"),
                        guide = guide_legend(override.aes = list(size = 10))) +
  scale_fill_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_size("# samples", breaks = seq(2500,20000,2500), limits = c(0,12000), guide = guide_legend(override.aes = list(alpha = 0.3) )) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none", fill = "none") + 
  facet_grid(pathogen ~ drug)+ 
  theme(strip.background = element_rect(fill= "white", size=1.5))

g3a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Enterococcus faecalis","Enterococcus faecium")), aes(x=age, y = proportion, group = gender)) + 
  geom_point(aes(size = total, shape = gender),alpha = 0.2, colour = "black") + #geom_point(aes(size = total, alpha = 0.5, colour = family_aware)) + 
  geom_smooth(aes(fill = family_aware, linetype = gender), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_linetype_manual(name="Sex", values = c(1, 2), labels = c("Female","Male"),
                        guide = guide_legend(override.aes = list(size = 10))) +
  scale_fill_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_size("# samples", breaks = seq(2500,20000,2500), limits = c(0,12000), guide = guide_legend(override.aes = list(alpha = 0.3) )) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none", fill = "none") + 
  facet_grid(pathogen ~  drug)+ 
  theme(strip.background = element_rect(fill= "white", size=1.5))

g4a <- ggplot(data_bybugdrug %>% filter(pathogen %in% c("Staphylococcus aureus", "Streptococcus pneumoniae")), aes(x=age, y = proportion, group = gender)) + 
  geom_point(aes(size = total, shape = gender),alpha = 0.2, colour = "black") + #geom_point(aes(size = total, alpha = 0.5, colour = family_aware)) + 
  geom_smooth(aes(fill = family_aware, linetype = gender), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_linetype_manual(name="Sex", values = c(1, 2), labels = c("Female","Male"),
                        guide = guide_legend(override.aes = list(size = 10))) +
  scale_fill_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_size("# samples", breaks = seq(2500,20000,2500), limits = c(0,12000), guide = guide_legend(override.aes = list(alpha = 0.3) )) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none", fill = "none") + 
  facet_grid(pathogen ~  drug) + 
  theme(strip.background = element_rect(fill= "white", size=1.5))

((g1a + g2a + g3a + g4a + plot_layout(guides = "collect", width = c(1,1.2))) +
  plot_layout(guides = "collect", width = c(1,1.2))) + 
  plot_annotation("European level (2015-2019)", tag_levels = "A") 

# p1a <- (g1a + g2a  + plot_layout(guides = "collect", width = c(1,1.2))) + 
#   plot_annotation("European level (2015-2019)", subtitle = "Gram negative")
# 
# p2a <- (g3a + g4a + plot_layout(guides = "collect", width = c(1,1.2))) + 
#   plot_annotation(subtitle = "Gram positive")
# 
# ((p1a) / (p2a)) + plot_layout(guides = "collect")

ggsave("plots/Fig1.tiff",dpi = 300, width = 22, height = 10)

ggplot(data_bybugdrug %>% filter(pathogen %in% c("Pseudomonas aeruginosa","Acinetobacter spp")), 
       aes(x=age, y = proportion, group = gender)) + 
  geom_point(aes(size = total, shape = gender),alpha = 0.2, colour = "black") + #colour = family_aware, 
  geom_smooth(aes(fill = family_aware, linetype = gender), col = "black") + 
  scale_x_continuous("Age") + 
  scale_y_continuous("Proportion of isolates tested\n that are resistant") + 
  scale_color_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_linetype_manual(name="Sex", values = c(1, 2), labels = c("Female","Male"),
                        guide = guide_legend(override.aes = list(size = 10))) +
  scale_fill_manual("Drug family\n(AWaRE\ngrouping)",values = familydrugawarepal, drop = FALSE) + 
  scale_size("# samples", breaks = seq(2500,20000,2500), limits = c(0,12000), guide = guide_legend(override.aes = list(alpha = 0.3) )) + 
  theme(axis.text.x = element_text(angle = 90),strip.text.y = element_text(face = "italic")) + 
  guides(colour = guide_legend(override.aes = list(size=10)), alpha = "none", fill = "none") + 
  facet_grid(pathogen ~ drug)+ 
  theme(strip.background = element_rect(fill= "white", size=1.5))

