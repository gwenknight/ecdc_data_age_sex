# comparing across models
library(data.table)
library(ggplot2)
library(readr)
colours_to_use <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
#### LOAD IN THE DATA #####
# combinations
path_drug_combos <- data.frame(read.csv2("path_drug_combos.csv", sep = ",", header = F))
country_reference <- as.data.table(read.csv("countries_reference.csv"))
colnames(country_reference) <- c("country", "country_long")

##### the overall ones ####

# search for files
files_to_search <- list.files("output")
# select those named right
relevant_files <- files_to_search[grep(pattern = "raw_comp", 
                                       x = files_to_search)]

comparisons_all <- data.table()
# read each one and rbind
for(i in 1:length(relevant_files)){
  temp <- read_rds(file = paste0("output/",relevant_files[i]))
  comparisons_all <- rbind(comparisons_all, temp)
}

comparisons_all[, bug_drug := paste0(bug, "_", resistance)]


#### country specific #### 

# search for files
files_to_search <- list.files("output")
# select those named right
relevant_files <- files_to_search[grep(pattern = "comparison_", 
                                       x = files_to_search)]

comparisons_1_100 <- data.table()
# read each one and rbind
for(i in 1:length(relevant_files)){
  temp <- read_rds(file = paste0("output/",relevant_files[i]))
  comparisons_1_100 <- rbind(comparisons_1_100, temp)
  print(i)
}

comparisons_1_100[, bug_drug := paste0(bug, "_", resistance)]

#### needed function ###

# function for use later (to calcualte comparisons by charactersic)
calculate_comparisons <- function(data_in, summary_level, subset_required){
  
  f <- paste0(summary_level, " ~ type")
  
  if(summary_level == "bug"){
    data_in <- data_in[resistance == subset_required,]} else if(
      summary_level == "resistance"
    ) {
      data_in <- data_in[bug == subset_required,]
    }
  
  summary_calculations <- data_in[, quantile(female_difference, probs = c(0.025, 0.5, 0.975)), by = c(summary_level)]
  colnames(summary_calculations)[2] <- "female_difference"
  summary_calculations$male_difference <- data_in[, quantile(male_difference, probs = c(0.025, 0.5, 0.975)), by = c(summary_level)]$V1
  summary_calculations$type <- rep(c("lower", "median", "upper"), nrow(unique(summary_calculations[,..summary_level])))
  female <- dcast.data.table(summary_calculations, f, value.var = "female_difference")
  male <- dcast.data.table(summary_calculations, f, value.var = "male_difference")
  female$gender <- "female"
  male$gender <- "male"
  
  combo <- rbind(female, male)
  combo <- combo[order(gender, median)]
  
  return(combo)
}


comparisons_1_100[country_reference, on = c("country"), country_long := i.country_long]


###### LOOK AT COUNTRY SPECIFICS #####

# all_countries <- c("BG", "ES", "FR", "HU", "IE", "IT", "PT", "SK", "UK", "AT",
#                    "EL", "NL", "HR", "LV", "PL", "SE", "CY", "NO", "BE", "CZ",
#                    "LU", "SI", "DK", "FI", "DE", "MT", "RO", "EE", "IS")

for(i in unlist(country_reference[,"country_long"])){
  
  country_of_interest <- i
  
  COUNTRY_PLOT <- ggplot(comparisons_1_100[country_long == country_of_interest],
                         aes(x = bug_drug, y = median, colour = gender)) + 
    geom_pointrange(aes(ymin= lower, ymax = upper),
                    position = position_dodge(width = 0.3)) + 
    labs(y = "Change in proportion between age 1 and 100",
         x = "bacteria - antibiotic",
         title = country_of_interest) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    geom_hline(yintercept = 0, linetype ="dashed") + 
    theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    lims(y=c(-0.5, 0.5))
  
  tiff(paste0("output_figures/summaries_country/country_plots_", country_of_interest, ".tiff"), 
       width = 3250, height = 2000, res = 300)
  print(COUNTRY_PLOT)
  dev.off()
  
  print(paste0("Country ", country_of_interest, " completed"))
  
}

##### LOOK AT BUG SPECIFICS #####

for(i in unique(path_drug_combos[,1])){
  
  bug_of_interest <- i
  
  BUG_PLOT <- ggplot(comparisons_1_100[bug == bug_of_interest],
                     aes(x = resistance, y = median, colour = gender, group = country)) + 
    geom_pointrange(aes(ymin= lower, ymax = upper),
                    position = position_dodge(width = 0.8)) + 
    labs(y = "Change in proportion between age 1 and 100",
         x = "antibiotic",
         title = bug_of_interest) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(gender ~.) + 
    geom_hline(yintercept = 0, linetype ="dashed") + 
    theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    lims(y=c(-0.5, 0.5))
  
  tiff(paste0("output_figures/summaries_country/bug_plots_", bug_of_interest, ".tiff"), 
       width = 3250, height = 2000, res = 300)
  print(BUG_PLOT)
  dev.off()
  
  combo_bug <- calculate_comparisons(comparisons_all, "resistance", bug_of_interest)
  
  
  BUG_COMPARE <- ggplot(combo_bug,
                        aes(x = resistance, y = median, colour = gender)) + 
    geom_pointrange(aes(ymin= lower, ymax = upper),
                    position = position_dodge(width = 0.3)) + 
    labs(y = "Change in proportion between age 1 and 100",
         x = "antibiotic",
         title = bug_of_interest) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    geom_hline(yintercept = 0, linetype ="dashed") + 
    theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    lims(y=c(-0.5, 0.5))
  
  tiff(paste0("output_figures/summaries_combined/bug_compare_", bug_of_interest ,".tiff"), 
       width = 3250, height = 2000, res = 300)
  print(BUG_COMPARE)
  dev.off()
  
  
  
  print(paste0("Bug ", bug_of_interest, " completed"))
  
}
##### LOOK AT DRUG SPECIFICS #####


for(i in unique(path_drug_combos[,2])){
  
  drug_of_interest <- i
  
  DRUG_PLOT <- ggplot(comparisons_1_100[resistance == drug_of_interest],
                      aes(x = bug, y = median, colour = gender, group = country_long)) + 
    geom_pointrange(aes(ymin= lower, ymax = upper),
                    position = position_dodge(width = 0.8)) + 
    labs(y = "Change in proportion between age 1 and 100",
         x = "bacteria",
         title = drug_of_interest) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    facet_grid(gender ~.) + 
    geom_hline(yintercept = 0, linetype ="dashed") + 
    theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    lims(y=c(-0.5, 0.5))
  
  
  tiff(paste0("output_figures/summaries_country/drug_plots", drug_of_interest, ".tiff"), 
       width = 3250, height = 2000, res = 300)
  print(DRUG_PLOT)
  dev.off()
  
  combo_drug <- calculate_comparisons(comparisons_all, "bug", drug_of_interest)
  
  
  DRUG_COMPARE <- ggplot(combo_drug,
                         aes(x = bug, y = median, colour = gender)) + 
    geom_pointrange(aes(ymin= lower, ymax = upper),
                    position = position_dodge(width = 0.3)) + 
    labs(y = "Change in proportion between age 1 and 100",
         x = "antibiotic",
         title = drug_of_interest) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    geom_hline(yintercept = 0, linetype ="dashed") + 
    theme_linedraw() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
    lims(y=c(-0.5, 0.5))
  
  tiff(paste0("output_figures/summaries_combined/drug_compare_", drug_of_interest ,".tiff"), 
       width = 3250, height = 2000, res = 300)
  print(DRUG_COMPARE)
  dev.off()
  
  
  
  
  print(paste0("Drug ", drug_of_interest, " completed"))
  
}

##### look at the correlations between fixed effects #####
fixed_effs_storage <- readRDS("output/fixed_effs_storage.RDS")

AGES_POSTERIORS <- ggplot(fixed_effs_storage, aes(x = age_s, y = age_squared_s, colour = bug_drug)) + 
  geom_point(alpha = 0.1, shape =16) + 
  theme_linedraw() + 
  geom_hline(yintercept = 0, linetype="dashed") + 
  geom_vline(xintercept = 0, linetype="dashed") 

tiff(paste0("output_figures/summaries_combined/age_posteriors.tiff"), 
     width = 3250, height = 3000, res = 300)
print(AGES_POSTERIORS)
dev.off()


fixed_effs_summmary <- fixed_effs_storage[, quantile(Intercept, probs = c(0.025, 0.5, 0.975)), by = bug_drug]
colnames(fixed_effs_summmary)[2] <- "Intercept"
fixed_effs_summmary$age <- fixed_effs_storage[, quantile(age_s, probs = c(0.025, 0.5, 0.975)), by = bug_drug]$V1
fixed_effs_summmary$age_squared <- fixed_effs_storage[, quantile(age_squared_s, probs = c(0.025, 0.5, 0.975)), by = bug_drug]$V1
fixed_effs_summmary$genderm <- fixed_effs_storage[, quantile(genderm, probs = c(0.025, 0.5, 0.975)), by = bug_drug]$V1
fixed_effs_summmary$year_s <- fixed_effs_storage[, quantile(year_s, probs = c(0.025, 0.5, 0.975)), by = bug_drug]$V1
fixed_effs_summmary$age_genderm <- fixed_effs_storage[, quantile(`age_s:genderm`, probs = c(0.025, 0.5, 0.975)), by = bug_drug]$V1
fixed_effs_summmary$type <- rep(c("lower", "median", "upper"), 33)
fixed_effs_summmary_m <- melt.data.table(fixed_effs_summmary, id.vars = c("bug_drug", "type"), measure.vars = c("Intercept", 
                                                                                                                "age", 
                                                                                                                "age_squared", 
                                                                                                                "genderm", 
                                                                                                                "year_s", 
                                                                                                                "age_genderm") )
fixed_effs_summmary_c <- dcast.data.table(fixed_effs_summmary_m, bug_drug + variable ~ type, 
                                          value.var = "value")

fixed_effs_summmary_c[, bug := strsplit(bug_drug, "[_]")[[1]][1]]
fixed_effs_summmary_c[, c("bug", "drug", "leftover1", "leftover2") := tstrsplit(bug_drug, "_", fixed=TRUE)]
fixed_effs_summmary_c[,c("leftover1", "leftover2") := NULL]

FIXED_EFFS_VALUE <- ggplot(fixed_effs_summmary_c, aes( y =  bug, x = median, colour = drug, group = bug_drug))+ 
  geom_vline(xintercept = 0, linetype = "dashed") + 
  # geom_point(position = position_dodge(width=0.75)) + 
  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width=0.5), size = 0.2) + 
  facet_wrap(variable~., scales= "free_x") +
  labs(x = "Bacteria", y = "Posterior value", colour = "Drug") 


tiff(paste0("output_figures/summaries_combined/fixed_effects.tiff"), 
     width = 2000, height = 2000, res = 300)
FIXED_EFFS_VALUE
dev.off()

##### plotting tests #####

temp <- comparisons_1_100[bug_drug == "staaur_mrsa_R",]
new_order <- unique(temp[order(gender, median)]$country_long)

comparisons_1_100$country_long <- factor(comparisons_1_100$country_long , 
                                         levels = new_order)

translate_bug_drugs <-as.data.table(read.csv("data/translate_drug_bugs.csv", header = F))
colnames(translate_bug_drugs) <- c("bug_long", "bug", "resistance", "V4", "V5", "drug_long")
comparisons_1_100[translate_bug_drugs, on ="bug", bug_long := i.bug_long]
comparisons_1_100[translate_bug_drugs, on ="resistance", resistance_long := i.drug_long]
comparisons_1_100[,bug_drug_long := paste0(bug_long, "\n", resistance_long)]


SPECIFIC <- ggplot(comparisons_1_100[gender == "female"], aes(x = country_long, y = median, colour = gender)) + 
  facet_wrap(bug_drug_long~., ncol = 12) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.3), 
                  size = 0.2) +
  geom_point() + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4), 
        legend.position = "none") + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  scale_color_manual(values = colours_to_use) +
  lims(y=c(-0.5, 0.5)) + 
  labs(x = "Country", y = "Change in proportion for ages 1-100", 
       colour = "gender")


tiff(paste0("output_figures/summaries_combined/bug_drug_comps.tiff"), 
     width = 6000, height = 3000, res = 300)
SPECIFIC
dev.off()



COUNTRY_SPECIFIC <- ggplot(comparisons_1_100[gender =="female"], aes(y = bug, x = median, colour = resistance)) + 
  facet_wrap(.~country_long, ncol = 10) + 
  geom_pointrange(aes(xmin= lower, xmax = upper),
                  position = position_dodge(width = 0.7), 
                  size = 0.2) +
  theme_linedraw() + 
  geom_point() + 
  # theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_vline(xintercept = 0, linetype ="dashed") + 
  lims(x=c(-0.5, 0.5)) + 
  labs(x = "bacteria - antibiotic", y = "Change in proportion ages 1-100 (female only)", 
       colour = "resitance")

tiff(paste0("output_figures/summaries_combined/by_country.tiff"), 
     width = 5000, height = 3000, res = 300)
COUNTRY_SPECIFIC
dev.off()

###### SUBREGION RANKING ####


sub_regions <- data.table(read.csv("subregions.csv"))

for_ranking <- comparisons_1_100[, quantile(male_difference, probs = 0.5), by = c("bug_drug", "country")]
for_ranking <- comparisons_1_100[, abs(quantile(male_difference, probs = 0.5)), by = c("bug_drug", "country")]
for_ranking[, ranking := rank(V1), by = bug_drug]

ranking_out <- for_ranking[, sum(ranking), by = "country"]


ranking_out[sub_regions, on="country", subregion := i.subregion]
ranking_out[country_reference, on = "country", country_long := i.country_long]
new_order <- ranking_out[order(V1, country_long)]$country_long
ranking_out$country_long <- factor(ranking_out$country_long, levels = new_order)

ggplot(ranking_out, aes(x = country_long, y = V1, group = country, colour = subregion)) + 
  geom_point() + theme_linedraw() + 
  labs(y = "Combined ranking across bug_drugs (absolute)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


combine_country_only <- data.table(read.csv( "output/incidence_country_level.csv"))
cc_2019 <- combine_country_only[year == 2019]

ranking_out[cc_2019, on = "country", incidence := i.incidence]


ggplot(ranking_out, aes(x = incidence, y = V1, group = country, colour = subregion)) + 
  geom_point() + theme_linedraw() + 
  labs(y = "Combined ranking across bug_drugs (absolute)", x = "incidence") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



###### correlation between age trends ######


testing <- dcast.data.table(comparisons_all, Sample + country  ~ bug_drug , value.var = "female_100")
testing <- testing[Sample <= 6000,3:33] # get the same number of samples for all.
cor_mat <- data.table(cor(testing))
cor_mat[upper.tri(cor_mat)] <- NA
cor_mat$bug_drug <- colnames(cor_mat)
cor_mat_m <- melt.data.table(cor_mat, id.vars = "bug_drug")
cor_mat_m <- na.omit(cor_mat_m)


FEMALE_COR <- ggplot(cor_mat_m, aes(x = bug_drug, y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#f1a340", high = "#998ec3", mid = "white", 
                       midpoint = 0, ) + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label = round(value, 2)), size = 1.2) + 
  labs(x= "", y ="", fill ="Correlation", title = "female") 


testing <- dcast.data.table(comparisons_all, Sample + country  ~ bug_drug , value.var = "male_100")
testing <- na.omit(testing)
testing <- testing[,3:35]
cor_mat <- data.table(cor(testing))
cor_mat[upper.tri(cor_mat)] <- NA
cor_mat$bug_drug <- colnames(cor_mat)
cor_mat_m <- melt.data.table(cor_mat, id.vars = "bug_drug")
cor_mat_m <- na.omit(cor_mat_m)


MALE_COR <- ggplot(cor_mat_m, aes(x = bug_drug, y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#f1a340", high = "#998ec3", mid = "white", 
                       midpoint = 0, ) + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label = round(value, 2)), size = 1.2) + 
  labs(x= "", y ="", fill ="Correlation", title = "male") 

library(gtable)
library(gridExtra)
legend = gtable_filter(ggplot_gtable(ggplot_build(MALE_COR)), "guide-box")

tiff(paste0("output_figures/summaries_combined/correlations.tiff"), 
     width = 3250, height = 3000, res = 300)

print( grid.arrange(
  FEMALE_COR + theme(legend.position="none"),
  MALE_COR + theme(legend.position="none"),
  legend, 
  layout_matrix = rbind(c(1,1,1,1,1,3),
                        c(2,2,2,2,2,3))
  
))
dev.off()


# but what I really want to know are which bugs/drugs are correlation across countries!
# i.e. it's the by country bit i need to get at somehow
# d0nt care about sample - "FEMALE"
testing <- dcast.data.table(comparisons_all, Sample + country  ~ bug_drug , value.var = "female_100")
testing <- testing[Sample <= 4,]
testing <- melt.data.table(testing, id.vars = c("country", "Sample"))
temp <- nrow(testing)

expanded_version <- testing[rep(seq_len(nrow(testing)), each = length(unique(testing$variable))), ]
colnames(expanded_version)[3] <- "variable1"
expanded_version$variable <- rep(unique(testing$variable), temp)
# now link the values back in
expanded_version[testing, on = c("country", "Sample", "variable"), value2 := i.value]
colnames(expanded_version) <- c("country", "sample", "variable1", "value1", 
                                "variable2", "value2")

#remove duplicated
de_duped <- expanded_version[ !duplicated(apply(expanded_version[,c("country","sample","variable1", "variable2")], 1, sort), MARGIN = 2), ]

# calculate correlation
correlation_output <- de_duped[,  cor(value1, value2), by = c("country", "variable1", "variable2")]

correlation_output[, c("bug1", "drug1", "leftover1", "leftover2") := tstrsplit(variable1, "_", fixed=TRUE)]
correlation_output[,c("leftover1", "leftover2") := NULL]
correlation_output[, c("bug2", "drug2", "leftover1", "leftover2") := tstrsplit(variable2, "_", fixed=TRUE)]
correlation_output[,c("leftover1", "leftover2") := NULL]

ggplot(correlation_output, aes(x = variable1, y = country, fill= V1)) + 
  facet_wrap(variable2~.) + 
  geom_tile() +
  #  geom_jitter(alpha = 0.5, shape =16) +
  scale_fill_gradient2(low = "#f1a340", high = "#998ec3", mid = "white", 
                       midpoint = 0, ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  #geom_text(aes(label = round(V1, 2)), size = 1.2) + 
  labs(x= "", y ="", fill ="Correlation", title = "FEMALE") 


