# model comparisons from raw predictions
library(data.table)
library(ggplot2)
library(readr)
library(gridExtra)

colours_to_use <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
familydrugpal <- c('#c51b7d','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd')
# familydrugawarepal <- c('#c51b7d','#e6f598',
#                         '#FF0000','#fc8d59',"#FFCC00",'#FFF000',
#                         '#000CCC','#99d594','#3288bd',"#654CFF")
familydrugawarepal <- c('#4575b4','#ffffbf','#d73027','#f46d43','#fdae61','#fee090','#74add1','#000CCC',"#654CFF")
familydrugawarepal <- c('#c51b7d',
                        '#FF0000','#fc8d59',"#FFCC00",'#FFF000',
                        '#000CCC','#99d594','#3288bd','#e6f598')

translate_bug_drugs <-as.data.table(read.csv("data/translate_drug_bugs.csv", header = F))
translate_bugs <- unique(translate_bug_drugs[,c("V2", "V1")])
colnames(translate_bugs) <- c("bug", "bug_long")
translate_drugs <- unique(translate_bug_drugs[,c("V3", "V6")])
colnames(translate_drugs) <- c("resistance", "drug")

translate_bugs[bug_long == "Pseudomonas aeruginosa", Gram := "Negative"]
translate_bugs[bug_long == "Klebsiella pneumoniae", Gram := "Negative"]
translate_bugs[bug_long == "Streptococcus pneumoniae", Gram := "Positive"]
translate_bugs[bug_long == "Escherichia coli", Gram := "Negative"]
translate_bugs[bug_long == "Enterococcus faecium", Gram := "Positive"]
translate_bugs[bug_long == "Staphylococcus aureus", Gram := "Positive"]
translate_bugs[bug_long == "Enterococcus faecalis", Gram := "Positive"]
translate_bugs[bug_long == "Acinetobacter species", Gram := "Negative"]




translate_drugs[resistance %in% c("amika_R", "aminogl_R","genta_high"), drug_family := "aminoglycoside"]
translate_drugs[resistance %in% c("carbapen_R", "cefIII_entero_R", "ert_R", "ureidopen_R", "ceftaz_R", 
                                  "mrsa_R", "penic_RI", "aminopen_R"), drug_family := "beta-lactam"]
translate_drugs[resistance %in% c("fq_pseudo_R", "fq_ent_R", "fq_staaur_R"), drug_family :="fluoroquinolone"]
translate_drugs[resistance %in% c("vanco_R"), drug_family :="glycopeptide"]
translate_drugs[resistance %in% c("macrol_R"), drug_family :="macrolides"]
translate_drugs[resistance %in% c("rifamp_R"), drug_family :="rifamycin"]

translate_drugs[drug  %in% c("Penicillins", "Aminoglycosides"), aware := "Access + Watch"]
translate_drugs[drug  %in% c("Third-generation cephalosporins"), aware := "Watch + Reserve"]
translate_drugs[drug  %in% c("Amikacin","Aminopenicillins","High-level aminoglycoside", "Oxacillin or cefoxitin"), aware := "Access"]
translate_drugs[is.na(aware), aware := "Watch"]

translate_drugs[drug_family == "beta-lactam", final_grouping := paste0(paste(drug_family, aware,sep="\n("),")")]
translate_drugs[drug_family != "beta-lactam", final_grouping := drug_family]


###### LOAD DATA #####

# search for files
files_to_search <- list.files("output")
# select those named right
relevant_files <- files_to_search[grep(pattern = "raw_comp", 
                                       x = files_to_search)]

comparisons_1_100 <- data.table()
# read each one and rbind
for(i in 1:length(relevant_files)){
  temp <- read_rds(file = paste0("output/",relevant_files[i]))
  comparisons_1_100 <- rbind(comparisons_1_100, temp)
}

comparisons_1_100[, bug_drug := paste0(bug, "_", resistance)]


## load function. 

calculate_comparisons <- function(data_in, summary_level){
  
  f <- paste0(summary_level, " ~ type")
  
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


###### BY BUG #####

combo_bug <- calculate_comparisons(comparisons_1_100, "bug")

BUG_COMPARE <- ggplot(combo_bug,
                      aes(x = bug, y = median, colour = gender)) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.3)) + 
  labs(y = "Change in proportion between age 1 and 100",
       x = "antibiotic",
       title = "comparison across bugs") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  lims(y=c(-0.5, 0.5)) + 
  scale_color_manual(values=colours_to_use[1:2]) + 
  theme(axis.text.y=element_text(face="italic"))

tiff(paste0("output_figures/summaries_combined/bug_compare.tiff"), 
     width = 3250, height = 2000, res = 300)
print(BUG_COMPARE)
dev.off()

summary_calculations <- comparisons_1_100[, quantile(female_difference, probs = c(0.025, 0.5, 0.975)), by = c("bug", "resistance")]
colnames(summary_calculations)[3] <- "female_difference"
summary_calculations$male_difference <- comparisons_1_100[, quantile(male_difference, probs = c(0.025, 0.5, 0.975)), by = c("bug", "resistance")]$V1
summary_calculations$type <- rep(c("lower", "median", "upper"), nrow(unique(summary_calculations[,c("bug", "resistance")])))
female <- dcast.data.table(summary_calculations, bug + resistance ~ type, value.var = "female_difference")
male <- dcast.data.table(summary_calculations, bug + resistance ~ type, value.var = "male_difference")
female$gender <- "female"
male$gender <- "male"

combo2 <- rbind(female, male)
combo2 <- combo2[order(gender, median)]

combo2[translate_bugs, on =c("bug"), bug_long := i.bug_long]
combo2[translate_drugs, on =c("resistance"), drug_long := i.drug]
combo2[translate_drugs, on =c("resistance"), final_grouping := i.final_grouping]
combo2[translate_bugs, on =c("bug"), Gram_stain := i.Gram]
# 
# combo2$final_grouping <- factor(combo2$final_grouping, levels = c("aminoglycoside", 
#                                                             "beta-lactam",
#                                                             "fluoroquinolone", 
#                                                             "glycopeptide", 
#                                                             "ansamycin", 
#                                                             "macrolides"))

BUG_COMPARE_DRUG <- ggplot(combo2,
                           aes(x = bug_long, y = median, colour = Gram_stain, group = drug_long)) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.6)) + 
  labs(y = "Change in proportion between age 1 and 100",
       x = "Bacteria",
       title = "A", colour = "Gram stain") + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  theme_linedraw() + 
  facet_grid(.~gender)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  theme(axis.text.y=element_text(face="italic"))+
  # scale_colour_manual(values = familydrugpal) + 
  coord_flip()

DRUG_COMPARE_BUG <- ggplot(combo2,
                           aes(x = drug_long, y = median, colour = final_grouping, group = bug_long)) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.6)) + 
  labs(y = "Change in proportion between age 1 and 100",
       x = "Antibiotic",
       title = "B", 
       colour = "Drug Family") + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  theme_linedraw() + 
  facet_grid(.~gender)+
  scale_colour_manual(values = familydrugawarepal) + 
  coord_flip() + 
  theme(legend.spacing.y = unit(0.2, 'cm')) +
  guides(colour = guide_legend(byrow = TRUE))

#legend = gtable_filter(ggplot_gtable(ggplot_build(DRUG_COMPARE_BUG)), "guide-box")


tiff(paste0("output_figures/summaries_combined/Fig5_original.tiff"), 
     width = 3250, height = 3000, res = 300)
print( grid.arrange(
  BUG_COMPARE_DRUG, # + theme(legend.position = "none"),
  DRUG_COMPARE_BUG, #+ theme(legend.position = "none") ,
  layout_matrix = rbind(c(1,1,1,1,1), 
                        c(2,2,2,2,2))
  
) )

dev.off()

#### Split by gram stain status
combo2 <- as.data.frame(combo2)
combo2$bug_long <- factor(combo2$bug_long, levels = 
                            c("Pseudomonas aeruginosa","Klebsiella pneumoniae","Escherichia coli","Acinetobacter species",
                              "Streptococcus pneumoniae", "Staphylococcus aureus","Enterococcus faecium", "Enterococcus faecalis"))

combo2$final_grouping <- factor(combo2$final_grouping)

gp <- ggplot(combo2 %>% filter(Gram_stain == "Positive"),aes(x = bug_long, y = median, colour = final_grouping, group = drug_long)) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.6)) + 
  facet_grid(.~gender)+
  labs(y = "Change in proportion between age 1 and 100",
       x = "Bacteria",
       title = "Gram positive", 
       colour = "Drug Family") + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  theme_linedraw() +
  coord_flip() + 
  theme(legend.spacing.y = unit(0.2, 'cm')) +
  theme(axis.text.y=element_text(face="italic"))+
  guides(colour = guide_legend(byrow = TRUE)) + 
  scale_colour_manual(values = familydrugawarepal, breaks = c("aminoglycoside", "beta-lactam\n(Access + Watch)", "beta-lactam\n(Access)", 
                                                              "beta-lactam\n(Watch + Reserve)", "beta-lactam\n(Watch)", 
                                                              "fluoroquinolone", "glycopeptide", "macrolides","rifamycin"), drop = FALSE)

gn <- ggplot(combo2 %>% filter(Gram_stain == "Negative"),aes(x = bug_long, y = median, colour = final_grouping, group = drug_long)) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.6)) + 
  facet_grid(.~gender)+
  labs(y = "Change in proportion between age 1 and 100",
       x = "Bacteria",
       title = "Gram negative", 
       colour = "Drug Family") + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  theme_linedraw() + 
  coord_flip() + 
  theme(legend.spacing.y = unit(0.2, 'cm')) +
  theme(axis.text.y=element_text(face="italic"))+
  guides(colour = guide_legend(byrow = TRUE)) + 
  scale_colour_manual(values = familydrugawarepal, breaks = c("aminoglycoside", "beta-lactam\n(Access + Watch)", "beta-lactam\n(Access)", 
                                                              "beta-lactam\n(Watch + Reserve)", "beta-lactam\n(Watch)", 
                                                              "fluoroquinolone", "glycopeptide", "macrolides","rifamycin"), drop = FALSE) 

gn / gp + plot_layout(guides = "collect") 
ggsave("output_figures/Fig5.tiff",dpi = 300, width = 10, height = 10)

###### BY DRUG #####

combo_drug <- calculate_comparisons(comparisons_1_100, "resistance")

DRUG_COMPARE <- ggplot(combo_drug,
                       aes(x = resistance, y = median, colour = gender)) + 
  geom_pointrange(aes(ymin= lower, ymax = upper),
                  position = position_dodge(width = 0.3)) + 
  labs(y = "Change in proportion between age 1 and 100",
       x = "antibiotic",
       title = "comparison across drugs") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept = 0, linetype ="dashed") + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  lims(y=c(-0.5, 0.5))

tiff(paste0("output_figures/summaries_combined/drug_compare.tiff"), 
     width = 3250, height = 2000, res = 300)
print(DRUG_COMPARE)
dev.off()

