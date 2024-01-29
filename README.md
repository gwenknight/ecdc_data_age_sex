# ecdc_data_age_sex
Analysis of ECDC EARS-NET data for age and sex patterns
##### Title: Antimicrobial resistance prevalence in bloodstream infection in 29 European countries by age and sex: an observational study 
##### Authors: Naomi R Waterlow, Ben S Cooper, Julie V Robotham, Gwenan M Knight
##### Code developed by Naomi Waterlow and Gwenan M Knight

##### Data needed: 
This analysis used European Antimicrobial Resistance Surveillance Network (EARS-Net) patient level data from the Surveillance System (TESSy), 
provided by Austria, Belgium, Bulgaria, Cyprus, Czechia, Germany, Denmark, Estonia, Greece, Spain, Finland, France, Croatia, Hungary, 
Ireland, Iceland, Italy, Luxembourg, Latvia, Malta, Netherlands, Norway, Poland, Portugal, Romania, Sweden, Slovenia, 
Slovakia and the United Kingdom and released by ECDC. 

This is available upon application to: https://www.ecdc.europa.eu/en/publications-data/european-surveillance-system-tessy

# File structure needed: 
data/

output/

plots/

## Stages: 
Step 1_
Data was cleaned in line with the steps in the methods and S1 Appendix to determine the final cleaned dataset. 

Step 2_ 
Calculations of variation in incidence of total infection numbers 
First results figure 

Step 3_
Fitting models in R by age and sex. 

Step 4_
Analysis of missing data and sensitivity analysis
