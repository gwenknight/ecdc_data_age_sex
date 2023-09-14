#### Load in libraries and set plotting rules
library(tidyverse)
library(data.table)
library(lme4)
library(boot)
library(lmtest)
theme_set(theme_bw(base_size = 11))

####### Read in data ######
#data_orig <- read_csv("data/ecdc_ears_net_2002.csv")[,-1] # Apply to TESSY for access
dim(data_orig) # Original = 3,549,617 isolates
# Make sure numeric values correct
data <- data_orig
data$age <- as.numeric(data$age)
data <- data %>% filter(!is.na(age)) %>% filter(gender %in% c("f","m"))

str(data_orig)

# convert to data.table as easier for me!
data_dt <- as.data.table(data)


##### Initial MRSA investigation on hospital and lab ids #####

# choose one bug to look at for now. # mrsa_resistance. 
data_dt <- data_dt[,c("DateUsedForStatisticsYear", "reportingcountry", "gender",
                      "hospitalid", "hospitalunittype", "laboratorycode", "pathogen", "mrsa_R", 
                      "age")]
# omit NAS, only mrsa and only years of interest
data_dt_mrsa <- na.omit(data_dt[pathogen == "staaur" & 
                                  DateUsedForStatisticsYear %in% c("2019", "2018", 
                                                                   "2017", "2016", 
                                                                   "2015"),])

# by lab id
data_by_lab <- data_dt_mrsa[, mean(mrsa_R), by = c("laboratorycode", "reportingcountry", "DateUsedForStatisticsYear")]
data_by_lab$N <-  data_dt_mrsa[, .N, by = c("laboratorycode", "reportingcountry", "DateUsedForStatisticsYear")]$N

ggplot(data_by_lab, aes(x = laboratorycode, y = V1, colour = log(N))) + 
  geom_point() + 
  facet_wrap(reportingcountry~., scales = "free_x") + 
  theme(axis.text.x = element_blank()) + 
  labs(x = "laboratory ID", y = "Proporiton resistant to MRSA")

## definetly looks like it varies by laboratory ID!

# by hospital
data_by_hospital <- data_dt_mrsa[, mean(mrsa_R), by = c("hospitalid", "reportingcountry", "DateUsedForStatisticsYear")]
data_by_hospital$N <-  data_dt_mrsa[, .N, by = c("hospitalid", "reportingcountry", "DateUsedForStatisticsYear")]$N

ggplot(data_by_hospital, aes(x = hospitalid, y = V1, colour = log(N))) + 
  geom_point() + 
  facet_wrap(reportingcountry~., scales = "free_x") + 
  theme(axis.text.x = element_blank()) + 
  labs(x = "Hospital ID", y = "Proporiton resistant to MRSA")

# also looks quite varied!

####### Run the initial country variance model #####
# for now just using the same packages as in course

mrsa_basic <- glm(mrsa_R ~ 1  , 
                  data = data_dt_mrsa, family = binomial(link = "logit"))
summary(mrsa_basic)

mrsa_country <- glmer(mrsa_R ~ 1  + (1 | reportingcountry), 
                      data = data_dt_mrsa, family = binomial(link = "logit"))
summary(mrsa_country)

# Perform likelihood ratio test
lrtest(mrsa_basic, mrsa_country)

# Calculate the VPC=ICC
rpm2 <- as.data.frame(VarCorr(mrsa_country))

rpm2$vcov[rpm2$grp == "reportingcountry"] / (rpm2$vcov[rpm2$grp == "reportingcountry"] + pi^2/3)
# so 30% of the variance is explained by country? Lots of clustering!

#Extract the estimated intercept
beta0 <- fixef(mrsa_country)[[1]]

# Extract the state effects using ranef()
u0 <- as.data.frame(ranef(mrsa_country, condVar = TRUE))
# The point estimates are stored in 'condval'. The SEs are stored in 'condsd'.

# Calculate predicted state probabilities of voting republican
u0$pstate <- inv.logit(beta0 + u0$condval)

# Calculate the lower and upper limits of the 95% confidence intervals.
u0$plower <- inv.logit(beta0 + u0$condval - 1.96 * u0$condsd)
u0$pupper <- inv.logit(beta0 + u0$condval + 1.96 * u0$condsd)

# Rank states according to conditional probabilities of voting republican.
u0$rank <- rank(u0$pstate)

# Plot a caterpillar plot of the predicted state effects with their 95%
# confidence intervals
ggplot(u0, aes(x = rank, y = pstate)) +
  geom_hline(yintercept = mean(u0$pstate)) +
  geom_point() +
  geom_errorbar(aes(ymin = plower, ymax = pupper)) +
  geom_text(aes(label = grp),
            position = position_dodge(width = 0.1),
            vjust = 6,
            size = 2) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, by = 0.2)) +
  labs(x = "Rank", y = expression(inv.logit(intercept + u[j]))) + 
  labs( y="predicted probability of testing resistant")+
  theme(axis.title = element_text(size = 17), 
        title = element_text(size = 17))


###### run variance by lab model ######

mrsa_lab_country <- glmer(mrsa_R ~ 1  + (1 | laboratorycode) + (1 | reportingcountry), 
                          data = data_dt_mrsa, family = binomial(link = "logit"))
summary(mrsa_lab_country)

# Perform likelihood ratio test
lrtest(mrsa_country, mrsa_lab_country)

# Calculate the VPC=ICC
rpm_lc <- as.data.frame(VarCorr(mrsa_lab_country))

# data frame using the VarCorr() function of the lme4 package.
rpm_lc$vcov[rpm_lc$grp == "reportingcountry"] / (sum(rpm_lc$vcov) + pi^2/3)
# this is the amount explained by country?

# lab level VPC (but not ICC)
rpm_lc$vcov[rpm_lc$grp == "laboratorycode"] / (sum(rpm_lc$vcov) + pi^2/3)

# country and lab combined VPC (or ICC)
(rpm_lc$vcov[rpm_lc$grp == "reportingcountry"] + rpm_lc$vcov[rpm_lc$grp == "laboratorycode"]) /
  (sum(rpm_lc$vcov)+ pi^2/3)

res <- as.data.frame(ranef(mrsa_lab_country, condVar = TRUE))
beta0_lc <- fixef(mrsa_lab_country)[[1]]
# Next we restrict our focus to the lab random effects from res
u0_lab <- res[res$grpvar == "laboratorycode", ]

# Calculate predicted state probabilities of voting republican
u0_lab$plab <- inv.logit(beta0_lc + u0_lab$condval)

# Calculate the lower and upper limits of the 95% confidence intervals.
u0_lab$plower <- inv.logit(beta0_lc + u0_lab$condval - 1.96 * u0_lab$condsd)
u0_lab$pupper <- inv.logit(beta0_lc + u0_lab$condval + 1.96 * u0_lab$condsd)
# Next rank the lab random effects
u0_lab$rank <- rank(u0_lab$condval)

# Plot a caterpillar plot of the lab effects using confidence interval
# values.
LAB_CATERPILLAR <- ggplot(u0_lab, aes(x = rank, y = plab, ymin = plower, ymax = pupper)) +
  geom_errorbar() +
  geom_point(colour = "darkgreen") +
  geom_hline(yintercept = mean(u0_lab$plab), colour = "darkblue") +
  labs(x = "Rank", y = "Predicted probability of resistance", title = "B: Laboratory variation") + 
  theme(axis.title = element_text(size = 17))

LAB_CATERPILLAR

# Next we restrict our focus to the lab random effects from res
u0_country <- res[res$grpvar == "reportingcountry", ]

# Calculate predicted state probabilities of voting republican
u0_country$pcountry <- inv.logit(beta0_lc + u0_country$condval)

# Calculate the lower and upper limits of the 95% confidence intervals.
u0_country$plower <- inv.logit(beta0_lc + u0_country$condval - 1.96 * u0_country$condsd)
u0_country$pupper <- inv.logit(beta0_lc + u0_country$condval + 1.96 * u0_country$condsd)
# Next rank the lab random effects
u0_country$rank <- rank(u0_country$condval)

# Plot a caterpillar plot of the lab effects using confidence interval
# values.
COUNTRY_CATERPILLAR <- ggplot(u0_country, aes(x = rank, y = pcountry, ymin = plower, ymax = pupper)) +
  geom_errorbar() +
  geom_point(colour = "darkgreen") +
  geom_hline(yintercept = mean(u0_lab$plab), colour = "darkblue") +
  labs(x = "Rank", y = "Predicted probability of resistance", title = "A: Country variation")+ 
  theme(axis.title = element_text(size = 17)) + 
  geom_text(aes(label = grp, y = plower, vjust = 1.5))

COUNTRY_CATERPILLAR

tiff(paste0("plots/variance_components_mrsa.tiff"), height = 3000, width = 3000, res = 300)
grid.arrange(COUNTRY_CATERPILLAR, LAB_CATERPILLAR, ncol = 1)
dev.off()

