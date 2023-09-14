# File for loading things

# libraries
library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(posterior)
library(patchwork)
library(rstanarm)
library(rstan)
library(loo)
library(boot) # for inv.logit
library(Hmisc) # for binomial confidence intervals
library(rcartocolor) # for colours 
library(data.table)
library(grid)
library(gridExtra)
library(gtable)
theme_set(theme_bw(base_size = 11))
colours_to_use <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

how_many_cores_free_default <- 3

if(!exists("how_many_cores_free")){
  how_many_cores_free <- how_many_cores_free_default
} 

rstan_options(auto_write = TRUE)
options(mc.cores = (parallel::detectCores()-how_many_cores_free))