## light penetration in ambient vs tie-back plots
## Author: Emily H
## Created: June 6, 2024
## Last edited: August 13, 2024

#install.packages("tidyverse")
#install.packages("lme4)
#install.packages("lmerTest")
#install.packages("car")
#install.packages("emmeans")


library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(emmeans)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

#### import data ####
#read in PAR data
PAR21.new <- read_rds("output/PAR21 new.rds")
str(PAR21.new)

# read in trt data
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X)
str(trts)

# merge dfs
dat <- merge(trts, PAR21.new, by = "Plot") %>% 
  filter(Thin == "Not thinned")
str(dat)

# calculate means and SEs
sum.dat <- dat %>%
  dplyr::group_by(Light) %>%
  dplyr::summarise(mean.light.penetration = mean(True.light.penetration),
                   se.light.penetration = sd(True.light.penetration)/sqrt(length(True.light.penetration))) %>%
  mutate_if(is.numeric, round, digits = 3)
View(sum.dat)

# test for significant differences
light.model <- lmer(True.light.penetration ~ Light + (1|Block), data = dat)
Anova(light.model)
emmeans(light.model, pairwise ~ Light, adjust = "Tukey") #differences between tie-backs and ambient/shade plots

