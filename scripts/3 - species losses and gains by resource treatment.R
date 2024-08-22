## Spp loss/gain by resource treatment models
## Author: Emily H
## Created: November 28, 2023
## Last edited: August 20, 2024

#install.packages("tidyverse")
#install.packages("lme4)
#install.packages("lmerTest")
#install.packages("car")
#install.packages("nnet")
#install.packages("ggsignif")

library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(nnet)
library(ggsignif)
library(emmeans)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

se <- function(x){
  sd(x, na.rm = TRUE)/sqrt(length(unique(x)))
}

#### import data ####
## import temporal beta diversity df 
gain.loss <- read_rds("output/bcd matrix.rds")

## import 2023 richness data
no.spp <- read_rds("output/richness23.rds")
str(no.spp)

## import treatment data 
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X) %>%
  merge(gain.loss, by = "Plot") %>% #merge with temporal beta diversity df
  merge(no.spp, by = "Plot") %>% # merge with richness values
  filter(Thin == "Not thinned") #use only not thinned plots b/c not the focus of MS
str(trts)

hist(trts$spp.losses)

#### models ####
spp.gains <- lmer(spp.gains ~ Nutrients*Light + (1|Block), data = trts)
summary(spp.gains)
Anova(spp.gains)
qqnorm(resid(spp.gains))
qqline(resid(spp.gains))
shapiro.test(resid(spp.gains))

spp.losses <- lmer(spp.losses ~ Nutrients*Light + (1|Block), data = trts)
summary(spp.losses)
Anova(spp.losses)
emmeans(spp.losses, pairwise ~ Light, adjust = "Tukey")
qqnorm(resid(spp.losses))
qqline(resid(spp.losses))
shapiro.test(resid(spp.losses))

spp.turnover.resources <- lmer(turnover ~ Nutrients*Light + (1|Block), data = trts)
summary(spp.turnover.resources)
Anova(spp.turnover.resources)
qqnorm(resid(spp.turnover.resources))
qqline(resid(spp.turnover.resources))
shapiro.test(resid(spp.turnover.resources))

biomass.resources <- lmer(AB.biomass ~ Nutrients*Light + (1|Block), data = trts)
summary(biomass.resources)
Anova(biomass.resources)
qqnorm(resid(productivity.resources))
qqline(resid(productivity.resources))
shapiro.test(resid(biomass.resources))
emmeans(biomass.resources, pairwise ~ Light, adjust = "Tukey") #difference appears to be between
                                                                  #ambient and shade plots

## mean light values
light.data.summary <- trts %>%
  pivot_longer(cols = c(spp.gains, spp.losses, turnover, AB.biomass), names_to = "index", values_to = "index.value") %>%
  group_by(Light, index) %>%
  summarise(mean = mean(index.value,na.rm=TRUE),
            se = se(index.value))
light.data.summary

## biomass values for resource combination
biomass.resource.summary <- trts %>%
  group_by(Nutrients, Light) %>%
  summarise(mean.biomass = mean(AB.biomass))
biomass.resource.summary
