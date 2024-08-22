## calculate social metrics and run models
## Author: Emily H
## Created: November 30, 2023
## Last edited: August 20, 2024

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#install.packages("cowplot")
#install.packages("sjPlot")
#install.packages("sjmisc")

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(car)
library(cowplot)
library(sjPlot)
library(sjmisc)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

#### import data ####
#import cooccurrence data
result.kin <- read_rds("output/aug 2021 kin cooccurrence matrix.rds") 
str(result.kin)

## import species richness and productivity data
no.spp <- read_rds("output/richness23.rds")
str(no.spp)

## import treatment data 
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X)
str(trts)

#import abundance data
kin <- read.csv("data/July 2021 % cover.csv")
kin[kin == "t"] <- '1'
kin <- kin %>%
  mutate_at(c(6:76), as.numeric) %>%
  dplyr::select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  dplyr::select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  #select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  column_to_rownames("Plot") %>%
  dplyr::select(.,order(colnames(.))) #put columns in alphabetical order
## convert to percent cover 
kin_percent <- kin / rowSums(kin)
str(kin_percent)

#### count how many positive vs. negative cooccurrences ####
cooccur.count <- result.kin %>%
  group_by(sp1_name, sign) %>%
  summarise(count = length(sign)) %>%
  na.omit() %>%
  group_by(sp1_name) %>%
  mutate(node.degree = sum(count)) %>% #total up positive and negative cooccurrences 
  pivot_wider(names_from = "sign", values_from = "count", values_fill = 0) %>%
  rename(spp = sp1_name,
         negative.cooccurrences = N,
         positive.cooccurrences = P) %>%
  filter(spp != "Che.alb") %>% #not in kin2 df
  column_to_rownames(var = "spp")
str(cooccur.count)
##export file
write_rds(cooccur.count, "output/species social metrics.rds")

#### calculate weighted mean of abundance * interactions ####
# find which spp are in interactions df
df <- data.frame(Spp = colnames(kin_percent),
                 interactions = colnames(kin_percent) %in% rownames(cooccur.count))
drops <- as.character(subset(df, interactions=="FALSE")$Spp) # There is no cooccurrence data for these taxa.
drops

#remove spp without cooccurrence data
kin2<-kin_percent[ , !(names(kin_percent) %in% drops)]
str(kin2) ## 23 spp - it worked!

#transpose
kin2 <- kin2 %>% t()
#put spp names in alphabetical order
kin2 <- kin2[ order(row.names(kin2)), ]
str(kin2)

## multiply %cover df by interactions for spp
weighted.degree <- as.data.frame(kin2*cooccur.count$node.degree)
weighted.degree <- weighted.degree %>% 
  t() %>%
  as.data.frame() %>%
  mutate(total.weighted.degree = rowSums(across(where(is.numeric)))) %>% #total all interactions within the plot
  rownames_to_column(var = "Plot") %>%
  dplyr::select(Plot, total.weighted.degree)

positive.weighted.cooccurrence <- kin2*cooccur.count$positive.cooccurrences
positive.weighted.cooccurrence <- positive.weighted.cooccurrence %>% 
  t() %>%
  as.data.frame() %>%
  mutate(total.pos.cooccurrence = rowSums(across(where(is.numeric)))) %>% #total all interactions within the plot
  rownames_to_column(var = "Plot") %>% 
  dplyr::select(Plot, total.pos.cooccurrence)

negative.weighted.cooccurrence <- kin2*cooccur.count$negative.cooccurrences
negative.weighted.cooccurrence <- negative.weighted.cooccurrence %>% 
  t() %>%
  as.data.frame() %>%
  mutate(total.neg.cooccurrence = rowSums(across(where(is.numeric)))) %>% #total all interactions within the plot
  rownames_to_column(var = "Plot") %>% 
  dplyr::select(Plot, total.neg.cooccurrence)

#### models ####
## import temporal beta diversity df 
gain.loss <- read_rds("output/bcd matrix.rds") %>% #import temporal beta diversity df
  merge(trts, by = "Plot") %>%#merge with treatment data
  merge(no.spp, by = "Plot") %>% ## merge with richness and productivity data
  merge(weighted.degree, by = "Plot") %>%## merge with weighted node degree data
  merge(positive.weighted.cooccurrence, by = "Plot") %>% ## merge with pos cooccurrence data
  merge(negative.weighted.cooccurrence, by = "Plot") %>% ## merge with neg cooccurrence data
  filter(Thin == "Not thinned") #use only not thinned plots b/c not the focus of MS
str(gain.loss)
## export data
write_rds(gain.loss, "output/CWM coocurrences.rds")

#### overall lmms ####
social.losses.lmm <- lmer(spp.losses ~ total.pos.cooccurrence*total.neg.cooccurrence + (1|Block), data = gain.loss)
summary(social.losses.lmm)
Anova(social.losses.lmm)
qqnorm(resid(social.losses.lmm))
qqline(resid(social.losses.lmm))
shapiro.test(resid(social.losses.lmm))

social.gains.lmm <- lmer(spp.gains ~ total.pos.cooccurrence*total.neg.cooccurrence + (1|Block), data = gain.loss)
summary(social.gains.lmm)
Anova(social.gains.lmm)
qqnorm(resid(social.gains.lmm))
qqline(resid(social.gains.lmm))
shapiro.test(resid(social.gains.lmm))

social.turnover.lmm <- lmer(turnover ~ total.pos.cooccurrence*total.neg.cooccurrence + (1|Block), data = gain.loss)
summary(social.turnover.lmm)
Anova(social.turnover.lmm)
qqnorm(resid(social.turnover.lmm))
qqline(resid(social.turnover.lmm))
shapiro.test(resid(social.turnover.lmm))

social.biomass.lmm <- lmer(AB.biomass ~ total.pos.cooccurrence*total.neg.cooccurrence + (1|Block), data = gain.loss)
summary(social.biomass.lmm)
Anova(social.biomass.lmm)
qqnorm(resid(social.biomass.lmm))
qqline(resid(social.biomass.lmm))
shapiro.test(resid(social.biomass.lmm)) 

#### figures ####
cooccurrence.forplot <- gain.loss %>%
  as.data.frame() %>%
  rename('Gregariousness' = total.pos.cooccurrence,
         'Reclusiveness' = total.neg.cooccurrence)
str(cooccurrence.forplot)
##export 
#write_rds(cooccurrence.forplot, "output/unweighted cooccurrence data for plot.rds")

#### plot ####
##losses 
greg.losses.plot <- ggplot(cooccurrence.forplot, aes(x = Gregariousness, y = spp.losses)) +
  geom_point() + 
  labs(x = " ",
       y = "Relative contribution of \nspecies losses") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.850")), hjust = 1.1, vjust = 1.1, size = 5)
greg.losses.plot

recl.losses.plot <- ggplot(cooccurrence.forplot, aes(x = Reclusiveness, y = spp.losses)) +
  geom_point() + 
  labs(x = " ",
       y = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.431")), hjust = 1.1, vjust = 1.1, size = 5)
recl.losses.plot

##gains
greg.gains.plot <- ggplot(cooccurrence.forplot, aes(x = Gregariousness, y = spp.gains)) +
  geom_point() + 
  labs(x = " ",
       y = "Relative contribution of \nspecies gains") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.048")), hjust = 1.1, vjust = 1.1, size = 5)
greg.gains.plot

recl.gains.plot <- ggplot(cooccurrence.forplot, aes(x = Reclusiveness, y = spp.gains)) +
  geom_point() + 
  labs(x = " ",
       y = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.915")), hjust = 1.1, vjust = 1.1, size = 5)
recl.gains.plot

##turnover
greg.turnover.plot <- ggplot(cooccurrence.forplot, aes(x = Gregariousness, y = turnover)) +
  geom_point() + 
  labs(x = " ",
       y = "Species \nturnover") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.850")), hjust = 1.1, vjust = 1.1, size = 5)
greg.turnover.plot

recl.turnover.plot <- ggplot(cooccurrence.forplot, aes(x = Reclusiveness, y = turnover)) +
  geom_point() + 
  labs(x = " ",
       y = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.431")), hjust = 1.1, vjust = 1.1, size = 5)
recl.turnover.plot

##productivity
greg.biomass.plot <- ggplot(cooccurrence.forplot, aes(x = Gregariousness, y = AB.biomass)) +
  geom_point() + 
  labs(x = "Gregariousness",
       y = 'Aboveground biomass \n(g)') +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.537")), hjust = 1.1, vjust = 1.1, size = 5)
greg.biomass.plot

recl.biomass.plot <- ggplot(cooccurrence.forplot, aes(x = Reclusiveness, y = AB.biomass)) +
  geom_point() + 
  labs(x = "Reclusiveness",
       y = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.806")), hjust = 1.1, vjust = 1.1, size = 5)
recl.biomass.plot

#### put it together ####
social.plots <- cowplot::plot_grid(greg.gains.plot, recl.gains.plot,
                                   greg.losses.plot, recl.losses.plot,
                                  greg.turnover.plot, recl.turnover.plot,
                                  greg.biomass.plot, recl.biomass.plot,
                                   ncol = 2,
                                  labels = "auto")
social.plots
ggsave(filename = "social metric plots.png", 
       social.plots,
       path = "figures/",
       width = 12,
       height = 15,
       units = "in"
)
