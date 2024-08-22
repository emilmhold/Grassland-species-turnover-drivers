## Which species are gained and lost?
## Author: Emily H
## Created: July 9. 2024
## Last edited: August 20, 2024

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("lme4")
#install.packages("lmerTest")
#install.packages("car")
#install.packages("DHARMa")
#install.packages("cowplot")

library(vegan)
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(DHARMa)
library(cowplot)
library(glmmTMB)
library(sjPlot)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

options(scipen = 999)

#### import data ####
## import July 2023 data
kin23 <- read.csv("data/July 2023 % cover.csv") %>%
  mutate_at(c(7:84), as.numeric) %>%
  dplyr::select(!(c(Date.Collected, Waypoint, Block, Plot.1, Treatment, Bare.Ground, Litter, Cow.pat, Rock))) %>%
  dplyr::select(!(starts_with("Unk"))) %>%
  rename(Rub.ida2 = Arb.ida) %>%
  mutate(Rub.ida1 = Rub.ida + Rub.ida2) %>% #collapse double observations
  dplyr::select(!c(Rub.ida, Rub.ida2)) %>%
  rename(Lat.ven = Other.lathyrus,
         Rub.ida = Rub.ida1,
         Cir.und = Cir.umb) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  column_to_rownames("Plot") %>%
  dplyr::select(order(colnames(.))) %>% #put columns in alphabetical order
  decostand(method = "pa") %>% # convert data to presence/absence
  rownames_to_column(var = "Plot")
str(kin23)

## import July 2021 data
##calculate total abundance for later
kin21 <- read.csv("data/July 2021 % cover.csv")
kin21[kin21 == "t"] <- '1'
kin21 <- kin21 %>%
  mutate_at(c(6:76), as.numeric) %>%
  dplyr::select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  dplyr::select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  dplyr::select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  column_to_rownames("Plot") %>%
  dplyr::select(.,order(colnames(.)))
kin21_total.abundance <- kin21 %>%
  summarize_if(is.numeric, sum, na.rm=TRUE) %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(Total.abundance = V1) %>%
  rownames_to_column(var = "Spp")
str(kin21_total.abundance)

kin21 <- read.csv("data/July 2021 % cover.csv")
kin21[kin21 == "t"] <- '1'
kin21 <- kin21 %>%
  mutate_at(c(6:76), as.numeric) %>%
  dplyr::select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  dplyr::select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  dplyr::select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  column_to_rownames("Plot") %>%
  dplyr::select(.,order(colnames(.))) %>% #put columns in alphabetical order
  decostand(method = "pa") %>% # convert data to presence/absence
  rownames_to_column(var = "Plot")
str(kin21)

kin21_total.abundance <- read.csv("data/July 2021 % cover.csv")
kin21_total.abundance[kin21_total.abundance == "t"] <- '1'
kin21_total.abundance <- kin21_total.abundance %>%
  mutate_at(c(6:76), as.numeric) %>%
  dplyr::select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  dplyr::select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  dplyr::select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  column_to_rownames("Plot") %>%
  dplyr::select(.,order(colnames(.)))

kin21_total.abundance <- colSums(kin21_total.abundance) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Spp") %>%
  rename(Total.abundance = ".")
str(kin21_total.abundance)

## import treatment data 
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X) %>%
  mutate(across(where(is.integer), as.character))
str(trts)

## import trait database spp average traits
traits <- read_rds("data/Cahill trait database species averages.rds") %>% # trait data: species are rows and traits are columns
  rownames_to_column(var = "Species") %>%
  dplyr::select(Species, Max.Height, SLA.Leaf, RTD) %>%
  dplyr::rename(Spp = Species,
                SLA = SLA.Leaf)
traits$Spp <- paste0(substr(traits$Spp, start = 1, stop = 3),
                         ".",
                         substr(traits$Spp, start = 4, stop = 6)) #add "." into spp codes
str(traits)

## import social metric data 
cooccur.count <- read_rds("output/species social metrics.rds") %>%
  rownames_to_column(var = "Spp") %>% 
  dplyr::mutate(across(where(is.integer),as.numeric))
str(cooccur.count)

#### merge p/a dfs ####
## convert dfs to long format
kin23_long <- kin23 %>%
  pivot_longer(cols = 2:72, names_to = "Spp", values_to = "Pres_23")
str(kin23_long)

kin21_long <- kin21 %>%
  pivot_longer(cols = 2:61, names_to = "Spp", values_to = "Pres_21")
str(kin21_long)

#merge dfs
pa_long <- kin21_long %>%
  full_join(kin23_long, by = c("Plot","Spp")) %>%
  full_join(trts, by = "Plot") %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  mutate(gained = if_else(Pres_21 == 0 & Pres_23 == 1, 1, 0)) %>% #code spp gains as 1
  mutate(lost = if_else(Pres_21 == 1 & Pres_23 == 0, 1, 0)) %>% # codes spp loss as 1
  filter(Thin == "Not thinned") #drop thinned plots
str(pa_long)

#### count gains and losses by spp ####
spp_summary <- pa_long %>%
  group_by(Spp) %>%
  summarise(present_21 = sum(Pres_21, na.rm = TRUE),
            present_23 = sum(Pres_23, na.rm = TRUE),
            gained = sum(gained, na.rm = TRUE),
            lost = sum(lost, na.rm = TRUE)) %>%
  left_join(traits, by = "Spp") %>% ## add height data
  dplyr::filter(present_21 > 24 | present_23 > 24) %>% #filter data to remove spp found in less than 10% of plots
  dplyr::mutate(prop.gained = gained/present_21,
                prop.lost = lost/present_21) %>% #calculate proportion of times a species was present 
                                                  #in a plot that it went missing/not
  left_join(kin21_total.abundance, by = "Spp") %>% # add total abundance data
  left_join(cooccur.count, by = "Spp") # add social metric data
str(spp_summary)
hist(spp_summary$prop.gained)

## export counts of species gain and loss for table S3
tableS3 <- spp_summary %>%
  dplyr::select(Spp, present_21, present_23, gained, lost) %>%
  dplyr::rename(Species = Spp,
                "Count of plots where present in 2021" = present_21,
                "Count of plots where present in 2023" = present_23,
                "Number of times gained" = gained,
                "Number of times lost" = lost)
write_csv(tableS3, "output/species gain and loss summary table.csv")

#### models ####
prop.loss.mod <- lm(prop.lost~Total.abundance, data = spp_summary)
qqnorm(resid(prop.loss.mod))
qqline(resid(prop.loss.mod))
shapiro.test(resid(prop.loss.mod))
Anova(prop.loss.mod)
summary(prop.loss.mod)

prop.gained.mod <- glm(prop.gained~Total.abundance, data = spp_summary, family = Gamma)
qqnorm(resid(prop.gained.mod))
qqline(resid(prop.gained.mod))
shapiro.test(resid(prop.gained.mod))
Anova(prop.gained.mod)
summary(prop.gained.mod)

traits.loss.mod <- lm(prop.lost~Max.Height*SLA*RTD, data = spp_summary)
qqnorm(resid(traits.loss.mod))
qqline(resid(traits.loss.mod))
shapiro.test(resid(traits.loss.mod))
Anova(traits.loss.mod)
summary(traits.loss.mod)

traits.gained.mod <- lm(prop.gained~Max.Height*SLA*RTD, data = spp_summary)
qqnorm(resid(traits.gained.mod))
qqline(resid(traits.gained.mod))
shapiro.test(resid(traits.gained.mod))
Anova(traits.gained.mod)
summary(traits.gained.mod)

social.loss.mod <- lm(prop.lost~negative.cooccurrences*positive.cooccurrences, data = spp_summary)
qqnorm(resid(social.loss.mod))
qqline(resid(social.loss.mod))
shapiro.test(resid(social.loss.mod))
Anova(social.loss.mod)
summary(social.loss.mod)

social.gained.mod <- glm(prop.gained~negative.cooccurrences*positive.cooccurrences, data = spp_summary, family = Gamma)
qqnorm(resid(social.gained.mod))
qqline(resid(social.gained.mod))
shapiro.test(resid(social.gained.mod))
Anova(social.gained.mod)
summary(social.gained.mod)

#### figures ####
loss.abundance.plot <- ggplot(spp_summary, aes(x = Total.abundance, y = prop.lost)) +
  geom_point() + 
  labs(y = "Proportion of times lost when \nspecies initially present",
       x = "Total % cover in 2021") +
  geom_smooth(method = "lm", se = FALSE) +
  ylim(0,0.8) +
  theme_classic(base_size = 20)
loss.abundance.plot

gains.abundance.plot <- ggplot(spp_summary, aes(x = Total.abundance, y = prop.gained)) +
  geom_point() + 
  labs(y = "Proportion of times gained when \nspecies not initially present",
       x = "Total % cover in 2021") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 20)
gains.abundance.plot

losses.height.plot <- ggplot(spp_summary, aes(x = Max.Height, y = prop.lost)) +
  geom_point() + 
  labs(y = "Proportion of times lost when \nspecies initially present",
       x = "Species mean maximum height (cm)") +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 20)
losses.height.plot

gains.height.plot <- ggplot(spp_summary, aes(x = Max.Height, y = prop.gained)) +
  geom_point() + 
  labs(y = "Proportion of times gained when \nspecies not initially present",
       x = "Species mean maximum height (cm)") +
  theme_classic(base_size = 20)
gains.height.plot

##put together
turnover.abundance.plot <- cowplot::plot_grid(loss.abundance.plot, gains.abundance.plot,
                                              losses.height.plot, gains.height.plot,
                                              nrow = 2,
                                              labels = "auto")
turnover.abundance.plot
ggsave(filename = "turnover abundance plots.png", 
       turnover.abundance.plot,
       path = "figures/",
       width = 20,
       height = 10,
       units = "in"
)