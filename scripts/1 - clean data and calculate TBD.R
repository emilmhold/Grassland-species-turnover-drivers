## Calculate species plot-wise gains and losses
## Author: Emily H
## Created: November 23, 2023
## Last edited: August 18, 2024

#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("adespatial")
#install.packages("XML")
#install.packages("arsenal")

library(vegan)
library(adespatial)
library(arsenal)
library(tidyverse)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

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
  dplyr::select(order(colnames(.))) #put columns in alphabetical order
str(kin23)

#import biomass data
biomass23 <- read.csv("data/August 2023 aboveground biomass & litter weights.csv") %>%
  rename(AB.biomass = Aboveground.biomass..g.) %>%
  mutate(Productivity = AB.biomass/(0.2*0.5)) %>% #divide biomass values by surface area of clipping quadrat (25x50 cm)
  select(Plot, AB.biomass, Productivity)
  #filter(Plot != 73) #remove outlier to improve normality of residuals
str(biomass23)

# calculate richness for later
richness23 <- kin23 %>%
  mutate(richness = apply(. > 0, 1, sum)) %>%
  rownames_to_column(var = "Plot") %>%
  dplyr::select(Plot, richness) %>%
  merge(biomass23, by = "Plot") %>% 
  merge(trts, by = "Plot") %>%
  filter(Thin != "Thinned")
str(richness23)
write_rds(richness23, "output/richness23.rds")

#import July 2021 data
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
  dplyr::select(.,order(colnames(.))) #put columns in alphabetical order
str(kin21)

#import treatment data
trts <- read.csv("data/Plot treatments.csv") %>%
  dplyr::select(!X)
str(trts)

# calculate richness
richness21 <- kin21 %>%
  mutate(richness = apply(. > 0, 1, sum)) %>%
  rownames_to_column(var = "Plot") %>%
  dplyr::select(Plot, richness) %>%
  merge(trts, by = "Plot") %>%
  filter(Thin == "Not thinned")

# find which columns are not in each dataframe
summary(comparedf(kin23, kin21))

#summarize which spp are not in each df
not_in_kin21 <- c("And.sep", "Ant.neg","Ara.div", "Ara.hir", "Ast.fle", "Bro.ine", "Che.alb", "Cir.und", "Col.lin", 
                  "Dodecatheon.sp.", "Gen.ama", "Lat.ven", "Pot.con", "Sil.lat", "Son.arv", "Vio.ped", "Ziz.apt")
not_in_kin23 <- c("All.tex", "Cal.mon", "Fes.rub","Pot.div", "Spi.bet", "Sti.com")
#create new columns for these spp
kin23[not_in_kin23] <- 0
kin23 <- kin23 %>% select(order(colnames(.))) #put columns back in alphabetical order
kin21[not_in_kin21] <- 0
kin21 <- kin21 %>% select(order(colnames(.)))

#check that I haven't missed any columns
summary(comparedf(kin23, kin21))

#### calculate temporal beta diversity (Legendre 2019) to quantify spp turnover ####
tbi <- TBI(kin21, kin23, method = "jaccard", pa.tr = TRUE, nperm = 999, BCD = TRUE)
bcd.mat <- tbi$BCD.mat %>%
  mutate(Plot = 1:204, .before = 1) %>%
  rename(spp.losses = `B/(A+B+C)`,
         spp.gains = `C/(A+B+C)`,
         turnover = `D=(B+C)/(A+B+C)`)
## When B > C, the site has lost species or abundances-per-species between time 1 and time 2; 
#this is indicated by a "-" sign in column Change. On the contrary, if B < C, the site has gained species or 
#abundances-per-species between time 1 and time 2; this is indicated by a "+" sign in that column. 
#Sites with equal amounts of losses and gains are marked with a "0".

## look @ distribution of variables ###
hist(bcd.mat$spp.losses)
bcd.mat$ln.losses <- log(bcd.mat$spp.losses)
hist(bcd.mat$ln.losses)

hist(bcd.mat$spp.gains)
bcd.mat$ln.gains <- log(bcd.mat$spp.gains)
hist(bcd.mat$ln.gains)

bcd.mat <- bcd.mat %>% mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 0, .x))) #get rid of infinite values
str(bcd.mat)

#### Table S4 ####
bcd.mat.summary <- bcd.mat %>%
  left_join(trts, by = "Plot") %>% 
  filter(Thin == "Not thinned") %>%
  group_by(Change) %>%
  count()
bcd.mat.summary

#### export files ####
write_rds(bcd.mat, "output/bcd matrix.rds")
bcd.summary <- tbi$BCD.summary 
bcd.summary
t.test.bc <- tbi$t.test_B.C
