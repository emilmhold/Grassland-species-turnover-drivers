## model spp cooccurrences against resource gradients
## Author: Emily H
## Created: January 30, 2024
## Last edited: July 11, 2024

#install.packages("cooccur")
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("patchwork")
#install.packages("cowplot")

library(cooccur)
library(tidyverse)
library(vegan)
library(patchwork)
library(cowplot)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

#### import data ####
kin <- read.csv("data/July 2021 % cover.csv")
kin[kin == "t"] <- '1'
kin <- kin %>%
  mutate_at(c(6:76), as.numeric) %>%
  select(!c(Date.Collected, Block, Waypoint, Treatment, Bare.Ground, Litter, moss, Cow.pat, X)) %>%
  rename(Hel.hoo2 = Heli.hoo) %>%
  mutate(Hel.hoo1 = Hel.hoo + Hel.hoo2,
         Ery.che1 = Ery.che + rough.mustard) %>% #collapse double observations
  select(!c(Hel.hoo, Hel.hoo2, fuzz.auricle.grass, slender.wheat.grass, Ery.che, rough.mustard)) %>%
  select(!(starts_with("un"))) %>%
  rename(Ago.gla = Tradub.like,
         Hel.hoo = Hel.hoo1,
         Ery.che = Ery.che1,
         Pas.smi = Pas.mi) %>%
  mutate_all(~replace(., is.na(.), 0)) %>% #replace empty cells with 0
  select(.,order(colnames(.))) #put columns in alphabetical order
str(kin)

## import N data 
N21 <- read.csv("data/2021 Soil N values.csv") %>%
  mutate(NH4.change = August.NH4.g - May.NH4.g) %>% #calculate change over season
  mutate(NO3.change = August.NO3.g - May.NO3.g) %>% #calcualte change over season
  mutate(TON = August.NH4.g + August.NO3.g) %>%
  dplyr::select(Plot, NH4.change, NO3.change, TON)
str(N21)

## import PAR data
PAR21.new <- read_rds("output/PAR21 new.rds")
str(PAR21.new)

#### create cooccurrence matrices ####
#prep data
kin.net <- trts %>%
  merge(kin, by = "Plot") %>%
  filter(Light == "Ambient" & Nutrients == "No fert" & Thin == "Not thinned") %>% #select control plots
  select(!c(Block, Light, Nutrients, Thin)) %>% #remove metadata
  column_to_rownames(var = "Plot") %>%
  t() %>% 
  decostand(method = "pa") #convert to presence/absence
str(kin.net)

## make cooccurrence matrix 
cooccur.kin <- cooccur(kin.net, type = "spp_site", thresh = TRUE, spp_names = TRUE) #note: thresh = TRUE means spp pairs 
#expected to have less than one co-occurrence are filtered from the dataset
summary(cooccur.kin)
plot(cooccur.kin) 
result.kin <- prob.table(cooccur.kin) %>%
  select(sp1_name, sp2_name, prob_cooccur)
str(result.kin)

#### model species abundances by resource gradients ####
## Note: only modelling for species pairs predicted to cooccur >1 time
kin <- kin %>% column_to_rownames(var = "Plot")

## convert abundance data to proportional cover
prop.kin <- kin / rowSums(kin)

df <- data.frame(Spp = colnames(prop.kin),
                 Cooccurrences = colnames(prop.kin) %in% result.kin$sp1_name)
drops <- as.character(subset(df, Cooccurrences=="FALSE")$Spp) # There is no trait data for these taxa.
df[df$Spp=="Vic.ame", "Cooccurrences"] <- TRUE #I have cooccurrence data for Vic ame
table(df$Cooccurrences) #24 spp are expected to cooccur >1 time
prop.kin<-prop.kin[ , !(names(prop.kin) %in% drops)] #remove spp without >1 significant cooccurrences

resources.cover <- prop.kin %>% 
  #rename_with(~ paste(., "cover", sep = "_")) %>% #distinguish cover columns
  mutate(Plot = 1:204) %>%
  mutate_if(is.character, as.numeric) %>%
  full_join(N21, by = "Plot") %>%
  #select(!NO3.change) %>%
  full_join(PAR21.new, by = "Plot") %>%
  pivot_longer(1:23, names_to = "spp", values_to ="prop.cover")
str(resources.cover)

#### plot ####
p1 <- ggplot(resources.cover, aes(x=NH4.change, y=prop.cover, colour = spp)) +
  geom_smooth(linewidth=1.5, se = FALSE) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Ammonium change (mg/g)",
       y = "Proportional cover",
       color = "Species") + 
  theme_classic(base_size = 18)
p1

p2 <- ggplot(resources.cover, aes(x=True.light.penetration, y=prop.cover, colour = spp)) +
  geom_smooth(linewidth=1.5, se = FALSE) +
  labs(x = "Light penetration (PAR)",
       y = " ",
       color = "Species") + 
  theme_classic(base_size = 18)
p2

p3 <- ggplot(resources.cover, aes(x=TON, y=prop.cover, colour = spp)) +
  geom_smooth(linewidth=1.5, se = FALSE) +
  labs(x = "Total organic nitrogen (mg/g)",
       y = "Proportional cover",
       color = "Species") + 
  theme_classic(base_size = 18)
p3

p4 <- ggplot(resources.cover, aes(x=NO3.change, y=prop.cover, colour = spp)) +
  geom_smooth(linewidth=1.5, se = FALSE) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "Nitrate change (mg/g)",
       y = " ",
       color = "Species") + 
  theme_classic(base_size = 18)
p4

#put panels together
niche.plots <- cowplot::plot_grid(
  p1 + theme(legend.position = "none"), 
  p4 + theme(legend.position = "none"),
  #p3 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  nrow = 1,
  labels = 'auto')
#extract legend
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)
#add legend to plot
niche.plots2 <- cowplot::plot_grid(niche.plots, legend, rel_widths = c(5, 1))
niche.plots2

ggsave(filename = "figures/niche plots.png", 
       niche.plots2,
       width = 16,
       height = 8,
       units = "in"
)

#### rda ####
##create abiotic matrix
cor(N21[1:4]) #N measures are correlated
resource.df <- N21 %>%
  full_join(PAR21.new, by = "Plot") %>% 
  column_to_rownames(var = "Plot") %>%
  dplyr::select(NH4.change, True.light.penetration) %>% 
  dplyr::rename(Light.penetration = True.light.penetration)
str(resource.df)
cor(resource.df) #not correlated

kin.rda <- rda(prop.kin ~ ., data = resource.df)
summary(kin.rda) #model explains only 4.6% of variation

#open figure
png("figures/species in resource space.png", width = 600, height = 400)
sc_sp <- scores(kin.rda, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(kin.rda, display="bp", choices=c(1, 2), scaling=1)
plot(kin.rda,  
     type = "points",
     col = "grey70",
     xlim = c(-2,2), 
     ylim = c(-2,2))
text(sc_sp + c(0.2, 0.1), # adjust text coordinates to avoid overlap with points 
     labels = rownames(sc_sp), 
     col = "black", 
     cex = 0.6)
arrows(0,0, # start them from (0,0)
       sc_bp[,1], sc_bp[,2], # end them at the score value
       col = "blue", 
       lwd = 1)
text(x = sc_bp[,1] -0.1, # adjust text coordinate to avoid overlap with arrow tip
     y = sc_bp[,2] - 0.03, 
     labels = rownames(sc_bp), 
     col = "red", 
     cex = 1, 
     font = 2)
#close the file
dev.off()

anova.cca(kin.rda, step = 1000) ## p = 0.001
anova.cca(kin.rda, step = 1000, by = "axis") ## RDA one is sig
anova.cca(kin.rda, step = 1000, by = "term") ## light seems to be driving the significance

##permanova
adonis2(prop.kin ~ NH4.change + True.light.penetration, data = resource.df, permutations = 999)


ordiplot(kin.rda, scaling = 1, type = "text", display = "species")
#linestack(scores(kin.rda, display = "sites"), cex = 1.1,)
#linestack(scores(kin.rda, display = "species"), cex = 1.1)
