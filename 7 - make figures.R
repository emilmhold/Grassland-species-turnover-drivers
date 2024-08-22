## Make figures for paper
## Author: Emily H
## Created: June 6, 2024
## Last edited: August 13, 2024

#install.packages("tidyverse")
#install.packages("scales")
#install.packages("cowplot")

library(tidyverse)
library(scales)
library(cowplot)

setwd("~/Documents/GitHub/Analyses/Resource, traits, and social networks drive spp turnover")

#### import data ####
## import social metric file
social.metrics <- read_rds("output/CWM coocurrences.rds") %>%
  rename(gregariousness = total.pos.cooccurrence,
         isolationism = total.neg.cooccurrence) %>%
  dplyr::select(Plot, gregariousness, isolationism, spp.losses, spp.gains, turnover, AB.biomass)
str(social.metrics)

##import trait .rds file
gain.loss.new <- read_rds("output/gain.loss.rds") %>%
  dplyr::select(Plot, Block, Light, Nutrients, Thin, Max.Height, SLA, RTD) %>%
  full_join(social.metrics, by = "Plot") %>%
  dplyr::mutate(Nutrients = replace(Nutrients, Nutrients == 'No fert', 'Not fertilized'))
str(gain.loss.new)

nutrients.for.plots <- gain.loss.new %>%
  group_by(Nutrients) %>%
  summarise(mean.gains = mean(spp.gains),
            se.gains = sd(spp.gains)/sqrt(length(spp.gains)),
            mean.losses = mean(spp.losses),
            se.losses = sd(spp.losses)/sqrt(length(spp.losses)),
            mean.turnover = mean(turnover),
            se.turnover = sd(turnover)/sqrt(length(turnover)),
            mean.ab.biomass = mean(AB.biomass),
            se.ab.biomass = sd(AB.biomass)/sqrt(length(AB.biomass))
            )
str(nutrients.for.plots)
##find how many times more productive fertilized plots were
42.6/36.1

light.for.plots <- gain.loss.new %>%
  group_by(Light) %>%
  summarise(mean.gains = mean(spp.gains),
            se.gains = sd(spp.gains)/sqrt(length(spp.gains)),
            mean.losses = mean(spp.losses),
            se.losses = sd(spp.losses)/sqrt(length(spp.losses)),
            mean.turnover = mean(turnover),
            se.turnover = sd(turnover)/sqrt(length(turnover)),
            mean.ab.biomass = mean(AB.biomass),
            se.ab.biomass = sd(AB.biomass)/sqrt(length(AB.biomass))
            )
str(light.for.plots)

#### gains figures ####
height.gains.plot <- ggplot(gain.loss.new, aes(x = Max.Height, y = spp.gains)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  #geom_smooth(method = "lm", se = FALSE) +
  scale_fill_manual(values = "#1F618D") +
  theme_classic(base_size = 20) + 
  annotate("text", x = -Inf, y = Inf, label = bquote(bold("p = 0.235")), hjust = -0.1, vjust = 1.1, size = 5)
height.gains.plot

SLA.gains.plot <- ggplot(gain.loss.new, aes(x = SLA, y = spp.gains)) +
  geom_point() +
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.674")), hjust = 1.1, vjust = 1.1, size = 5)
SLA.gains.plot

RTD.gains.plot <- ggplot(gain.loss.new, aes(x = RTD, y = spp.gains)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) +
  annotate("text", x = -Inf, y = Inf, label = bquote(bold("p = 0.258")), hjust = -0.1, vjust = 1.1, size = 5)
RTD.gains.plot

nutrients.gains.plot <- ggplot(nutrients.for.plots, aes(x = Nutrients, y = mean.gains, fill=Nutrients)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.gains-se.gains, ymax=mean.gains+se.gains), position = "dodge") + 
  labs(x = " ",
       y = "Relative contribution of \nspecies gains") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#1F618D", "#DAF7A6")) +
  theme_classic(base_size = 20) + 
  annotate("text", x = -Inf, y = Inf, label = bquote(bold("p = 0.569")), hjust = -0.1, vjust = 1.1, size = 5)
nutrients.gains.plot

light.gains.plot <- ggplot(light.for.plots, aes(x = Light, y = mean.gains, fill = Light)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.gains-se.gains, ymax=mean.gains+se.gains), position = "dodge") + 
  labs(x = " ",
       y = " ") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#FFC300", "#FF5733", "#900C3F")) +
  #geom_signif(stat="identity",aes(x=signif, xend = signif, y=y, yend = y, annotation="*"), size = 5) +
  theme_classic(base_size = 20) +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.903")), hjust = 1.1, vjust = 1.1, size = 5)
light.gains.plot

#### losses figures ####
height.losses.plot <- ggplot(gain.loss.new, aes(x = Max.Height, y = spp.losses)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.307")), hjust = 1.1, vjust = 1.1, size = 5)
height.losses.plot

SLA.losses.plot <- ggplot(gain.loss.new, aes(x = SLA, y = spp.losses)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.507")), hjust = 1.1, vjust = 1.1, size = 5)
SLA.losses.plot

RTD.losses.plot <- ggplot(gain.loss.new, aes(x = RTD, y = spp.losses)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.786")), hjust = 1.1, vjust = 1.1, size = 5)
RTD.losses.plot

nutrients.losses.plot <- ggplot(nutrients.for.plots, aes(x = Nutrients, y = mean.losses, fill = Nutrients)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.losses-se.losses, ymax=mean.losses+se.losses), position = "dodge") + 
  labs(x = " ",
       y = "Relative contribution of \nspecies losses") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#1F618D", "#DAF7A6")) +
  #geom_signif(stat="identity",aes(x=signif, xend = signif, y=y, yend = y, annotation="*"), size = 5) +
  theme_classic(base_size = 20) + 
  annotate("text", x = 1, y = nutrients.for.plots$mean.losses[1] + 0.04, label = "a", size = 6) +
  annotate("text", x = 2, y = nutrients.for.plots$mean.losses[2] + 0.04, label = "b", size = 6) +
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p < 0.001")), hjust = 1, vjust = 1.1, size = 5)
  #annotate("text", x = 1.5, y = max(nutrients.for.plots$mean.losses) + 0.03, label = "*", size = 6, fontface = "bold")
nutrients.losses.plot

light.losses.plot <- ggplot(light.for.plots, aes(x = Light, y = mean.losses, fill = Light)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.losses-se.losses, ymax=mean.losses+se.losses), position = "dodge") + 
  labs(x = " ",
       y = " ") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#FFC300", "#FF5733", "#900C3F")) +
  theme_classic(base_size = 20) + 
  annotate("text", x = -Inf, y = Inf, label = bquote(bold("p = 0.289")), hjust = -0.1, vjust = 1, size = 5)
light.losses.plot

#### turnover figures ####
height.turnover.plot <- ggplot(gain.loss.new, aes(x = Max.Height, y = turnover)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = -Inf, label = bquote(bold("p = 0.032")), hjust = 1.1, vjust = -0.5, size = 5)
height.turnover.plot

SLA.turnover.plot <- ggplot(gain.loss.new, aes(x = SLA, y = turnover)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = -Inf, label = bquote(bold("p = 0.693")), hjust = 1.1, vjust = -0.5, size = 5)
SLA.turnover.plot

RTD.turnover.plot <- ggplot(gain.loss.new, aes(x = RTD, y = turnover)) +
  geom_point() + 
  labs(y = " ",
       x = " ") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = -Inf, label = bquote(bold("p = 0.268")), hjust = 1.1, vjust = -0.5, size = 5)
RTD.turnover.plot

nutrients.turnover.plot <- ggplot(nutrients.for.plots, aes(x = Nutrients, y = mean.turnover, fill = Nutrients)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.turnover-se.turnover, ymax=mean.turnover+se.turnover), position = "dodge") + 
  labs(x = " ",
       y = "Species \nturnover") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#1F618D", "#DAF7A6")) +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.006")), hjust = 1, vjust = 1.1, size = 5) + 
  annotate("text", x = 1, y = nutrients.for.plots$mean.turnover[1] + 0.05, label = "a", size = 6) +
  annotate("text", x = 2, y = nutrients.for.plots$mean.turnover[2] + 0.05, label = "b", size = 6)
nutrients.turnover.plot

light.turnover.plot <- ggplot(light.for.plots, aes(x = Light, y = mean.turnover, fill = Light)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.turnover-se.turnover, ymax=mean.turnover+se.turnover), position = "dodge") + 
  labs(x = " ",
       y = " ") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#FFC300", "#FF5733", "#900C3F")) +
  #geom_signif(stat="identity",aes(x=signif, xend = signif, y=y, yend = y, annotation="*"), size = 5) +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.633")), hjust = 1, vjust = 1, size = 5)
light.turnover.plot

#### biomass figures ####
height.biomass.plot <- ggplot(gain.loss.new, aes(x = Max.Height, y = AB.biomass)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = ' ',
       x = "CWM Maximum height (cm)") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.025")), hjust = 1.1, vjust = 1.1, size = 5)
height.biomass.plot

SLA.biomass.plot <- ggplot(gain.loss.new, aes(x = SLA, y = AB.biomass)) +
  geom_point() + 
  labs(y = " ",
       x = "CWM SLA (cm\u00B2/g)") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.946")), hjust = 1.1, vjust = 1.1, size = 5)
SLA.biomass.plot

RTD.biomass.plot <- ggplot(gain.loss.new, aes(x = RTD, y = AB.biomass)) +
  geom_point() + 
  labs(y = " ",
       x = "CWM RTD (mg/cm\u00B2)") +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.149")), hjust = 1.1, vjust = 1.1, size = 5)
RTD.biomass.plot

nutrients.biomass.plot <- ggplot(nutrients.for.plots, aes(x = Nutrients, y = mean.ab.biomass, fill = Nutrients)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.ab.biomass-se.ab.biomass, ymax=mean.ab.biomass+se.ab.biomass), position = "dodge") + 
  labs(x = "Nutrient treatment",
       y = "Aboveground biomass \n(g)") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#1F618D", "#DAF7A6")) +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p < 0.001")), hjust = 1, vjust = 1, size = 5) + 
  annotate("text", x = 1, y = nutrients.for.plots$mean.gains[1] + 48, label = "a", size = 6) +
  annotate("text", x = 2, y = nutrients.for.plots$mean.gains[2] + 40, label = "b", size = 6)
nutrients.biomass.plot

light.biomass.plot <- ggplot(light.for.plots, aes(x = Light, y = mean.ab.biomass, fill = Light)) +
  geom_col() +
  geom_errorbar(aes(ymin=mean.ab.biomass-se.ab.biomass, ymax=mean.ab.biomass+se.ab.biomass), position = "dodge") + 
  labs(x = "Light treatment",
       y = " ") + 
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_fill_manual(values = c("#FFC300", "#FF5733", "#900C3F")) +
  theme_classic(base_size = 20) + 
  annotate("text", x = Inf, y = Inf, label = bquote(bold("p = 0.039")), hjust = 1, vjust = 1, size = 5) +
  annotate("text", x = 1, y = light.for.plots$mean.ab.biomass[1] + 10, label = "a", size = 6) +
  annotate("text", x = 2, y = light.for.plots$mean.ab.biomass[2] + 6, label = "b", size = 6) +
  annotate("text", x = 3, y = light.for.plots$mean.ab.biomass[3] + 7, label = "ab", size = 6)
light.biomass.plot

#### put it together ####
traits.resources.multiplot <- cowplot::plot_grid(
  nutrients.turnover.plot + theme(legend.position = "none"), light.turnover.plot + theme(legend.position = "none"), height.turnover.plot, SLA.turnover.plot, RTD.turnover.plot, 
  nutrients.gains.plot + theme(legend.position = "none"), light.gains.plot + theme(legend.position = "none"), height.gains.plot, SLA.gains.plot, RTD.gains.plot, 
  nutrients.losses.plot + theme(legend.position = "none"), light.losses.plot + theme(legend.position = "none"), height.losses.plot, SLA.losses.plot, RTD.losses.plot, 
  nutrients.biomass.plot + theme(legend.position = "none"), light.biomass.plot + theme(legend.position = "none"), height.biomass.plot, SLA.biomass.plot, RTD.biomass.plot, 
  nrow = 4,
  labels = 'auto')
traits.resources.multiplot

ggsave(filename = "traits and resource plots.png", 
       traits.resources.multiplot,
       path = "figures/",
       width = 30,
       height = 20,
       units = "in"
)

