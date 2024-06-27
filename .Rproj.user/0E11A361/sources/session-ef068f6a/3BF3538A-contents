################################################################################
################################################################################
#########################   GRIN - Disturbance    ##############################
#########################    Functional Groups    ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggsignif")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(extrafont)
#font_import()
loadfonts(device = "win")
library(tidyverse)
library(vegan)
library(agricolae)
library(ggsignif)

##########################     Read in 2021-2023 Data       #########################
GRIN = read.csv("Data/GRIN_FUN-2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Remove Year 3 - Introduction of Fire # 
GRIN = filter(GRIN, Year != 3)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
  grepl(0, Coverage) ~ 0,
  grepl(1, Coverage) ~ 0.1,
  grepl(2, Coverage) ~ 0.5,
  grepl(3, Coverage) ~ 1.5,
  grepl(4, Coverage) ~ 3.5,
  grepl(5, Coverage) ~ 7.5,
  grepl(6, Coverage) ~ 17.5,
  grepl(7, Coverage) ~ 37.5,
  grepl(8, Coverage) ~ 62.5,
  grepl(9, Coverage) ~ 85,
  grepl(10, Coverage) ~ 97.5
))

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

# Filter out years 2021 - 2023 #

# Composition of Treatments using mean coverage per species. 
df = group_by(GRIN, Year, Treatment, Group) %>% 
  dplyr::summarize(sum = sum(Coverage))
df = filter(df, sum > 1) # Species <5% Coverage not included. 

# Label Years for ggplot #
year_names <- as_labeller(
  c(`1` = "2021", `2` = "2022"))

Veg_Bar = 
  ggplot(df, aes(x = Treatment, y = sum, fill = Group)) +
  geom_col(position = "fill", color = "black", alpha = 0.5) +
  facet_wrap(vars(Year), labeller = year_names) +
  #geom_text(aes(label = Species), stat = "identity", 
  #          size = 5.25, position=position_fill(0.5), colour = "black") +
  scale_fill_manual(breaks=c('Bare', 'Forb', 'Grass', 'Sedge', 'Woody'),
                    values = c("tan", "#660066", "#339966", "#E7B800", "brown"), 
                    labels=c("Bare", "Forb", "Grass", "Sedge", 'Woody')) +
  geom_signif(y_position = c(1.01,1.01), xmin = c(0.7,2), xmax = c(1.9,3),
              annotation=c("***", "***"), size = 0.8, textsize = 15, 
              tip_length = 0.00001) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 1),
        text=element_text(size=25),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),    
        axis.title.y = element_text(size=25, face="bold", colour = "black"),   
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        strip.text.x=element_text(size = 25,colour = "black",face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +   
  labs(x = "Treatment", y = "Total Proportional Coverage")
Veg_Bar
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2021-2022_CompBar2.png", 
       width = 12, height = 8)
