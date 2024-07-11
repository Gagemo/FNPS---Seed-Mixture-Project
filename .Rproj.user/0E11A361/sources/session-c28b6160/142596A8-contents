################################################################################
################################################################################
#########################   GRIN - Disturbance    ##############################
#########################    HEAT - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2021 - 2023         ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "labdsv", "pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)

##########################     Read in  Data       #########################
GRIN = read.csv("Data/GRIN - 2020-2023.csv")

GRIN$Coverage = as.numeric(GRIN$Coverage)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
  grepl(10, Coverage) ~ 97.5,
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
  
))

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

GRIN$YID <- paste(GRIN$Year,GRIN$ID)

GRIN$ID_ <- paste(GRIN$Treatment, GRIN$ID)

# Seperate Pre & Post-Treatment Data #
GRIN_21 = filter(GRIN, Year == "1")
GRIN_22 = filter(GRIN, Year == "2")

# Orders years and treatments so that they display in same sequence in graphs #
GRIN$Year = factor(GRIN$Year, levels=c('1','2'))
GRIN$Treatment = factor(GRIN$Treatment, levels=c('C','Tsp','Tw'))

#Renames values in fire treatments for heat map later #
GRIN$Treatment <- recode(GRIN$Treatment, 
                         C ="No-Till", Tw = "Winter", Tsp = "Late-Spring")

# Create Species Pivot Table with revisions to analyses #
# All treatments #
Spp <- dplyr::select(GRIN, YID, Species, Coverage) %>% 
  matrify()

# Select 2021 Data matrix #
Spp_21 <- dplyr::select(GRIN_21, YID, Species, Coverage) %>% 
  matrify()

# Select 2022 Treatment Data matrix #
Spp_22 <- dplyr::select(GRIN_22, YID, Species, Coverage) %>% 
  matrify()

#################### Species abundances ########################################
# Creates and joins  data year 22 & 23 to make long data format #
Two_Abundance <- GRIN[which(GRIN$Year == "1"),]
Three_Abundance <- GRIN[which(GRIN$Year == "2"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance, 
                         by = c('ID_', "Treatment", 'Species'))
Abundance_w = arrange(Abundance_w, Treatment)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(ID_, Treatment, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(ID_, Treatment, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x) %>%
  filter(Change_abundance != 0)

Change_Abundance_H = filter(Change_Abundance, Change_abundance >= 25)
Change_Abundance_L = filter(Change_Abundance, Change_abundance <= -25)

Change_Abundance = full_join(Change_Abundance_H, Change_Abundance_L)

Treat = ungroup(Change_Abundance) %>% 
  dplyr::select(ID_, Treatment) %>%
  group_by(Treatment, ID_) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID_')


GRIN_change <- ungroup(Change_Abundance) %>%
  dplyr::select(ID_, Species, Change_abundance) %>%
  as.data.frame() %>%
  matrify()

GRIN_change = as.matrix(GRIN_change[, -1])

ann_colors = list(
  Treatment = c("No-Till" = "#FF99CC", "Late-Spring" = "#FFFF99", 
                "Winter" = "#99CCFF"))

speciesHEAT = pheatmap(GRIN_change, show_rownames=F, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat, 
                       annotation_colors = ann_colors, fontsize = 20,
                       border_color = "black", display_numbers = FALSE,
                       cellheight=10, cellwidth = 20)
speciesHEAT

png(file = "Figures/Chapter 1 - Soil Disturbance Seasonality/Heat.png", 
    units="cm", width=20, height=30, res=1500)
speciesHEAT
dev.off()

