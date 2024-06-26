################################################################################
################################################################################
#########################   FNPS - Seed Mixture   ##############################
#########################    HEAT - Community     ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2022 - 2023        ##############################
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
Data = read.csv("Data/FNPS - Seed Mixture Project - 2021-2023.csv")

Data$Coverage = as.numeric(Data$Coverage)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
Data <- mutate(Data, Coverage = case_when(
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

str(Data)
summary(Data)

Data$YID <- paste(Data$Year,Data$ID)

Data$ID_ <- paste(Data$Treatment, Data$ID)

# Separate Pre & Post-Treatment Data #
Data_22 = filter(Data, Year == "1")
Data_23 = filter(Data, Year == "2")

# Orders years and treatments so that they display in same sequence in graphs #
Data$Year = factor(Data$Year, levels=c('1','2'))
Data$Treatment = factor(Data$Treatment, levels=c('C','W', 'BH','BM', 'LH', 'LM'))

#Renames values in Treatment treatments for heat map later #
Data$Treatment <- recode(Data$Treatment, 
                    C ="Control", W = "Wiregrass", BH = "Broomsedge High", 
                    BM = "Broomsedge Medium", LH = "Lovegrass High", 
                    LM = "Lovegrass Medium")

# Create Species Pivot Table with revisions to analyses #
# All treatments #
Spp <- dplyr::select(Data, YID, Species, Coverage) %>% 
  matrify()

# Select 2022 Data matrix #
Spp_22 <- dplyr::select(Data_22, YID, Species, Coverage) %>% 
  matrify()

# Select 2022 Treatment Data matrix #
Spp_23 <- dplyr::select(Data_23, YID, Species, Coverage) %>% 
  matrify()

#################### Species abundances ########################################
# Creates and joins  data year 22 & 23 to make long data format #
Two_Abundance <- Data[which(Data$Year == "1"),]
Three_Abundance <- Data[which(Data$Year == "2"),]

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

Change_Abundance_H = filter(Change_Abundance, Change_abundance >= 5)
Change_Abundance_L = filter(Change_Abundance, Change_abundance <= -5)

Change_Abundance = full_join(Change_Abundance_H, Change_Abundance_L)

Treat = ungroup(Change_Abundance) %>% 
  dplyr::select(ID_, Treatment) %>%
  group_by(Treatment, ID_) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID_')


Data_change <- ungroup(Change_Abundance) %>%
  dplyr::select(ID_, Species, Change_abundance) %>%
  as.data.frame() %>%
  matrify()

Data_change = as.matrix(Data_change[, -1])

#ann_colors = list(
#  Treatment = c("No Burn" = "#333333", "Late-Spring" = "#FF9900", 
#           "Winter" = "#3366FF"))


speciesHEAT = pheatmap(Data_change, show_rownames=T, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat, fontsize = 20,
                       border_color = "black", display_numbers = FALSE,
                       cellheight=16, cellwidth = 20,
                       color=colorRampPalette(c("blue", "white", "red"))(50))
speciesHEAT

png(file = "Figures/Chapter 2 - Treatment/Heat.png", 
    units="cm", width=20, height=20, res=100)
speciesHEAT
dev.off()
