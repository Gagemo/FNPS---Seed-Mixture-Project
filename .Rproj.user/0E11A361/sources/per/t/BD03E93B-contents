################################################################################
################################################################################
#########################   FNPS - Seed Mixture   ##############################
#########################    NMDS - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2022 - 2023         ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggsignif", "multcompView", "ggpubr", "rstatix",
                      "vegan", "labdsv", "pairwiseAdonis", "devtools")
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
library(multcompView)
library(ggpubr)
library(rstatix)
library(vegan)
library(labdsv)
library(devtools)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

##########################     Read in 2022-2023 Data       ####################

Data = read.csv("Data/FNPS - Seed Mixture Project - 2021-2023.csv")
Data$Coverage = as.numeric(Data$Coverage)

str(Data)
summary(Data)

Data_22 = filter(Data, Year==1)
Data_23 = filter(Data, Year==2)

########################### 2022 Data ##########################################

# Create species pivot table #
Spp_22 = dplyr::select(Data_22, ID, Species, Coverage) %>% matrify()

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(Data_22, ID, Treatment, Sub_Plot) %>% summarise()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot three dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS_22 = metaMDS(Spp_22, distance = "bray", k=4)
MDS_22$stress
stressplot(MDS_22) 
goodness(MDS_22)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(Data_22, Species, Group) %>% summarise()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- 
  as.data.frame(vegan:::scores.metaMDS(MDS_22, display = c("species")))

# create a column of species, from the row names of species.scores  #   
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #                                                      
species.scores$Group <- species_groups$Group

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS_22 = data.frame(ID = Treat$ID, MDS_22 = MDS_22$points, Treatment = Treat$Treatment,
                  Sub_Plot = Treat$Sub_Plot)

# NMDS Graphs
NMDS_graph_22 = 
  ggplot() +
  geom_point(data = NMDS_22, aes(x = MDS_22.MDS1, y = MDS_22.MDS2, fill = Treatment),
             alpha = 0.7, size = 5, shape = 21) +
  ylim(-0.8,0.8) +
  # geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species)) +
  annotate("text", x = -1, y = 0.5, 
           label = paste0("Stress: ", format(MDS_22$stress, digits = 2)), 
           hjust = 0, size = 8) +
  ggtitle("2022") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black", 
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_graph_22

ggsave("Figures/NMDS_22.PNG", 
       width = 10, height = 7)

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Spp_22 ~ Treatment, data = NMDS_22, method="bray")
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp_22 ~ Treatment, data = NMDS_22)
pairwise.adonis

##########################     2023 Data       #################################

# Create species pivot table #
Spp_23 = dplyr::select(Data_23, ID, Species, Coverage) %>% matrify()

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(Data_23, ID, Treatment) %>% summarise()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot two dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS_23 = metaMDS(Spp_23, distance = "bray", trymax = 500, maxit = 999, k=3, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS_23$stress
stressplot(MDS_23) 
goodness(MDS_23)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(Data_23, Species, Group) %>% summarise()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- 
  as.data.frame(vegan:::scores.metaMDS(MDS_23, display = c("species")))

# create a column of species, from the row names of species.scores  #                                                              
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #
species.scores$Group <- species_groups$Group

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS_23 = data.frame(MDS_23 = MDS_23$points, Treatment = Treat$Treatment, 
                  Plot = Treat$ID)

# NMDS Graphs
NMDS_graph_23 = 
  ggplot() +
  geom_point(data = NMDS_23, aes(x = MDS_23.MDS1, y = MDS_23.MDS2, fill = Treatment),
             alpha = 0.7, size = 5, shape = 21) +
  ylim(-1,1) +
  # geom_text(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, label = Plot)) +
  annotate("text", x = -1, y = 0.5, 
           label = paste0("Stress: ", format(MDS_23$stress, digits = 2)), 
           hjust = 0, size = 8) +
  ggtitle("2023") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black",
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),    
        axis.title.y = element_text(size=25, face="bold", colour = "black"),   
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Treatment", 
       fill = "Treatment")
NMDS_graph_23

ggsave("Figures/NMDS_23.png", 
       width = 10, height = 7)

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Spp_23 ~ NMDS_23$Treatment, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp_23 ~ Treatment, data = NMDS_23)
pairwise.adonis

################## Save Figures Above using ggarrange ##########################
NMDS_22_23 = 
  ggarrange(NMDS_graph_22, NMDS_graph_23, ncol = 2, nrow = 1, 
            common.legend = TRUE, legend="bottom")
NMDS_22_23

ggsave("Figures/22-23_NMDS.png", 
       width = 18, height = 10)

