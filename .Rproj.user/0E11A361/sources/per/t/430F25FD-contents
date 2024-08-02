################################################################################
################################################################################
#########################     GRIN - Fire         ##############################
#########################    NMDS - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2021 - 2023         ##############################
################################################################################
################################################################################

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", "ggrepel",
                      "ggsignif", "multcompView", "ggpubr", "rstatix",
                      'rmarkdown', "labdsv", "pairwiseAdonis", "devtools", 
                      "knitr", "tables", "openxlsx", "labdsv")
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
library(rmarkdown)
library(ggsignif)
library(multcompView)
library(ggpubr)
library(rstatix)
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(vegan)
library(labdsv)
library(devtools)
library(knitr)
library(tables)
library(openxlsx)
library(labdsv)


install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

################################################################################
######################## Read in 2020-2023 Data ################################
################################################################################

GRIN = read.csv("Data/GRIN - 2020-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)

str(GRIN)
summary(GRIN)

# Remove Tilling Treatments # 
GRIN = filter(GRIN, Treatment == 'S')

GRIN$Treatment = factor(GRIN$Treatment, levels=c('C', 'Sp', 'W'))

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
  grepl(10, Coverage) ~ 97.5,
  grepl(1, Coverage) ~ 0.1,
  grepl(2, Coverage) ~ 0.5,
  grepl(3, Coverage) ~ 1.5,
  grepl(4, Coverage) ~ 3.5,
  grepl(5, Coverage) ~ 7.5,
  grepl(6, Coverage) ~ 17.5,
  grepl(7, Coverage) ~ 37.5,
  grepl(8, Coverage) ~ 62.5,
  grepl(9, Coverage) ~ 85,
  grepl(0, Coverage) ~ 0,
))

# Separate Years #
GRIN_20 = filter(GRIN, Year == 0)
GRIN_21 = filter(GRIN, Year == 1)
GRIN_22 = filter(GRIN, Year == 2)
GRIN_23 = filter(GRIN, Year == 3)

################################################################################
################################################################################
########################### 2023 Data ##########################################
################################################################################
################################################################################

# Create species pivot table #
Spp = dplyr::select(GRIN_23, ID, Species, Coverage) %>% matrify()

# Transforms Data #
Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_23, ID, Fire, Plot) %>% summarise()

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
MDS = metaMDS(Veg_Spp, distance = 'bray', trymax = 500, maxit = 999, k=3, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN_23, Species, Group, Native) %>% dplyr::summarize()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #                                                       
species.scores$Group <- species_groups$Group

# create a column for functional groups for NMDS graph #                                                       
species.scores$Native <- species_groups$Native

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Fire, Plot = Treat$ID)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################

envfit_result <- envfit(MDS, Spp)

# Extract vectors and/or factors
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
vectors_pvals <- as.data.frame(envfit_result$vectors$pvals)

if (!is.null(envfit_result$factors)) {
  factors_df <- as.data.frame(envfit_result$factors$centroids)
  factors_pvals <- as.data.frame(envfit_result$factors$pvals)
}
vectors_combined <- cbind(vectors_df, p_value = vectors_pvals)

write.csv(
  vectors_combined, 
  "Figures/Chapter 2 - Fire/NMDS_FitValues_23.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)

# Remove absent species
Spp <- Spp[, -absent_species]

# Perform Indicator Species Analysis
indicator_results <- indval(Spp, Treat$Fire)

# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)

# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_Tsp = indval_df$Sp,
  Indicator_Value_Tw = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)

# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species$Rank <- 
  ave(significant_species$Max_Class, significant_species$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Figures/Chapter 2 - Fire/indicator_species_analysis_results_23.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################

simper_result <- simper(Spp, Treat$Fire, permutations = 99)
simper_summary <- summary(simper_result)
simper_summary

# Extract the first pairwise comparison (adjust as needed)
pairwise_data <- simper_result[[1]]

# Create the data frame with p-values
table_data <- data.frame(
  Species = names(pairwise_data$average),
  Contribution = pairwise_data$average,
  SD = pairwise_data$sd,
  P_value = pairwise_data$p,
  Cumulative_Contribution = pairwise_data$cusum
)

# Sort by contribution
table_data <- 
  table_data[order(table_data$Contribution, decreasing = TRUE), ]

# Add rank column
table_data$Rank <- seq_len(nrow(table_data))

head(table_data)

write.csv(
  table_data, 
  file = "Figures/Chapter 2 - Fire/SIMPER_23.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################

# NMDS Graphs
NMDS_graph_23 = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = Treatment),
             alpha = 0.7, size = 5, shape = 21) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, fill = Treatment, 
                   color = Treatment), linetype = "solid", 
               show.legend = T, alpha = 0.15) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  annotate("text", x = -1, y = 1, 
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
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
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_graph_23

NMDS_Spp_graph_23 = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2), 
                  label = species.scores$species, colour = "black",
                  size = 4, fontface = "bold") +
  annotate("text", x = -1, y = 1,               
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
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
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_Spp_graph_23

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray")
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis

#save tables
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, startRow = start_row, colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis[[name]], startRow = start_row)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    # Write the character string
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis[[name]], startRow = start_row, colNames = FALSE)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/Chapter 2 - Fire/pairwise_adonis_23_same_sheet.xlsx", 
             overwrite = TRUE)

################################################################################
################################################################################
################################################################################
########################### 2022 Data ##########################################
################################################################################
################################################################################
################################################################################

# Create species pivot table #
Spp = dplyr::select(GRIN_22, ID, Species, Coverage) %>% matrify()

# Transforms Data #
Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_22, ID, Fire) %>% summarise()

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
MDS = metaMDS(Veg_Spp, distance = "bray", trymax = 500, maxit = 999, k=3, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS$stress
stressplot(MDS) 
goodness(MDS)

# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN_22, Species, Group, Native) %>% 
  dplyr::summarize()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #                                                       
species.scores$Group <- species_groups$Group

# create a column for functional groups for NMDS graph #                                                       
species.scores$Native <- species_groups$Native

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Fire, Plot = Treat$ID)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################

envfit_result <- envfit(MDS, Spp)

# Extract vectors and/or factors
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
vectors_pvals <- as.data.frame(envfit_result$vectors$pvals)

if (!is.null(envfit_result$factors)) {
  factors_df <- as.data.frame(envfit_result$factors$centroids)
  factors_pvals <- as.data.frame(envfit_result$factors$pvals)
}
vectors_combined <- cbind(vectors_df, p_value = vectors_pvals)

write.csv(
  vectors_combined, 
  "Figures/Chapter 2 - Fire/NMDS_FitValues_22.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)

# Remove absent species
Spp <- Spp[, -absent_species]

# Perform Indicator Species Analysis
indicator_results <- indval(Spp, Treat$Fire)

# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)

# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_Tsp = indval_df$Sp,
  Indicator_Value_Tw = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)

# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species$Rank <- 
  ave(significant_species$Max_Class, significant_species$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Figures/Chapter 2 - Fire/indicator_species_analysis_results_22.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################


simper_result <- simper(Spp, Treat$Fire, permutations = 99)
simper_summary <- summary(simper_result)
simper_summary

# Extract the first pairwise comparison (adjust as needed)
pairwise_data <- simper_result[[1]]

# Create the data frame with p-values
table_data <- data.frame(
  Species = names(pairwise_data$average),
  Contribution = pairwise_data$average,
  SD = pairwise_data$sd,
  P_value = pairwise_data$p,
  Cumulative_Contribution = pairwise_data$cusum
)

# Sort by contribution
table_data <- 
  table_data[order(table_data$Contribution, decreasing = TRUE), ]

# Add rank column
table_data$Rank <- seq_len(nrow(table_data))

head(table_data)

write.csv(
  table_data, 
  file = "Figures/Chapter 2 - Fire/SIMPER_22.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################

NMDS_graph_22 = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = Treatment),
             alpha = 0.7, size = 5, shape = 21) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, 
                   fill = Treatment, color = Treatment), 
               linetype = "solid", show.legend = T, alpha = 0.15) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  annotate("text", x = -1, y = 1, 
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
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
        axis.title.y = element_text(size=25, face="bold", colour = "black"),   
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_text(size=25, face = "bold", color = "black"),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_graph_22

NMDS_Spp_graph_22 = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2), 
                  label = species.scores$species, colour = "black",
                  size = 4, fontface = "bold") +
  annotate("text",                
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
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
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_Spp_graph_22

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ NMDS$Treatment, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, startRow = start_row, colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis[[name]], startRow = start_row)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    # Write the character string
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis[[name]], startRow = start_row, colNames = FALSE)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/Chapter 2 - Fire/pairwise_adonis_22_same_sheet.xlsx", 
             overwrite = TRUE)

################################################################################
################################################################################
################################################################################
##########################     2021 Data       #################################
################################################################################
################################################################################
################################################################################

# Create species pivot table #
Spp = dplyr::select(GRIN_21, ID, Species, Coverage) %>% matrify()

# Transforms Data #
Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_21, ID, Fire) %>% summarise()

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
MDS = metaMDS(Veg_Spp, distance = "bray", trymax = 500, maxit = 999, k=2, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS$stress
stressplot(MDS) 
goodness(MDS)

# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN_22, Species, Group, Native) %>% 
  dplyr::summarize()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Fire, Plot = Treat$ID)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################

envfit_result <- envfit(MDS, Spp)

# Extract vectors and/or factors
vectors_df <- as.data.frame(envfit_result$vectors$arrows)
vectors_pvals <- as.data.frame(envfit_result$vectors$pvals)

if (!is.null(envfit_result$factors)) {
  factors_df <- as.data.frame(envfit_result$factors$centroids)
  factors_pvals <- as.data.frame(envfit_result$factors$pvals)
}
vectors_combined <- cbind(vectors_df, p_value = vectors_pvals)

write.csv(
  vectors_combined, 
  "Figures/Chapter 2 - Fire/NMDS_FitValues_21.csv", 
  row.names = TRUE)

################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Perform Indicator Species Analysis
indicator_results <- indval(Spp, Treat$Fire)

# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)

# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_Tsp = indval_df$Sp,
  Indicator_Value_Tw = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)

# Filter for significant species (p < 0.05)
significant_species1 <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species1$Rank <- 
  ave(significant_species1$Max_Class, significant_species1$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species1, 
  file = "Figures/Chapter 2 - Fire/indicator_species_analysis_results_21.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################


simper_result <- simper(Spp, Treat$Fire, permutations = 99)
simper_summary <- summary(simper_result)
simper_summary

# Extract the first pairwise comparison (adjust as needed)
pairwise_data <- simper_result[[1]]

# Create the data frame with p-values
table_data <- data.frame(
  Species = names(pairwise_data$average),
  Contribution = pairwise_data$average,
  SD = pairwise_data$sd,
  P_value = pairwise_data$p,
  Cumulative_Contribution = pairwise_data$cusum
)

# Sort by contribution
table_data <- 
  table_data[order(table_data$Contribution, decreasing = TRUE), ]

# Add rank column
table_data$Rank <- seq_len(nrow(table_data))

head(table_data)

write.csv(
  table_data, 
  file = "Figures/Chapter 2 - Fire/SIMPER_21.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################

NMDS_graph_21 = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = Treatment),
             alpha = 0.7, size = 5, shape = 21) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, 
                   fill = Treatment, color = Treatment), 
               linetype = "solid", show.legend = T, alpha = 0.15) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  annotate("text", x = -1, y = 1, 
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
           hjust = 0, size = 8) +
  ggtitle("2021") +
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
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_graph_21

NMDS_Spp_graph_21 = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2), 
                  label = species.scores$species, colour = "black",
                  size = 4, fontface = "bold") +
  annotate("text",                
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
           hjust = 0, size = 8) +
  ggtitle("2021") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black", 
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),    
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_Spp_graph_21

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ NMDS$Treatment, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis

# Create a new workbook
wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", 
            x = name, startRow = start_row, colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], startRow = start_row)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
    # Write the character string
    writeData(wb, sheet = "All_Tables", 
              x = pairwise.adonis[[name]], startRow = start_row, colNames = FALSE)
    
    # Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/Chapter 2 - Fire/pairwise_adonis_21_same_sheet.xlsx", 
             overwrite = TRUE)

################################################################################
################################################################################
################################################################################
##########################     2020 Data       #################################
################################################################################
################################################################################
################################################################################

# Create species pivot table #
Spp = dplyr::select(GRIN_20, ID, Species, Coverage) %>% matrify()

# Transforms Data #
Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_20, ID, Fire) %>% summarise()

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
MDS = metaMDS(Veg_Spp, distance = "bray", trymax = 500, maxit = 999, k=2, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS$stress
stressplot(MDS) 
goodness(MDS)

# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN_22, Species, Group, Native) %>% 
  dplyr::summarize()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Fire, Plot = Treat$ID)

# NMDS Graphs
NMDS_20 = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = Treatment),
             alpha = 0.7, size = 5, shape = 21) +
  annotate("text", x = -1, y = 1, 
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
           hjust = 0, size = 8) +
  # geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species)) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, fill = Treatment, 
                   color = Treatment), linetype = "solid", 
               show.legend = T, alpha = 0.15) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  ggtitle("2020") +
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
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_20

NMDS_Spp_graph_20 = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2), 
                  label = species.scores$species, colour = "black",
                  size = 4, fontface = "bold") +
  annotate("text",             
           label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
           hjust = 0, size = 8) +
  ggtitle("2020") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black", 
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),    
        axis.text.x=element_text(size=25, face = "bold", color = "black"),
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.y=element_blank(),
        legend.text=element_text(size=25, face = "bold", color = "black"),
        legend.title=element_text(size=25, face = "bold", color = "black"),
        legend.position="bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_Spp_graph_20

# Perform adonis to test the significance of treatments #
adon.results <- adonis2(Spp ~ NMDS$Treatment, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp ~ Treatment, data = NMDS)
pairwise.adonis

# Create a new workbook

wb <- createWorkbook()

# Add a worksheet
addWorksheet(wb, "All_Tables")

# Initialize starting row
start_row <- 1

# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis)) {
  # Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, startRow = start_row, 
            colNames = FALSE)
  
  # Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
  
  # Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis[[name]], 
              startRow = start_row)
    
# Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis[[name]]) + 2
  } else if (is.character(pairwise.adonis[[name]])) {
# Write the character string
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis[[name]], startRow = start_row, colNames = FALSE)
    
# Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + 2
  }
}

# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/Chapter 2 - Fire/pairwise_adonis_20_same_sheet.xlsx", 
             overwrite = TRUE)

################## Save Figures Above using ggarrange ##########################
NMDS_22_23 = 
  ggarrange(NMDS_graph_22, NMDS_Spp_graph_22,
            NMDS_graph_23, NMDS_Spp_graph_23, ncol = 2, nrow = 2, 
            common.legend = TRUE, legend="bottom")

ggsave("Figures/Chapter 2 - Fire/22-23_NMDS.png", 
       width = 14, height = 12)

NMDS_21_22 = 
  ggarrange(NMDS_graph_21, NMDS_Spp_graph_21, ncol = 2, nrow = 1, 
            common.legend = TRUE, legend="bottom")

ggsave("Figures/Chapter 2 - Fire/21-22_NMDS.png", 
       width = 18, height = 10)
