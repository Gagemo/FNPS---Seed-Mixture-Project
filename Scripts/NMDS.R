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
                        "ggrepel","ggsignif", "multcompView", "ggpubr", 
                        "rstatix", 'rmarkdown', "labdsv", "pairwiseAdonis", 
                        "devtools", "knitr", "tables", "openxlsx", "labdsv")
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
  
  ##########################     Read in 2022-2023 Data       ####################
  Data = read.csv("Data/FNPS - Seed Mixture Project - 2021-2023.csv")
  Data$Coverage = as.numeric(Data$Coverage)
  
  str(Data)
  summary(Data)
  
  # Reclasifys coverage data (CV) from 1-10 scale to percent scale #
  Data <- mutate(Data, Coverage = case_when(
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
  
  Data$Treatment = factor(Data$Treatment, 
                          levels=c('BH', 'BM', 'LH', 'LM', 'W', "C"))
  
  Data_22 = filter(Data, Year==1)
  Data_23 = filter(Data, Year==2)
  
################################################################################
################################################################################
################################################################################
########################### 2022 Data ##########################################
################################################################################
################################################################################
################################################################################
  
# Create species pivot table #
Spp = dplyr::select(Data_22, ID, Species, Coverage) %>% matrify()
  
# Transforms Data #
Veg_Spp = vegdist(Spp, method = 'bray')
  
# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(Data, ID, Treatment) %>% dplyr::summarize()
  
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
species_groups = group_by(Data_22, Species, Type) %>% 
  dplyr::summarize()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment,Plot = Treat$ID)

################################################################################
############## NMDS Vector fitting with envfit##################################
################################################################################
  
envfit_result <- envfit(MDS, Spp)
envfit_result

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
    "Figures/NMDS_FitValues_22.csv", 
    row.names = TRUE)
  
################################################################################
######################### Indicator Species Analysis ###########################
################################################################################
  
# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)
  
# Remove absent species
Spp_filtered <- Spp[, -absent_species]
  
# Perform Indicator Species Analysis
indicator_results <- indval(Spp_filtered, Treat$Treatment)
  
# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)
  
# Combine data into a single data frame
combined_df <- data.frame(
    Species = species_names,
    Indicator_Value_C = indval_df$C,
    Indicator_Value_BH = indval_df$BH,
    Indicator_Value_BM = indval_df$BM,
    Indicator_Value_LH = indval_df$LH,
    Indicator_Value_LM = indval_df$LM,
    Indicator_Value_W = indval_df$W,
    Max_Class = indicator_results$maxcls,
    Indicator_Class = indicator_results$indcls,
    P_Value = indicator_results$pval
  )
  
# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]
  
# Add rank column based on indicator value in the max class
significant_species$Rank <- 
    ave(significant_species$Max_Class, 
        significant_species$Max_Class, FUN = rank)
  
# Save the combined data frame to CSV
write.csv(
    significant_species, 
    file = "Figures/indicator_species_analysis_results_22.csv", 
    row.names = TRUE)
  
################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################
  
simper_result <- simper(Spp, Treat$Treatment)
simper_summary <- summary(simper_result)
  
# Extract the first pairwise comparison (adjust as needed)
simper_result <- simper(Spp, Treat$Treatment, permutations = 99)
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
  file = "Figures/SIMPER_22.csv", 
  row.names = TRUE)
  
################################################################################
#############################NMDS Graphs########################################
################################################################################
  
# NMDS Graphs
NMDS_graph_22 = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
  scale_color_manual(labels=c('BH', 'BM', 'LH', 'LM', 'W', "C"),
                     values=c("#663333", "#FF9966", "#006600",
                              "#99FF99", "#CC0000", "#330099")) +
  scale_fill_manual(labels=c('BH', 'BM', 'LH', 'LM', 'W', "C"),
                    values=c("#663333", "#FF9966", "#006600",
                             "#99FF99", "#CC0000", "#330099")) +
  scale_shape_manual(labels=c('BH', 'BM', 'LH', 'LM', 'W', "C"),
                     values=c(22, 23, 24, 25, 21, 10)) +
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
  guides(
    #   fill = guide_legend(nrow = 1), # Removed fill
    shape = guide_legend(nrow = 1)
  ) +
  labs(x = "MDS1", y = "MDS2",
       color = "Treatment", # Added color to the labs
       fill = "Treatment") # added fill to the labs
NMDS_graph_22
  
NMDS_Spp_graph_22 = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2), 
                  label = species.scores$species, colour = "black",
                  size = 5, fontface = "bold") +
  annotate("text", label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
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
        axis.text.x = element_text(size=25, face = "bold", color = "black"),
        axis.text.y = element_text(size=25, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size=25, face = "bold", color = "black"),
        legend.title = element_text(size=25, face = "bold", color = "black"),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1)) +  # Adjust the number of rows to control orientation
  labs(x = "MDS1", y = "MDS2", color = "", 
       fill = "")
NMDS_Spp_graph_22
  
# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray")
print(adon.results)
write.csv.tabular(adon.results, "Figures/adonis_22.csv")
pairwise.adonis_22<-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis_22
  
#save tables
# Create a new workbook
  
wb <- createWorkbook()
  
 # Add a worksheet
  addWorksheet(wb, "All_Tables")
  
# Initialize starting row
  start_row <- 1
  
# Loop through the list of tables and add each to the same sheet
  for (name in names(pairwise.adonis_22)) {
# Add table name as a header
    writeData(wb, sheet = "All_Tables", x = name, 
              startRow = start_row, colNames = FALSE)
    
# Increment the starting row to leave a gap between the header and the table
    start_row <- start_row + 1
    
    # Check if the element is a data frame or a character string
    if (is.data.frame(pairwise.adonis_22[[name]])) {
      # Write the table
      writeData(wb, sheet = "All_Tables", 
                x = pairwise.adonis_22[[name]], startRow = start_row)
      
# Increment the starting row for the next table, adding a few extra rows for spacing
      start_row <- start_row + nrow(pairwise.adonis_22[[name]]) + 2
    } else if (is.character(pairwise.adonis_22[[name]])) {
      # Write the character string
      writeData(wb, sheet = "All_Tables", 
                x = pairwise.adonis_22[[name]], 
                startRow = start_row, colNames = FALSE)
      
# Increment the starting row for the next table, adding a few extra rows for spacing
      start_row <- start_row + 2
    }
  }
  
  # Save the workbook to an Excel file
  saveWorkbook(wb, "Figures/pairwise_adonis_22_same_sheet.xlsx", overwrite = TRUE)
  
################################################################################
################################################################################
################################################################################
########################### 2023 Data ##########################################
################################################################################
################################################################################
################################################################################

# Create species pivot table #
Spp = dplyr::select(Data_23, ID, Species, Coverage) %>% matrify()

# Transforms Data #
Veg_Spp = vegdist(Spp, method = 'bray')

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(Data_23, ID, Treatment) %>% dplyr::summarize()
Treat
  
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
species_groups = group_by(Data_23, Species, Type) %>% 
  dplyr::summarize()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- as.data.frame(wascores(MDS$points, Spp))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)
  
# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment,Plot = Treat$ID)

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
  "Figures/NMDS_FitValues_23.csv", 
  row.names = TRUE)
  
################################################################################
######################### Indicator Species Analysis ###########################
################################################################################

# Identify columns with only zeros
absent_species <- which(colSums(Spp) == 0)

# Remove absent species
Spp_filtered <- Spp[, -absent_species]
  
# Perform Indicator Species Analysis
indicator_results <- indval(Spp_filtered, Treat$Treatment)
  
# Extract data frames
indval_df <- indicator_results$indval
species_names <- rownames(indval_df)
  
# Combine data into a single data frame
combined_df <- data.frame(
  Species = species_names,
  Indicator_Value_C = indval_df$C,
  Indicator_Value_BH = indval_df$BH,
  Indicator_Value_BM = indval_df$BM,
  Indicator_Value_LH = indval_df$LH,
  Indicator_Value_LM = indval_df$LM,
  Indicator_Value_W = indval_df$W,
  Max_Class = indicator_results$maxcls,
  Indicator_Class = indicator_results$indcls,
  P_Value = indicator_results$pval
)
  
# Filter for significant species (p < 0.05)
significant_species <- combined_df[combined_df$P_Value < 0.05, ]

# Add rank column based on indicator value in the max class
significant_species$Rank <- 
  ave(significant_species$Max_Class, 
      significant_species$Max_Class, FUN = rank)

# Save the combined data frame to CSV
write.csv(
  significant_species, 
  file = "Figures/indicator_species_analysis_results_23.csv", 
  row.names = TRUE)

################################################################################
####################Similarity Percentages (SIMPER)#############################
################################################################################
  

simper_result <- simper(Spp, Treat$Treatment)
simper_summary <- summary(simper_result)

# Extract the first pairwise comparison (adjust as needed)
simper_result <- simper(Spp, Treat$Treatment, permutations = 99)
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
  file = "Figures/SIMPER_23.csv", 
  row.names = TRUE)

################################################################################
#############################NMDS Graphs########################################
################################################################################
  
# NMDS Graphs
NMDS_graph_23 <-
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2,
                              fill = Treatment, shape = Treatment),
             alpha = 0.7, size = 5) +
  scale_color_manual(labels=c('BH', 'BM', 'LH', 'LM', 'W', "C"),
                     values=c("#663333", "#FF9966", "#006600",
                              "#99FF99", "#CC0000", "#330099")) +
  scale_fill_manual(labels=c('BH', 'BM', 'LH', 'LM', 'W', "C"),
                    values=c("#663333", "#FF9966", "#006600",
                             "#99FF99", "#CC0000", "#330099")) +
  scale_shape_manual(labels=c('BH', 'BM', 'LH', 'LM', 'W', "C"),
                     values=c(22, 23, 24, 25, 21, 10)) +
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
  guides(
    #   fill = guide_legend(nrow = 1), # Removed fill
    shape = guide_legend(nrow = 1)
  ) +
  labs(x = "MDS1", y = "MDS2",
       color = "Treatment", # Added color to the labs
       fill = "Treatment") # added fill to the labs
NMDS_graph_23

  
NMDS_Spp_graph_23 = 
  ggplot() +
  geom_text_repel(data = species.scores, aes(x = MDS1, y = MDS2), 
                  label = species.scores$species, colour = "black",
                  size = 5, fontface = "bold") +
  annotate("text", label = paste0("Stress: ", format(MDS$stress, digits = 2)), 
           hjust = 0, size = 8) +
  ylim(-0.25, 0.50) +
  ggtitle("2023") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, color="black", 
                                  size=25, face="bold"),
        axis.title.x = element_text(size=25, face="bold", colour = "black"),    
        axis.text.x = element_text(size=25, face = "bold", color = "black"),
        axis.text.y = element_text(size=25, face = "bold", color = "black"),
        axis.title.y = element_blank(),
        legend.text = element_text(size=25, face = "bold", color = "black"),
        legend.title = element_text(size=25, face = "bold", color = "black"),
        legend.position = "bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "MDS1", y = "MDS2", color = "", 
       fill = "")
NMDS_Spp_graph_23
  
# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Veg_Spp ~ Treatment, data = NMDS, method="bray",
                        permutations = 999)
print(adon.results)
write.csv.tabular(adon.results, "Figures/adonis_23.csv")
pairwise.adonis_23<-pairwise.adonis2(Veg_Spp ~ Treatment, data = NMDS)
pairwise.adonis_23
  
#save tables
# Create a new workbook
  
wb <- createWorkbook()
  
# Add a worksheet
addWorksheet(wb, "All_Tables")
  
# Initialize starting row
start_row <- 1
  
# Loop through the list of tables and add each to the same sheet
for (name in names(pairwise.adonis_23)) {
# Add table name as a header
  writeData(wb, sheet = "All_Tables", x = name, 
            startRow = start_row, colNames = FALSE)
  
# Increment the starting row to leave a gap between the header and the table
  start_row <- start_row + 1
    
# Check if the element is a data frame or a character string
  if (is.data.frame(pairwise.adonis_23[[name]])) {
    # Write the table
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis_23[[name]], 
              startRow = start_row)
      
# Increment the starting row for the next table, adding a few extra rows for spacing
    start_row <- start_row + nrow(pairwise.adonis_23[[name]]) + 2
    } else if (is.character(pairwise.adonis_23[[name]])) {
      # Write the character string
    writeData(wb, sheet = "All_Tables", x = pairwise.adonis_23[[name]], startRow = start_row, colNames = FALSE)
      
# Increment the starting row for the next table, adding a few extra rows for spacing
      start_row <- start_row + 2
    }
  }
  
# Save the workbook to an Excel file
saveWorkbook(wb, "Figures/pairwise_adonis_23_same_sheet.xlsx", overwrite = TRUE)
  
# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Spp ~ NMDS$Treatment, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp ~ Treatment, data = NMDS)  
pairwise.adonis %>% as.data.frame() %>% write.csv(file = "Figures/adonis_2023.csv")
  
################## Save Figures Above using ggarrange ##########################
NMDS_23_Combine = ggarrange(NMDS_graph_23, NMDS_Spp_graph_23,
                            ncol = 2, nrow = 1, common.legend = TRUE,
                            legend = "bottom",
                            legend.grob = get_legend(NMDS_graph_23))  # Extract legend from NMDS_graph_23
NMDS_23_Combine


ggsave("Figures/NMDS_23_Combine.png", 
       width = 14, height = 8)

################## Save Figures Above using ggarrange ##########################
NMDS_Combine = ggarrange(NMDS_graph_22, NMDS_Spp_graph_22,
                         NMDS_graph_23, NMDS_Spp_graph_23,
                         ncol = 2, nrow = 2, common.legend = TRUE,
                         legend = "bottom",
                         legend.grob = get_legend(NMDS_graph_23))  # Extract legend from NMDS_graph_23
NMDS_Combine


ggsave("Figures/NMDS_Combine.png", 
       width = 14, height = 16)
  
