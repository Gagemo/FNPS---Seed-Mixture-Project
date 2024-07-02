################################################################################
################################################################################
#########################   FNPS - Seed Mixture   ##############################
#########################     Seeded Forb Cover   ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2022 - 2023        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", "plotrix", 
                      "ggsignif", "multcompView", "ggpubr", "rstatix", "labdsv",
                      "tables")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(agricolae)
library(extrafont)
library(ggsignif)
library(multcompView)
library(ggpubr)
library(plotrix)
library(rstatix)
library(tables)

##########################     Read in  Data       #############################
Data = read.csv("Data/FNPS - Seed Mixture Project - 2021-2023.csv")

Data$Coverage = as.numeric(Data$Coverage)

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
str(Data)
summary(Data)

Data = filter(Data, Species == "Dalea pinnata" | Species == "Liatris gracilis" |
                Species == "Pityopsis trayci")

# Creates data sets by year #
Data_22 = filter(Data, Year == 1)
Data_23 = filter(Data, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Data_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Data_22$Treatment= as.factor(Data_22$Treatment)
Data_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Data_22 %>% anova_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Data_22 %>% 
  tukey_hsd(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Data_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Data_23$Treatment= as.factor(Data_23$Treatment)
Data_23 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_23 = Data_23 %>% anova_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Data_23 %>% 
  tukey_hsd(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

################################################################################
################ Create Box Plot Across Years ##################################
################################################################################

## Lovegrass Coverage 2022 Box plot ##
Box22 = 
  ggplot(Data_22, aes(x = Treatment, y = Coverage), colour = Species) +
  geom_boxplot(aes(fill=Species), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% Coverage", title = "2022")
Box22

ggsave("Figures/forb_box22.png", 
       width = 12, height = 8)

## Lovegrass Coverage 2023 Boxplot ##
Box23 = 
  ggplot(Data_23, aes(x = Treatment, y = Coverage), colour = Species) +
  geom_boxplot(aes(fill=Species), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "% Coverage", title = "2023")
Box23

ggsave("Figures/Ptbox23.png", 
       width = 12, height = 8)

################## Save Figures Above using ggarrange ##########################
ggarrange(Box22, Box23, ncol = 2, nrow = 1)
ggsave("Figures/22-23_ForbBox.png", 
       width = 12, height = 10)

