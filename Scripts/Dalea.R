################################################################################
################################################################################
#########################   FNPS - Seed Mixture   ##############################
#########################         Dalea           ##############################
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
                      "ggsignif", "multcompView", "ggpubr", "rstatix",
                      "vegan", "labdsv")
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

####################### Read in 2021 - 2023 Data  ##############################
Data = read.csv("Data/FNPS - Seed Mixture Project - 2021-2023.csv")
Data$Coverage = as.numeric(Data$Coverage)
Data$Plot = as.character(Data$Plot)

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

DP = filter(Data, Species == "Dalea pinnata")
summary(DP)

# Creates data sets by year #
DP_22 = filter(DP, Year == 1)
DP_23 = filter(DP, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = DP_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
DP_22$Treatment= as.factor(DP_22$Treatment)
DP_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = DP_22 %>% anova_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- DP_22 %>% 
  tukey_hsd(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = DP_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
DP_23$Treatment= as.factor(DP_23$Treatment)
DP_23 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_23 = DP_23 %>% anova_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- DP_23 %>% 
  tukey_hsd(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

################################################################################
################ Create Box Plot for Pityopsis Across Years ####################
################################################################################

## Lovegrass Coverage 2022 Box plot ##
DPBox22 = 
  ggplot(DP_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
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
  labs(x = "", y = "D. pinnata % Coverage", title = "2022")
DPBox22

ggsave("Figures/DP_box22.png", 
       width = 12, height = 8)

## Lovegrass Coverage 2023 Boxplot ##
DPBox23 = 
  ggplot(DP_23, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "D. pinnata % Coverage", title = "2023")
DPBox23

ggsave("Figures/DPbox23.png", 
       width = 12, height = 8)

################## Save Figures Above using ggarrange ##########################
ggarrange(DPBox22, DPBox23, ncol = 2, nrow = 1)
ggsave("Figures/22-23_DPBox.png", 
       width = 12, height = 10)

