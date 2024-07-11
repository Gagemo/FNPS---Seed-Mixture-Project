################################################################################
################################################################################
#########################   FNPS - Seed Mixture   ##############################
#########################       NonGrasses        ##############################
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

Data$YID <- paste(Data$Year,Data$ID)
Data$ID_ <- paste(Data$Treatment, Data$ID)

# Orders years and treatments so that they display in same sequence in graphs #
Data$Year = factor(Data$Year, levels=c('1','2'))
Data$Treatment = factor(Data$Treatment, levels=c('C','W', 'BH','BM', 'LH', 'LM'))

#Renames values in Treatment treatments for heat map later #
#Data$Treatment <- recode(Data$Treatment, 
#                    C ="Control", W = "Wiregrass", BH = "Broomsedge High", 
#                    BM = "Broomsedge Medium", LH = "Lovegrass High", 
#                    LM = "Lovegrass Medium")

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
  mutate(Change_abundance = Coverage.y - Coverage.x)

##################################  COVER CAHNGES ##############################
PN = 
  Change_Abundance[which(Change_Abundance$Species == "Paspalum notatum"),]
PN<-as.data.frame(PN)
PN$Treatment<-factor(PN$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = PN)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_PN = PN %>% anova_test(Change_abundance ~ Treatment) %>% 
  add_significance()
anova_PN

lm(formula = Change_abundance ~ Treatment, PN)
tukey_PN <- PN %>% 
  tukey_hsd(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_PN

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=PN)
tmp

PN_change_Box = 
  ggplot(PN, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_PN,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_PN, detailed = TRUE),
       caption = get_pwc_label(tukey_PN)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Change in Coverage", title = "Paspalum notatum")
PN_change_Box
ggsave("Figures/PNchange.png", 
       width = 10, height = 7)

################################################################################
########################### Cynodon dactylon ###################################
################################################################################
CD = 
  Change_Abundance[which(Change_Abundance$Species == "Cynodon dactylon"),]
CD<-as.data.frame(CD)
CD$Treatment<-factor(CD$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = CD)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_CD = CD %>% anova_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_CD)

tukey_CD <- CD %>% 
  tukey_hsd(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_CD

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd), data=CD )
tmp

CD_change_Box = 
  ggplot(CD, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_CD,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_CD, detailed = TRUE),
       caption = get_pwc_label(tukey_CD)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Change in Coverage", title = "Cynodon dactylon")
CD_change_Box
ggsave("Figures/change_CD.png", 
       width = 10, height = 7)

################## Save Figures Above using ggarrange ##########################
Change = 
  ggarrange(PN_change_Box, CD_change_Box, ncol = 1, nrow = 2)
annotate_figure(Change, top = text_grob("", color = "black", 
                                        face = "bold", size = 25))
ggsave("Figures/Change_NonGrass.png", 
       width = 12, height = 8)

