################################################################################
################################################################################
#########################       GRIN - Fire       ##############################
#########################      Weather Data       ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "lubridate", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(agricolae)
library(lubridate)
library(ggpubr)

##########################     Read in  Data      ##############################
data = read.csv("Data/GRIN - Weather Data.csv")
data$Date <- mdy(data$Date)

data$Month <- months(as.Date(data$Date))
data$Year <- as.numeric(format(data$Date,'%Y'))

############################# Year 2023 ########################################
data = filter(data, Year == 2023)

rain = 
  ggplot(data, aes(x = Date)) +
  geom_line(aes(y = Precip.Total*25.4), color = "black", size = 1.25) +
  scale_x_date(date_breaks = "months", date_labels = "%b") +
  geom_vline(xintercept = as.Date(data$Date[36]), linetype="solid", 
               color = "blue", size=1.2) +
  geom_vline(xintercept = as.Date(data$Date[160]), linetype="solid", 
             color = "orange", size=1.2) +
  theme_classic() +
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        axis.title.x = element_text( size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="bold", colour = "black"),   
        axis.text.x=element_text(angle=90, hjust=1,
                                 size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Total Precipitation (mm)", title = "")
rain
ggsave("Figures/Chapter 2 - Fire/precip.png", 
       width = 10, height = 7)

temp = 
  ggplot(data, aes(x = Date)) +
  geom_line(aes(y = (Temp.Max-(32))*(5/9)), color = "red", size = 1.25) +
  geom_line(aes(y = (Temp.Ave-(32))*(5/9)), color = "black", size = 1.25) +
  geom_line(aes(y = (Temp.Min-(32))*(5/9)), color = "lightblue", size = 1.25) +
  scale_x_date(date_breaks = "months" , date_labels = "%b") +
  geom_vline(xintercept = as.Date(data$Date[36]), linetype="solid", 
             color = "blue", size=1.2) +
  geom_vline(xintercept = as.Date(data$Date[160]), linetype="solid", 
             color = "orange", size=1.2) +
  theme_classic() +
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        axis.title.x = element_text( size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="bold", colour = "black"),   
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.line.x = element_blank(),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Temperature C", title = "")
temp

ggsave("Figures/Chapter 2 - Fire/23_temp.png", 
       width = 14, height = 7)

################## Save Figures Above using ggarrange ##########################
ggarrange(temp, rain, ncol = 1, nrow = 2)
ggsave("Figures/Chapter 2 - Fire/23_rainTemp.png", 
       width = 14, height = 12)

##################### 2021 - 2022 Weather Data #################################
data = read.csv("Data/GRIN - Weather Data.csv")
data$Date <- mdy(data$Date)

data$Month <- months(as.Date(data$Date))
data$Year <- as.numeric(format(data$Date,'%Y'))

# Temp data sorting by month #
Max = data %>%
  group_by(Year, Month) %>%
  summarize(Max = max(Temp.Max))

Min = data %>%
  group_by(Year, Month) %>%
  summarize(Min = max(Temp.Min))

Ave = data %>%
  group_by(Year, Month) %>%
  summarize(Ave = mean(Temp.Ave))

prcp = data %>%
  group_by(Year, Month) %>%
  summarize(prcp = sum(Precip.Total))

df = prcp
df$Max = Max$Max
df$Min = Min$Min
df$Ave = Ave$Ave

############################# Weather Graph ####################################
rain = 
  ggplot(df, aes(x = Month)) +
  geom_bar(stat = "identity", aes(y = prcp*25.4), 
           fill = "orange") +
  scale_x_discrete(limits = month.name) +
  facet_wrap(~Year) +
  theme_classic() +
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        axis.title.x = element_text( size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="bold", colour = "black"),   
        axis.text.x=element_text(angle=90, hjust=1,
                                 size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Total Precipitation (mm)", title = "")
rain

ggsave("Figures/Chapter 2 - Fire/precip.png", 
       width = 10, height = 7)

temp = 
  ggplot(df, aes(x = Month)) +
  geom_line(aes(y = (Max-32)*(5/9), group = 1), color = 'red', size = 1.5) +
  geom_line(aes(y = (Min-32)*(5/9), group = 1), color = "blue", size = 1.5) +
  geom_point(aes(y = (Max-32)*(5/9)), size = 2) +
  geom_point(aes(y = (Min-32)*(5/9)), size = 2) +
  scale_x_discrete(limits = month.name)  +
  facet_wrap(~Year) +
  theme_classic() +
  theme(text = element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        axis.title.x = element_text( size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="bold", colour = "black"),   
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = element_text(size = 20, colour = "black", face = "bold"),
        axis.line.x = element_blank(),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Temperature C", title = "")
temp

ggsave("Figures/Chapter 2 - Fire/temp.png", 
       width = 10, height = 7)

################## Save Figures Above using ggarrange ##########################
ggarrange(temp, rain, ncol = 1, nrow = 2)
ggsave("Figures/Chapter 2 - Fire/21-22_rainTemp.png", 
       width = 8, height = 12)

