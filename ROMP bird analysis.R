
# load required packages ----
library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)
library(glmmTMB)
library(sjPlot)
library(psych)
library(tidyr)
library(car)
library(emmeans)
library(performance)
library(ggpubr)
library(odbc)
library(tidyverse)

# Set working directory - WERG drive
setwd("~/uomShare/wergProj/W12 - Revegetation/ROMP surveys") # change to match path on your computer

# load data
birddata <- read_xlsx("Birdlife_surveys_221024.xlsx")
sitedata <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Site_details")

str(birddata)
glimpse(birddata)
describe(birddata)

birddata <- birddata %>%
  rename(site = 'WptID',
         species = 'Scientific Name',
         type = 'sitetype',
         count = 'Individual Count',
         date = 'Start Date')
birddata <- birddata %>%
  mutate(birddata, type= recode(type, "'Remnant'='Remnant'; 'Works'='Revegetated'"))

birddata_grouped <- birddata %>%
  group_by(type, site, species) %>% summarise(abun = mean(count))
birddatasummary <- birddata_grouped %>%
  group_by(type, site) %>%
  summarise(
    richness = n_distinct(species),
    meanabun = sum(abun))

birdsummarycomplete <-  complete(birddatasummary, fill = list(richness= 0, abun = 0))

# graph native species richness
birdrichnessBOX <- ggplot(
  data = birdsummarycomplete, aes(x=type, y=richness)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site Type', y = "Bird Species Richness") +
  scale_fill_brewer(labels = c("Remnant", "Revegetated"))+
  theme_classic()
birdrichnessBOX
bird_richmodel <- glmmTMB(richness ~ type + (1|site), data = birdsummarycomplete, family = poisson(link = "log"))
summary(bird_richmodel)

birdabunBOX <- ggplot(
  data = birdsummarycomplete, aes(x=type, y=meanabun)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site Type', y = "Bird Species Abundance") +
  scale_fill_brewer(labels = c("Remnant", "Revegetated"))+
  theme_classic()
birdabunBOX
bird_abunmodel <- glmmTMB(meanabun ~ type + (1|site), data = birdsummarycomplete, family = poisson(link = "log"))
summary(bird_abunmodel)
ggsave(bird_abunmodel, filename = "~/uomShare/wergProj/W12 - Revegetation/2024_25 outputs/bird_abunmodel.tiff", width = 16, height = 12, units = "cm", dpi = 300)

########NMDS plot of native plants in ROMP sites#################
ROMPbeltdata24_nat <- filter(ROMPbeltdata24, origin != 'Introduced')
wpt_beltID_df <- ROMPbeltdata24_nat %>% mutate(wpt_beltID=paste(site, belt, sep = "."))
wpt_beltID_cast <- dcast(wpt_beltID_df, sitetype+catchment+wpt_beltID ~ species, mean, value.var="count", fill = 0)
str(wpt_beltID_cast)

wpt_beltID_cast2 <- wpt_beltID_cast[,3:77]
wpt_beltID_cast2$wpt_beltID <- as.numeric(wpt_beltID_cast2$wpt_beltID)
ROMP_NMDS.count <- metaMDS(wpt_beltID_cast2, distance = "bray", k = 2,trymax=100, autotransform = T)
stressplot(ROMP_NMDS.count)
plot(ROMP_NMDS.count)
plot(ROMP_NMDS.count$points)
plot(ROMP_NMDS.count, type = "n")
names_sitetype <- wpt_beltID_cast[,1:1]
orditorp(ROMP_NMDS.count, display = "sites", labels = F, pch = 15, col = c("green", "blue") [as.factor(names_sitetype)], cex = 1)

ROMPNMDS_df <- data.frame(ROMP_NMDS.count$points)
colnames(ROMPNMDS_df) <- c("NMDS1", "NMDS2")
ROMPNMDS_df$sitetype<- factor(wpt_beltID_cast$sitetype) #add group of interest - treatment and ecosystem
ROMPNMDS_df$catchment<- factor(wpt_beltID_cast$catchment)

ggplot(ROMPNMDS_df, aes(x = NMDS1, y = NMDS2, color = sitetype, shape=catchment)) + 
  geom_point(size = 2) +
  xlab("NMDS1") +
  ylab("NMDS2") + 
  ggtitle("Remnant & Revegetated Areas") +
  theme_classic()

data.scores.gg.count = as.data.frame(vegan::scores(ROMPNMDS_df, "sites"))
