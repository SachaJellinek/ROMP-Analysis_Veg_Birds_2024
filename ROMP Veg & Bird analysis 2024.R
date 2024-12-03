### THIS FILE SUMMARISES DATA COLLECTED FOR ROMP 2024 DATA ON REVEGETATION OUTCOMES 

# Project: ROMP 2024 STUDY ON REVEGETATION OUTCOMES  #
# Required data files:  XXXXXX
# 
# Purpose: statistical analysis of revegetation outcomes                  #
# Author: SACHA JELLINEK                                            #
# Date: 20240415                                             #


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
sitedata <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Site_details")
#beltdata_update <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Belt_trans_species_update")
beltdata24 <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Belt_trans_species")
#heightdata <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Tree Heights")
#dbhdata <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "DBH")
stratadata<- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Structure_poll")
canopydata <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Canopy")
quadratdata24 <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "Cover_estimates")
visitdata <- read_xlsx("ROMP_data_v4_05092024.xlsx", sheet = "iVisit")
floravic <- read_excel("ROMP_data_v4_05092024.xlsx", sheet = "Flora_vic", col_types = "text")
as.numeric(floravic$TAXON_ID)

str(quadratdata24)
glimpse(sitedata)
describe(sitedata)

# bring in site type into beltdata dataframe
beltdata24$sitetype <- sitedata$'Site_type'[match(beltdata24$'WptID', sitedata$'WptID')]
# bring in site pair into beltdata dataframe
beltdata24$pair <- sitedata$'Target_WptID'[match(beltdata24$'WptID', sitedata$'WptID')]
beltdata24$catchment <- sitedata$'Catchment'[match(beltdata24$'WptID', sitedata$'WptID')]
#beltdata$visit <- visitdata$'Visit_number'[match(beltdata$'VisitID', visitdata$'VisitID')]
beltdata24$origin <- floravic$'ORIGIN'[match(beltdata24$'Sp_name', floravic$'SCI_NAME')]
beltdata24$Month_Yr <- format(as.Date(beltdata24$Date), "%Y-%m")
beltdata24 <- mutate(beltdata24, 
                     origin = ifelse(origin == "Native but some stands may be alien", "Native", origin)) 
ROMPbeltdata24 <- beltdata24 %>%
  rename(site = 'WptID',
         species = 'Sp_name',
         type = 'Type (reveg, existing, recruit)',
         alive = 'Count_alive',
         dead = 'Count_dead',
         visit = 'Visit_number', 
         belt = 'Belt_ID')

# identify factors
ROMPbeltdata24$site <- as.factor(ROMPbeltdata24$site)
ROMPbeltdata24$species <- as.factor(ROMPbeltdata24$species)
ROMPbeltdata24$origin <- as.factor(ROMPbeltdata24$origin)
ROMPbeltdata24$sitetype <- as.factor(ROMPbeltdata24$sitetype)
ROMPbeltdata24$pair <- as.factor(ROMPbeltdata24$pair)
describeBy(ROMPbeltdata24, ROMPbeltdata24$sitetype)

# Check species without any records
missing_val <- ROMPbeltdata24 %>% filter(alive==0, dead==0)
#writexl::write_xlsx(missing_val, '~/uomShare/wergProj/W12 - Revegetation/ROMP surveys/2024_25 outputs/missing_val.xlsx')
# Check which columns have missing values
sapply(ROMPbeltdata24, function(x) sum(is.na(x)))
missing_types <- filter(ROMPbeltdata24, !complete.cases(ROMPbeltdata24$type))
#writexl::write_xlsx(missing_types, '~/uomShare/wergProj/W12 - Revegetation/ROMP surveys/2024_25 outputs/missing_types.xlsx')
listed_type_dead <- ROMPbeltdata24 %>% filter(type=="dead")
#writexl::write_xlsx(listed_type_dead, '~/uomShare/wergProj/W12 - Revegetation/ROMP surveys/2024_25 outputs/listed_type_dead.xlsx')

# Count of sites with number of belt transects
count_belt <- ROMPbeltdata24 %>% group_by(Month_Yr, site, Site_ID, catchment, visit, sitetype, belt) %>% count()
belt_by_visit <- count_belt %>% group_by(Month_Yr, site, Site_ID, catchment, visit) %>% count(visit)
belt_by_visit
site_visit <- belt_by_visit %>% group_by(visit) %>% count(visit)
type_visit <- count_belt %>% group_by(Month_Yr, site, Site_ID, catchment, visit, sitetype) %>% count(visit)
type_visit <- type_visit %>% group_by(Month_Yr, visit, catchment, sitetype) %>% count(visit)
writexl::write_xlsx(type_visit, '~/uomShare/wergProj/W12 - Revegetation/ROMP surveys/2024_25 outputs/type_visit.xlsx')
count_belt_2 <- ROMPbeltdata24 %>% group_by(visit, site) %>% count(site)
count_belt_3 <- count_belt_2 %>% group_by(site)%>% count(site)

#### Q1/2a. NATIVE TREE AND SHRUB Richness ----

ROMPbeltdata24 <- ROMPbeltdata24 %>% filter(!is.na(site))
ROMPbeltdata24 <- ROMPbeltdata24 %>% replace_na(list(dead = 0))
ROMPbeltdata24 <- ROMPbeltdata24 %>% replace_na(list(alive = 0))
ROMPbeltdata24 <- ROMPbeltdata24 %>% mutate(count= alive + dead)

ROMPbeltdata24 <- filter(ROMPbeltdata24, species != 'Unknown')
ROMPbeltdata24 <- filter(ROMPbeltdata24, species != 'Nil')

#beltdata %>% filter_at(vars(Sp_name), any_vars(. %in% c("(?i)^Allocasuarina"))) 
#belt1 <- beltdata %>% filter_at(vars(Sp_name),select(matches("^Allocasuarina")))

count_spabun <- ROMPbeltdata24 %>% group_by(site, sitetype, species, origin) %>% summarise(sumofalive = sum(alive))
count_spabun2 <- ROMPbeltdata24 %>% group_by(sitetype, species, origin) %>% summarise(sumofalive = sum(alive))
#Calculate overall species richness
count_sprich <- count_spabun2 %>% group_by(sitetype, origin) %>% summarize(richness=n())
ggplot(count_sprich, aes(x = sitetype, y = richness, fill = origin)) +
  geom_bar(stat="identity", position=position_dodge())

#Calculate mean Species Richness with standard error bars
mean_sprich <- count_spabun %>% 
  group_by(site, sitetype, origin) %>% summarize(richness=n())
mean_sprichSE <- mean_sprich%>% 
  group_by(sitetype, origin)%>%
  summarise(mean = mean(richness),num_obs=n(),sum_habitat_richness=sum(richness),
            sd_habitat_richness= sd(richness),se_mean=sd_habitat_richness/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)
mean_sprichCI <- mean_sprich%>% 
  group_by(sitetype, origin)%>%
  summarise(mean.rich = mean(richness, na.rm = TRUE),
            sd.rich = sd(richness, na.rm = TRUE),
            n.rich = n()) %>%
  mutate(se.rich = sd.rich / sqrt(n.rich),
         lower.ci.rich = mean.rich - qnorm(0.975) * se.rich,
         upper.ci.rich = mean.rich + qnorm(0.975) * se.rich)
#Plot sp rich with error bars
ggplot(data=mean_sprichSE, aes(x = reorder(sitetype, -mean),y=mean,fill=origin))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=sitetype, ymin= se_lower , ymax = se_upper), width=0.2, position=position_dodge(.9)) + 
  ylab("Species Richness") + 
  xlab('Site type')+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))
#Plot sp rich with 95% confidence intervals
ggplot(data=mean_sprichCI, aes(x = reorder(sitetype, -mean.rich),y=mean.rich,fill=origin))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=sitetype, ymin= lower.ci.rich , ymax = upper.ci.rich), width=0.2, position=position_dodge(.9)) + 
  ylab("Species Richness") + 
  xlab('Site type')+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))

# Mean species richness by origin and visit
count_spabunVisit <- ROMPbeltdata24 %>% group_by(site, sitetype, visit, species, origin) %>% summarise(sumofalive = sum(alive))
mean_sprichVisit <- count_spabunVisit %>% 
  group_by(site, sitetype, visit, origin) %>% summarize(richness=n())
mean_sprichVisit <- tidyr::unite(mean_sprichVisit,"vis_orig",visit, origin,remove = F)
mean_sprichVisitCI <- mean_sprichVisit%>% 
  group_by(sitetype, vis_orig)%>%
  summarise(mean.rich = mean(richness, na.rm = TRUE),
            sd.rich = sd(richness, na.rm = TRUE),
            n.rich = n()) %>%
  mutate(se.rich = sd.rich / sqrt(n.rich),
         lower.ci.rich = mean.rich - qnorm(0.975) * se.rich,
         upper.ci.rich = mean.rich + qnorm(0.975) * se.rich)
mean_sprich2SE <- mean_sprichVisit%>% 
  group_by(sitetype, vis_orig)%>%
  summarise(mean = mean(richness),num_obs=n(),sum_habitat_richness=sum(richness),
            sd_habitat_richness= sd(richness),se_mean=sd_habitat_richness/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)
#Plot sp rich with error bars
ggplot(data=mean_sprich2SE, aes(x = reorder(sitetype, -mean),y=mean,fill=vis_orig))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=sitetype, ymin= se_lower , ymax = se_upper), width=0.2, position=position_dodge(.9)) + 
  ylab("Species Richness") + 
  xlab('Site type')+
  labs(fill='Visit and Sp origin')  +
  scale_fill_manual(values=c("pink", "lightblue", "pink", "lightblue","pink", "lightblue"))+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))
ggplot(data=mean_sprichVisitCI, aes(x = reorder(sitetype, -mean.rich),y=mean.rich,fill=vis_orig))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=sitetype, ymin= lower.ci.rich , ymax = upper.ci.rich), width=0.2, position=position_dodge(.9)) + 
  ylab("Species Richness") + 
  xlab('Site type')+
  labs(fill='Visit and Sp origin')  +
  scale_fill_manual(values=c("pink", "lightblue", "pink", "lightblue","pink", "lightblue"))+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))

# Mean species abundance accross remnant and reveg areas
mean_spabun <- count_spabun %>% 
  group_by(site, sitetype, origin) %>% summarize(abun=sum(sumofalive))
mean_spabunSE <- mean_spabun%>% 
  group_by(sitetype, origin)%>%
  summarise(mean = mean(abun),num_obs=n(),sum_habitat_abun=sum(abun),
            sd_habitat_abun= sd(abun),se_mean=sd_habitat_abun/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)
ggplot(data=mean_spabunSE, aes(x = reorder(sitetype, -mean),y=mean,fill=origin))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=sitetype, ymin= se_lower , ymax = se_upper), width=0.2, position=position_dodge(.9)) + 
  ylab("Species Abundance") + 
  xlab('Site type')+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))

# The survival between planting and first surveys
beltdata24$catchment <- sitedata$'Catchment'[match(beltdata24$'site', sitedata$'WptID')]

surv <- filter(beltdata24, type == "reveg")
surv <- filter(surv, sitetype != 'Remnant')
surv <- filter(surv, visit != 1)
surv <- surv %>% 
  dplyr::mutate(alive = as.integer(alive), 
                count = as.integer(count))
surv2 <- surv %>% 
  group_by(site, visit) %>% summarize(alive=sum(alive), total=sum(count))
Surv_yeardiff <- surv2 %>% group_by(site) %>% mutate(filterout = sum(visit))
Surv_yeardiff1 <- Surv_yeardiff %>% filter(filterout==5)

#Surv_yeardiff1$visit <- as.factor(Surv_yeardiff1$visit)

survsum <- Surv_yeardiff1 %>% mutate(percent=alive/total)
mean_survSE <- survsum%>% 
  group_by(visit)%>%
  summarise(mean = mean(percent),num_obs=n(),sum_surv=sum(percent),
            sd_surv= sd(percent),se_mean=sd_surv/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)

mean_survivalCI <- survsum%>% 
  group_by(visit)%>%
  summarise(mean.perc = mean(percent, na.rm = TRUE),
            sd.perc = sd(percent, na.rm = TRUE),
            n.perc = n()) %>%
  mutate(se.perc = sd.perc / sqrt(n.perc),
         lower.ci.perc = mean.perc - qnorm(0.975) * se.perc,
         upper.ci.perc = mean.perc + qnorm(0.975) * se.perc)
ggplot(data=mean_survSE, aes(x = visit, y=mean))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=visit, ymin= se_lower , ymax = se_upper), width=0.2, position=position_dodge(.9)) + 
  ylab("% Survival") + 
  xlab('Visit')+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))

# Survival accross catchments
survcatch <- surv %>% 
  group_by(site, catchment, visit) %>% summarize(alive=sum(alive), total=sum(count))
survsum <- survcatch %>% mutate(percent=alive/total)
mean_survCatchCI <- survsum%>% 
  group_by(catchment, visit)%>%
  summarise(mean.perc = mean(percent, na.rm = TRUE),
            sd.perc = sd(percent, na.rm = TRUE),
            n.perc = n()) %>%
  mutate(se.perc = sd.perc / sqrt(n.perc),
         lower.ci.perc = mean.perc - qnorm(0.975) * se.perc,
         upper.ci.perc = mean.perc + qnorm(0.975) * se.perc)
ggplot(data=mean_survCatchCI, aes(x = visit, y=mean.perc, fill = catchment))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar (aes(x=visit, ymin= lower.ci.perc , ymax = upper.ci.perc), width=0.2, position=position_dodge(.9)) + 
  ylab("% Survival") + 
  xlab('Visit')+
  theme_bw()+
  theme(text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 45,vjust=0.5))

Surv_yeardiff <- surv2 %>% group_by(site) %>% mutate(filterout = sum(visit))
Surv_yeardiff1 <- Surv_yeardiff %>% filter(filterout==5)
Surv_yeardiff <- Surv_yeardiff1 %>%
  group_by(site) %>%
  mutate(overall_surv = first(alive) - alive)
mean_survSE <- filter(Surv_yeardiff, visit != 2)
mean_survSE <- Surv_yeardiff%>%
  summarise(mean = mean(overall_surv),num_obs=n(),sum_surv=sum(overall_surv),
            sd_habitat_richness= sd(overall_surv),se_mean=sd_habitat_richness/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)

##############################################################
# Joe code - create summary of species per sitetype by no of sites
# Using visit 2 data - Post reveg

speciessummarybysite <- ROMPbeltdata24 %>%
  group_by(species, visit, origin, sitetype, site) %>%
  summarise(noofstems = sum(count))

# Change the visit in remnnat sites to 2 and filter data by second visit only
speciessummarybysite$visit[speciessummarybysite$sitetype == 'Remnant'] <- 2
speciessummarybysite <- speciessummarybysite %>% filter(visit == 2)

speciessummary <- speciessummarybysite %>%
  group_by(species, origin) %>%
  summarise(
    noofremsites = sum(sitetype == "Remnant"),
    noofworkssites = sum(sitetype == "Revegetation"),
    meannoofstemsrem = mean(noofstems[sitetype == "Remnant"], na.rm = TRUE),
    meannoofstemsworks = mean(noofstems[sitetype == "Revegetation"], na.rm = TRUE))

speciessummary <- speciessummary %>%
  mutate(sitedifremminusworks = noofremsites - noofworkssites)
# write.csv(speciessummary, "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/specieslistforTableS1.csv", row.names=FALSE)
# bring in site type into beltdata dataframe
belt_cover_df <- left_join(beltdata24, cover_pellet, by = c("site", "visit", "belt", "sitetype"))
belt_cover_df$visit[belt_cover_df$sitetype == 'Remnant'] <- 2
belt_cover <- belt_cover_df %>% filter(visit == 2)
# summarise
beltdatasummary <- belt_cover %>%
  group_by(sitetype, pair, site, origin) %>%
  summarise(
    richness = n_distinct(species),
    nostems = sum(count),
    nativevegcoveragepercentage = mean(native))

beltsummarycomplete <-  complete(beltdatasummary, origin, fill = list(richness= 0, nostems = 0, norecruits = 0))

# filter to consider only native
beltdatasummarynative <- filter(beltdatasummary, origin == "Native")
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)

# graph native species richness
nativetreerichness <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site, -richness), y=richness, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Richness of native trees and shrubs") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreerichness
nativetreerichness <- nativetreerichness + theme(legend.position = "none")

#ggsave(nativetreerichness, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativetreerichness.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model it
nativerichnessbysitetypemodel <- glmmTMB(richness ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativerichnessbysitetypemodel)
r2(nativerichnessbysitetypemodel)


# BOX PLOT OF NATIVE RICHNESS

nativetreerichnessBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=richness, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Native woody richness") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativetreerichnessBOX
nativetreerichnessBOX <- nativetreerichnessBOX + theme(legend.position = "none")


#### Q1/2a. NATIVE TREE AND SHRUB ABUNDANCE ----

nativetreeabundance <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site,-nostems), y=nostems, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Abundance of native trees and shrubs") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreeabundance

# but need to consider TRANSECT length!!!
belt_lengthdata <- read_xlsx("ROMP_data_v3_20240718_SP.xlsx", sheet = "Belt_length")
belt_lengthdata$site <- as.factor(belt_lengthdata$site)
beltdatasummarynative <- left_join(beltdatasummarynative, belt_lengthdata, by = c("site"))
beltdatasummarynative %>% left_join(belt_lengthdata)

# incorporate belt area into count data
beltdatasummarynative$stemsperha <- round(beltdatasummarynative$nostems/beltdatasummarynative$belt_length*10000)
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)

# plot with new data
nativetreeabundance <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site,-stemsperha), y=stemsperha, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Native tree and shrub density (stems/ha)") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreeabundance

# model it - old analyses
nativeabundancebysitetypemodel <- glm(nostems ~ sitetype, data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativeabundancebysitetypemodel)
r2(nativeabundancebysitetypemodel)


# model it - new analyses (considers total surveyed area)
nativeabundanceperhabysitetypemodel <- glmmTMB(stemsperha ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativeabundanceperhabysitetypemodel)
r2(nativeabundanceperhabysitetypemodel)
1/exp(-0.566653)

# BOX PLOT OF NATIVE WOODY PLANT ABUNDANCE
str(beltdatasummarynative)
nativetreeabundanceBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=stemsperha, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Native woody plant density") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativetreeabundanceBOX
nativetreeabundanceBOX <- nativetreeabundanceBOX + theme(legend.position = "none")

### Survivorship analysis in revegetated areas
beltsurvdat <- filter(belt_cover_df, alive != 0)
beltdatasummary_surv <- beltsurvdat %>%
  group_by(sitetype, pair, site, visit, type, origin) %>%
  summarise(
    richness = n_distinct(species),
    nostems = sum(alive),
    nativevegcoveragepercentage = mean(native))

beltsummarycomplete_surv <-  complete(beltdatasummary_surv, origin, fill = list(richness= 0, nostems = 0, norecruits = 0))

# filter to consider only native
beltdatasummaryreveg <- filter(beltdatasummary_surv, type == "reveg")
describeBy(beltdatasummaryreveg, beltdatasummaryreveg$visit)

beltdatasummaryreveg$catchment <- sitedata$'Catchment'[match(beltdatasummaryreveg$'site', sitedata$'WptID')]
beltdatasummaryreveg <- filter(beltdatasummaryreveg, sitetype != 'Remnant')
beltdatasummaryreveg <- filter(beltdatasummaryreveg, visit != 1)
#beltdatasummaryreveg <- beltdatasummaryreveg %>% 
#  dplyr::mutate(alive = as.integer(alive), 
#                count = as.integer(count))
Yr1_survival<- beltdatasummaryreveg %>% 
  group_by(site, visit) %>% summarize(richness=sum(richness), stems=sum(nostems))
Yr1Surv_year_2_3 <- Yr1_survival %>% group_by(site) %>% mutate(filterout = sum(visit))
Yr1Surv_year_2_3 <- Yr1Surv_year_2_3 %>% filter(filterout==5)
#Yr1Surv_year_2_3$visit <- as.factor(Yr1Surv_year_2_3$visit)

abun_stemsperha_reveg <- left_join(Yr1Surv_year_2_3, belt_lengthdata, by = c("site"))
abun_stemsperha_reveg %>% left_join(belt_lengthdata)
# incorporate belt area into count data
abun_stemsperha_reveg$stemsperha <- round(abun_stemsperha_reveg$stems/abun_stemsperha_reveg$belt_length*10000)
describeBy(abun_stemsperha_reveg, abun_stemsperha_reveg$visit)
#survsum <- surv2 %>% mutate(percent=alive/total)
mean_surv_rich <- Yr1Surv_year_2_3%>% 
  group_by(visit)%>%
  summarise(mean = mean(richness),num_obs=n(),sum_surv=sum(richness),
            sd_surv= sd(richness),se_mean=sd_surv/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)

mean_surv_abun <- abun_stemsperha_reveg%>% 
  group_by(visit)%>%
  summarise(mean = mean(stemsperha),num_obs=n(),sum_surv=sum(stemsperha),
            sd_surv= sd(stemsperha),se_mean=sd_surv/sqrt(num_obs),
            se_upper=mean+se_mean,se_lower=mean-se_mean)

reveg_rich_mod <- glmmTMB(richness ~ visit + (1|site), data = Yr1Surv_year_2_3, family = poisson(link = "log"))
summary(reveg_rich_mod)
r2(reveg_rich_mod)
Yr1Surv_year_2_3$visit <- as_factor(Yr1Surv_year_2_3$visit)
surv_richBOX <- ggplot(
  data = Yr1Surv_year_2_3, aes(x=visit, y=richness, fill = visit)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Revegetation species richness") +
  scale_fill_brewer(palette="Set2", labels = c("Post planting", "Establishment"))+
  guides(colour=guide_legend(title="Visit"))+
  guides(fill=guide_legend(title="Visit"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
surv_richBOX
surv_richBOX <- surv_richBOX + theme(legend.position = "none")

# model it - new analyses (considers total surveyed area)
reveg_abun_mod <- glmmTMB(stemsperha ~ visit + (1|site), data = abun_stemsperha_reveg, family = poisson(link = "log"))
summary(reveg_abun_mod)
r2(reveg_abun_mod)
1/exp(-0.566653)

surv_abunBOX <- ggplot(
  data = beltdata_reveg, aes(x=visit, y=stemsperha, fill = visit)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Revegetation stems per ha") +
  scale_fill_brewer(palette="Set2", labels = c("Post planting", "Establishment"))+
  guides(colour=guide_legend(title="Visit"))+
  guides(fill=guide_legend(title="Visit"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
surv_abunBOX
surv_abunBOX <- surv_abunBOX + theme(legend.position = "none")

#### Q1/2b. 1m QUADRAT - Native Ground Cover species & other ground cover ----
# bring in site type into quadratdata dataframe
quadratdata24$sitetype <- sitedata$'Site_type'[match(quadratdata24$'WptID', sitedata$'WptID')]

# bring in site pair into beltdata dataframe
quadratdata24$pair <- sitedata$'Target_WptID'[match(quadratdata24$'WptID', sitedata$'WptID')]

# Rename the different coverage types
quadratdata <- quadratdata24 %>%
  rename(site = 'WptID',
         transectno = 'BeltID',
         beltpoint = 'Quadrat_ID',
         exotic = 'Exotic plants',
         native = 'Native plants',
         coarsewd = 'Course WD',
         finewd = 'Fine WD/litter',
         bareground = 'Bare Ground',
         matting = 'Matting/woodchips',
         rockcover = 'Rock cover',
         water = 'Water',
         macropod = 'Macropod_pellets',
         deer = 'Deer_pellets',
         rabbit = 'Rabbit_pellets',
         otherpel = 'Other_pellets',
         visit = 'Visit_number') 
quadratdata <- quadratdata %>% 
  dplyr::mutate(exotic = as.numeric(exotic), 
                native = as.numeric(native),
                coarsewd = as.numeric(coarsewd),
                finewd = as.numeric(finewd),
                matting = as.numeric(matting),
                rockcover = as.numeric(rockcover),
                water = as.numeric(water),
                deer = as.integer(deer),
                rabbit = as.integer(rabbit),
                macropod = as.integer(macropod),
                otherpel = as.integer(otherpel),)

# Only looking at visit 2 data
quadratdata$visit[quadratdata$sitetype == 'Remnant'] <- 3
quadratdata <- quadratdata %>% filter(visit == 3)
# Remove data with no WptID
quadratdata <- quadratdata %>% filter(!is.na(site))
# Remove data with no cover estimates
quadratdata <- quadratdata %>% filter(!is.na(exotic))
# Give pellet data 0 where nil is present
quadratdata <- quadratdata %>% 
  mutate_at(c('deer','rabbit','macropod','otherpel'), ~replace_na(.,0))
# summarise data
describe(quadratdata)
describeBy(quadratdata, quadratdata$sitetype)


# identify factors
quadratdata$site <- as.factor(quadratdata$site)
quadratdata$beltpoint <- as.factor(quadratdata$beltpoint)
quadratdata$sitetype <- as.factor(quadratdata$sitetype)


# identify variables
quadratdata$transectno <- as.numeric(quadratdata$transectno)
quadratdata$exotic <-as.numeric(quadratdata$exotic)
quadratdata$native <-as.numeric(quadratdata$native)
quadratdata$coarsewd <-as.numeric(quadratdata$coarsewd)
quadratdata$finewd <-as.numeric(quadratdata$finewd)
quadratdata$bareground <-as.numeric(quadratdata$bareground)
quadratdata$rockcover <-as.numeric(quadratdata$rockcover)
quadratdata$water <-as.numeric(quadratdata$water)


# summarise
quadratdatasummary <- quadratdata %>%
  group_by(sitetype, pair, site) %>%
  summarise(
    exotic = mean(exotic, na.rm = TRUE), 
    native = mean(native, na.rm = TRUE),
    coarsewd = mean(coarsewd,na.rm = TRUE ),
    finewd = mean(finewd,na.rm = TRUE ),
    bareground = mean(bareground, na.rm = TRUE ),
    rockcover = mean(rockcover, na.rm = TRUE ),
    water = mean(water, na.rm = TRUE)
  )
quadratdatasummary

# summarise data
describe(quadratdatasummary)
describeBy(quadratdatasummary, quadratdatasummary$sitetype)

# graph native ground cover percentages
nativegroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-native), y=native, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage native of ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

nativegroundcoverpercentage


#ggsave(nativegroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativegroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model native ground cover using glmm
quadratdata$nativebeta <- (quadratdata$native+1)/102
nativegroundcoverbysitetypeglmm <- glmmTMB(nativebeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(nativegroundcoverbysitetypeglmm)
r2(nativegroundcoverbysitetypeglmm)


# graph exotic ground cover percentages
exoticgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-exotic), y=exotic, fill = sitetype)) +
  geom_bar(stat="identity",colour = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Exotic groundcover (%)") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())

exoticgroundcoverpercentage

# ggsave(exoticgroundcoverpercentage, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/exoticgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 300)


# model exotic ground cover using glmm
quadratdata$exoticbeta <- (quadratdata$exotic+1)/102
exoticgroundcoverbysitetypeglmm <- glmmTMB(exoticbeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(exoticgroundcoverbysitetypeglmm)
r2(exoticgroundcoverbysitetypeglmm)


# graph finewd ground cover percentages
finewdgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-finewd), y=finewd, fill = sitetype)) +
  geom_bar(stat="identity",colour = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage of fine woody debris ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
finewdgroundcoverpercentage

#ggsave(finewdgroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/finewdgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model finewd ground cover using glmm
quadratdata$finewdbeta <- (quadratdata$finewd+1)/102
finewdgroundcoverbysitetypeglmm <- glmmTMB(finewdbeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(finewdgroundcoverbysitetypeglmm)


# graph coarsewd ground cover percentages
coarsewdgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-coarsewd), y=coarsewd, fill = sitetype)) +
  geom_bar(stat="identity",color = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage of coarse woody debris ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

coarsewdgroundcoverpercentage

# ggsave(coarsewdgroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/coarsewdgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model coarsewd ground cover using glmm 
quadratdata$coarsewdbeta <- (quadratdata$coarsewd+1)/102
coarsewdgroundcoverbysitetypeglmm <- glmmTMB(coarsewdbeta ~ sitetype + (1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(coarsewdgroundcoverbysitetypeglmm)


# graph  bareground ground cover percentages
baregroundgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-bareground), y=bareground, fill = sitetype)) +
  geom_bar(stat="identity",color = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage of bare ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
baregroundgroundcoverpercentage

# ggsave(baregroundgroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/baregroundgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model  bareground ground cover using glmm 
quadratdata$baregroundbeta <- (quadratdata$ bareground+1)/102
baregroundgroundcoverbysitetypeglmm <- glmmTMB(baregroundbeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(baregroundgroundcoverbysitetypeglmm)

#### Q1/2b 1m QUADRAT - COMPAIRING AVERAGE COVERAGE TYPE ----

#COMBINED GRAPH
# pivot longer
quadratdatasummarylong <- quadratdatasummary %>%
  pivot_longer(cols = exotic:bareground, names_to = "covertype", values_to = "cover")
quadratdatasummarylong$covertype <- factor(quadratdatasummarylong$covertype, levels = c("native", "exotic", "finewd", "coarsewd", "bareground"))

#summarise
quadratdatalongsummary <- quadratdatasummarylong %>%
  group_by(sitetype, covertype) %>%
  summarise(
    meancover = mean(cover),
    coverse = sd(cover)/sqrt(n()))

# graph it
covertypecomparisonfigure <- ggplot(
  data = quadratdatalongsummary, aes(x=covertype, y= meancover, fill = sitetype)) +
  geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin=meancover-coverse, ymax=meancover+coverse), position = position_dodge(0.9, preserve = "single"), width=0.1, size=0.5, color="black") +
  ylim(0, 100) +
  labs(x= NULL, y = "Understorey cover %") +
  scale_x_discrete(breaks=c("native","exotic","finewd", "coarsewd", "bareground"),
                   labels=c("Native plants","Exotic plants","Leaf litter", "Woody debris", "Bare ground")) +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site type", labels = c("Remnant", "Revegetated")))+
  guides(fill=guide_legend(title="Site type", labels = c("Remnant", "Revegetated")))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.ticks.x = element_blank())

covertypecomparisonfigure

#### Q1/2d. RECRUIT ABUNDANCE ----

# calculate recruits per belt area into count data
beltdatasummarynative$norecruitsperha <- round(beltdatasummarynative$norecruits/beltdatasummarynative$totalbeltarea*10000)
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)


# graph native and exotic species richness
recruitnativetreerichness <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site,-norecruits), y=norecruits, fill = sitetype)) +
  geom_bar(stat="identity", color = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Number of native woody recruits") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

recruitnativetreerichness

#ggsave(recruitnativetreerichness, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/recruitnativetreerichness.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
recruitrichnessbysitetypemodel <- glmmTMB(norecruitsperha ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(recruitrichnessbysitetypemodel)
r2(recruitrichnessbysitetypemodel)
1/exp(-0.88515)

# NATIVE RECRUITS BOXPLOT
nativerecruitsBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=norecruitsperha, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  expand_limits(y = 0) +
  labs(x = 'Site', y = "Native woody recruits (stems/ha)") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativerecruitsBOX
nativerecruitsBOX <- nativerecruitsBOX + theme(legend.position = "none")



#### Q1/3a. EXOTIC TREE AND SHRUB Richness----

# fill in missing cases for exotic taxa
# beltdatasummary <- beltdatasummary[complete.cases(beltdatasummary),]
#beltdatasummary <-  complete(beltdatasummary, origin, fill = list(richness= 0, nostems = 0, norecruits = 0))

# filter to only consider exotic
beltdatasummaryexotic <- filter(beltdatasummary, origin == "Introduced")
describeBy(beltdatasummaryexotic, beltdatasummaryexotic$sitetype)

# old figure of raw data
exotictreerichness <- ggplot(data = beltdatasummaryexotic, aes(x=reorder(site,-richness), y=richness, fill = sitetype)) +
  geom_bar(stat="identity", color ="black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Richness of exotic trees and shrubs") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

exotictreerichness

#ggsave(exotictreerichness, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/exotictreerichness.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
exoticrichnessbysitetypemodel <- glmmTMB(richness ~ sitetype +(1|pair), data = beltdatasummaryexotic, family = poisson(link = "log"))
summary(exoticrichnessbysitetypemodel)
r2(exoticrichnessbysitetypemodel)


# EXOTIC RICHNESS BOXPLOT
exotictreerichnessBOX <- ggplot(
  data = beltdatasummaryexotic, aes(x=sitetype, y=richness, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  expand_limits(y = 0) +
  labs(x = 'Site', y = "Exotic woody richness") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
exotictreerichnessBOX
exotictreerichnessBOX <- exotictreerichnessBOX + theme(legend.position = "none")

#### Q1/3a. EXOTIC TREE AND SHRUB ABUNDANCE ----

# incorporate belt area into exotic count data
beltdatasummaryexotic$totalbeltarea <- beltlength$belt_length[match(beltdatasummaryexotic$site, beltlength$site)]
beltdatasummaryexotic$stemsperha <- round(beltdatasummaryexotic$nostems/beltdatasummarynative$totalbeltarea*10000)
describeBy(beltdatasummaryexotic, beltdatasummaryexotic$sitetype)


# graph of raw data
exotictreeabundance <- ggplot(
  data =  beltdatasummaryexotic, aes(x=reorder(site,-nostems), y=nostems, fill = sitetype)) +
  geom_bar(stat="identity", color="black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Abundance of exotic trees and shrubs") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

exotictreeabundance


# ggsave(exotictreeabundance, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/exotictreeabundance.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model it
exoticabundancebysitetypemodel <- glmmTMB(nostems ~ sitetype + (1|pair), data = beltdatasummaryexotic, family = poisson(link = "log"))
summary(exoticabundancebysitetypemodel)
r2(exoticabundancebysitetypemodel)
exp(0.7529)

# EXOTIC Woody plant abundance BOXPLOT
exotictreeabundanceBOX <- ggplot(
  data = beltdatasummaryexotic, aes(x=sitetype, y=stemsperha, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  expand_limits(y = 0) +
  labs(x = 'Site', y = "Exotic woody density") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
exotictreeabundanceBOX
legend <- get_legend(exotictreeabundanceBOX)
exotictreeabundanceBOX <- exotictreeabundanceBOX + theme(legend.position = "none")


# combine belt transect data boxplots
boxmedley <- ggarrange(nativetreerichnessBOX, exotictreerichnessBOX, legend,  nativetreeabundanceBOX, exotictreeabundanceBOX, nativerecruitsBOX, align = "v", ncol = 3, nrow = 2, labels = c("A", "B", "", "C", "D", "E"))
boxmedley <- ggarrange(nativetreerichnessBOX, nativetreeabundanceBOX, nativerecruitsBOX, exotictreerichnessBOX, exotictreeabundanceBOX, legend, align = "v", ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "")) # different order
boxmedley

# ggsave(boxmedley, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/boxplotmedley.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)
# ggsave(boxmedley, filename = "boxplotmedley.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)


# ggsave(covertypecomparisonfigure, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/understoreycovers.tiff", width = 210, height = 120, units = c("mm"), bg = "white", dpi = 300)


#### Q1/2c STRATA HEIGHTS - CANOPY COVER VEGETATION HEIGHT DIVERSITY ----

#  CANOPY COVER 
# lets look at the canopy first
# bring in site type into canopydata dataframe
canopydata$sitetype <- sitedata$'Site type'[match(canopydata$'Site Code', sitedata$'Site code')]

# bring in site pair into beltdata dataframe
canopydata$pair <- sitedata$'Paired site code'[match(canopydata$'Site Name', sitedata$'Site Name')]

# rename 'messy' factor names
canopydata <- canopydata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         beltm = `Belt intercept`,
         canopy = 'Canopy')


# identify factors
canopydata$sitetype <- as.factor(canopydata$sitetype)
canopydata$pair <- as.factor(canopydata$pair)
canopydata$site <- as.factor(canopydata$site)
canopydata$transect <- as.factor(canopydata$transect)
canopydata$beltm <- as.factor(canopydata$beltm)
canopydata$canopy <- as.factor(canopydata$canopy)
canopydata$canopy <- factor(canopydata$canopy,levels=c("P","A"),labels=c('presence','absence'))

# exclude transects rows with NAs
canopydata <- canopydata[complete.cases(canopydata), ]


# create frequency table of canopy
dffreq <- canopydata %>%
  count(sitetype, pair, site, canopy) %>%
  group_by(site) %>%          
  mutate(prop = prop.table(n))
dffreq <- as.data.frame(dffreq)


# consider only presence
dffreqpres <- filter(dffreq, canopy == "presence")


# freq figure for presence of canopy by siteypte NB. Definitely probably worth separating out East and West sites.
frequencyofcanopy <- ggplot(
  data = dffreqpres, aes(x=reorder(site,-prop), y= prop*100, fill = sitetype)) +   geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  ylim(0, 100) +
  labs(x="Site", y = "Canopy cover %")+
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

frequencyofcanopy

# ggsave(frequencyofcanopy, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/frequencyofcanopy.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
canopybysitetypeglmm <- glmmTMB(prop ~ sitetype + (1|pair), data = dffreqpres, family=beta_family(link="logit"))
summary(canopybysitetypeglmm)

describe(dffreqpres)
describeBy(dffreqpres, dffreqpres$sitetype)


# Canopy presence BOXPLOT
canopyBOX <- ggplot(
  data = dffreqpres, aes(x=sitetype, y=prop, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  ylim(0,1) +
  labs(x = 'Site type', y = "Canopy presence %") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
canopyBOX
# legend <- get_legend(canopyBOX)
# canopyBOX <- canopyBOX + theme(legend.position = "none")

# combine with boxmedley to create Figure 3 for paper
Fig3 <- ggarrange(nativetreerichnessBOX, nativetreeabundanceBOX,nativerecruitsBOX, exotictreerichnessBOX, exotictreeabundanceBOX, canopyBOX, common.legend = TRUE, legend = "bottom", align = "v", ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F")) # different order
Fig3

# ggsave(Fig3, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/Fig3.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)


### ASSESS VEGETATION STRUCTURE - TREE HEIGHT AND SIZE ----

# bring in site type
heightdata$sitetype <- sitedata$'Site type'[match(heightdata$'Site Code', sitedata$'Site code')]

# bring in site pair into beltdata dataframe
heightdata$pair <- sitedata$'Paired site code'[match(heightdata$'Site Name', sitedata$'Site Name')]


# rename 'messy' factor names
heightdata <- heightdata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         species = `Species name`,
         height = 'Heights (first 5 of each Sp to 0.5m)')

head(heightdata)
species <- levels(as.factor(heightdata$species))
species


# identify factors
heightdata$sitetype <- as.factor(heightdata$sitetype)
heightdata$pair <- as.factor(heightdata$pair)
heightdata$site <- as.factor(heightdata$site)
heightdata$transect <- as.factor(heightdata$transect)
heightdata$species <- as.factor(heightdata$species)


# identify variables
heightdata$height <- as.numeric(heightdata$height)

head(heightdata)


# graph raw heights using boxplots?
heightfig <- ggplot(
  data = heightdata, aes(x=species, y= height, fill = sitetype)) +
  geom_boxplot() +
  labs(x= NULL, y = "Height (m)") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfig

# there are LOTS of different species... subset?

# graph raw heights using boxplots at sitetype level only?
heightfig2 <- ggplot(
  data = heightdata, aes(x=sitetype, y= height, fill = sitetype)) +
  geom_boxplot() +
  labs(x= NULL, y = "Height (m)") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfig2


# summarise heightdata
heightdatasummary <- heightdata %>%
  group_by(sitetype, pair, site, species) %>%
  summarise(
    count = n(),
    meanheight = mean(height))

describeBy(heightdatasummary, heightdatasummary$sitetype)
# new <- heightdatasummary[which(heightdatasummary$count >= 3),]


# filter out species for which <3 plants were surveyed at a site?
heightdatasummarytrunc <- filter(heightdatasummary, count > 2)

# at how many sites were those species found? 
speciesbysite <- heightdatasummarytrunc %>%
  group_by(sitetype, species) %>%
  summarise(
    count = n(),
    height = mean(meanheight))


# Matt's fancy code to select only species which occur at at least three site of a given sitetype
new <- speciesbysite[which(speciesbysite$count > 2),]

want <- unique(as.character(new$species))

test <- NULL
for(i in 1:NROW(want)) {
  
  get <- heightdata[which(as.character(heightdata$species) == want[i]),]
  test <- rbind(test, get) }

test <- filter(test, species != 'Eucalyptus ovata')


# graph raw heights using boxplots using truncated data?
heightfigtrunc <- ggplot(
  data = test, aes(x=species, y= height, fill = sitetype)) +
  geom_boxplot(position = position_dodge(preserve = 'single')) +
  labs(x= "Species", y = "Plant height (m)") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(fill=guide_legend(title="Site type"))+
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfigtrunc

# ggsave(heightfigtrunc, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/Fig S1.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)




# graph raw heights using boxplots at sitetype level only?
heightfigtrunc2 <- ggplot(
  data = test, aes(x=sitetype, y= height, fill = sitetype)) +
  geom_boxplot() +
  labs(x= NULL, y = "Height (m)") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfigtrunc2


# Model of plant heights by species for truncated data
heightslmm <- glmmTMB(height ~ sitetype*species + (1|site/pair), data = test, family=gaussian())
summary(heightslmm)
plot_model(heightslmm)




### STRATA DATA - ASSESS VEGETATION STRUCTURE - VEGETATION STRUCTURE ----
# load data
sitedata <- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Site info")
stratadata<- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Strata Heights")

# lets look at the canopy first
# bring in site type into strata dataframe
stratadata$sitetype <- sitedata$'Site type'[match(stratadata$'Site Code', sitedata$'Site code')]

# bring in site pair into strata dataframe
stratadata$pair <- sitedata$'Paired site code'[match(stratadata$'Site Name', sitedata$'Site Name')]


# rename 'messy' factor names
stratadata <- stratadata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         beltm = `Belt intercept`,
         height = 'Poll height')


# identify factors
stratadata$sitetype <- as.factor(stratadata$sitetype)
stratadata$pair <- as.factor(stratadata$pair)
stratadata$site <- as.factor(stratadata$site)
stratadata$transect <- as.factor(stratadata$transect)
stratadata$beltm <- as.factor(stratadata$beltm)
stratadata$height <- as.factor(stratadata$height)
canopydata$canopy <- factor(canopydata$canopy,levels=c("P","A"),labels=c('presence','absemce'))

# identify variables
stratadata$Tree <- as.numeric(stratadata$Tree)
stratadata$Shrub <- as.numeric(stratadata$Shrub)
stratadata$Fern <- as.numeric(stratadata$Fern)
stratadata$Herbaceous <- as.numeric(stratadata$Herbaceous)

# exclude transects rows with NAs
stratadata <- stratadata[complete.cases(stratadata), ] # OR
stratadata <- na.omit(stratadata)

# create count able of different plant types
stratacount <- stratadata %>%
  group_by(sitetype, pair, site, height) %>%          
  summarise(
    count = n(),
    tree = sum(Tree),
    shrub = sum(Shrub),
    fern = sum(Fern),
    herb = sum(Herbaceous))

stratacountworks <- filter(stratacount, sitetype == "Works")


describe(stratacountworks)

# create frequency table of different plant types
stratasum <- stratadata %>%
  group_by(sitetype, height) %>%          
  summarise(
    count = n(),
    tree = sum(Tree)/n(),
    shrub = sum(Shrub)/n(),
    fern = sum(Fern)/n(),
    herb = sum(Herbaceous)/n())


# convert to long form for graphing
stratasumlong <- stratasum %>%
  pivot_longer(cols = tree:herb, names_to = "planttype", values_to = "frequency")


# Figure for frequency of different plant types at different heights by sitetype


# OPTION 1
stratavegtypecomparison1 <- ggplot(
  data = stratasumlong, aes(x=height, y= frequency, fill = sitetype)) +
  geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  ylim(0, 1) +
  labs(x= NULL, y = "Frequency") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~planttype)

stratavegtypecomparison1

# ggsave(stratavegtypecomparison1, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/stratavegtypecomparison1.tiff", width = 16, height = 12, units = "cm", dpi = 600)



# OPTION 2 #### THE ONE I"M USING # NB. NO ERROR BARS! Need to create df from site summary (which is below) to create errorbars
stratasitesypecomparison <- ggplot(
  data = stratasumlong, aes(x=height, y= frequency, fill = planttype)) +
  geom_bar(stat = "identity", aes(fill = planttype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  ylim(0, 1) +
  labs(x= NULL, y = "Frequency") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Plant form"))+
  guides(fill=guide_legend(title="Plant form"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~sitetype, ncol = 1)

stratasitesypecomparison

#ggsave(stratasitesypecomparison, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/stratasitesypecomparison.tiff", width = 16, height = 12, units = "cm", dpi = 600)

stratacomparesitetypefig <- ggplot(data =stratasumlong, aes(y=height, x=planttype, fill=frequency, size=frequency, colour = frequency)) +
  geom_point(alpha=0.5, shape=21) +
  scale_size(range = c(1, 30), name="Frequency of hits") +
  facet_wrap(.~sitetype) +
  theme_classic()+
  labs(x= "Plant type", y="Pole Heights")
stratacomparesitetypefig


stratacomparesitetypefigurecircle <- ggplot(data =stratasumlong, aes(y=height, x=planttype, fill=frequency, size=frequency)) +
  geom_point(alpha=1, shape=21, color="black") +
  facet_wrap(.~sitetype) +
  labs(x= "Plant type", y="Pole Heights") +
  scale_fill_continuous(name = "Frequency",  limits=c(0, 1), breaks=seq(0, 1, by=0.2)) + 
  scale_size_continuous(name = "Frequency", range = c(1, 15), limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  scale_fill_steps2(
    low = "darkorange", 
    mid = "purple", 
    high = "darkgreen", 
    midpoint = 0.6)+
  guides(fill=guide_colourbar(reverse=TRUE)) + 
  guides(fill=guide_legend(), size = guide_legend()) +
  theme_classic()

stratacomparesitetypefigurecircle

#ggsave(stratacomparesitetypefigurecircle, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/stratacomparesitetypefigurecircle.tiff", width = 16, height = 12, units = "cm", dpi = 600)

#summarise 
stratasumbysite <- stratadata %>%
  group_by(sitetype, pair, site, height) %>%    #  NEED TO ADD SITE AS GROUPING FACTOR HERE     
  summarise(
    count = n(),
    tree = sum(Tree)/n(),
    shrub = sum(Shrub)/n(),
    fern = sum(Fern)/n(),
    herb = sum(Herbaceous)/n())


# convert to long form for graphing
stratasumbysitelong <- stratasumbysite %>%
  pivot_longer(cols = tree:herb, names_to = "planttype", values_to = "frequency")


#summarise by sitetype and planttype
stratasumbysitelongsum <- stratasumbysitelong %>%
  group_by(sitetype, height, planttype) %>%    #  NEED TO ADD SITE AS GROUPING FACTOR HERE     
  summarise(
    meanfreq = mean(frequency),
    sefreq = sd(frequency)/sqrt(n()))


# FIGURE FOR PAPER
stratasumbysitelongsum$planttype <- factor(stratasumbysitelongsum$planttype, levels=c("tree","shrub","herb","fern"),labels=c("Tree","Shrub","Herb","Fern"))

stratavegtypesitetypecomparison <- ggplot(
  data = stratasumbysitelongsum, aes(x=height, y= meanfreq, fill = sitetype)) +
  geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin=meanfreq-sefreq, ymax=meanfreq+sefreq), position = position_dodge(0.9, preserve = "single"), width=0.1, size=0.5, color="black") +
  ylim(0, 1) +
  labs(x= "Vegetation height (m)", y = "Frequency of hits") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site type"))+
  guides(fill=guide_legend(title="Site type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~planttype) +
  theme(strip.text.x = element_text(size = 14))
stratavegtypesitetypecomparison 
stratavegtypesitetypecomparison <- stratavegtypesitetypecomparison +coord_flip()
stratavegtypesitetypecomparison

# ggsave(stratavegtypesitetypecomparison, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/stratavegtypesitetypecomparison.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)


## STRATA MODELS
# if you need summary stats adapt this
stratasumbysitelongsumtemp <- filter(stratasumbysitelongsum, planttype == "Fern")
describeBy(stratasumbysitelongsumtemp, stratasumbysitelongsumtemp$sitetype)


# model it - tree frequency by height and site type
stratasumbysite$treebeta <- (stratasumbysite$tree+.01)
stratatreehitfrequencyglmm <- glmmTMB(treebeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(stratatreehitfrequencyglmm)
Anova(stratatreehitfrequencyglmm)

emmip(stratatreehitfrequencyglmm, height ~ sitetype)

#model it - fern frequency by height and site type
stratasumbysite$fernbeta <- (stratasumbysite$fern+.01)
stratafernhitfrequencyglmm <- glmmTMB(fernbeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(stratafernhitfrequencyglmm)
Anova(stratafernhitfrequencyglmm)

emmip(stratafernhitfrequencyglmm, height ~ sitetype)

stratasumbysite$shrubbeta <- (stratasumbysite$shrub+.01)
stratashrubhitfrequencyglmm <- glmmTMB(shrubbeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(stratashrubhitfrequencyglmm)
Anova(stratashrubhitfrequencyglmm)

emmip(stratashrubhitfrequencyglmm, height ~ sitetype)

stratasumbysite$herbbeta <- (stratasumbysite$herb+.01)/1.02
strataherbhitfrequencyglmm <- glmmTMB(herbbeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(strataherbhitfrequencyglmm)
Anova(strataherbhitfrequencyglmm)

emmip(strataherbhitfrequencyglmm, height ~ sitetype)

#### Q3 LANDSCAPE CONTEXT - NATIVE WOODY VEGETATION COVERAGE - RECRUIT ABUNDANCE ---- 

# consider only works sites
beltdatasummarynativeworks <- filter(beltdatasummarynative, sitetype == "Works")

# graph native and exotic species richness
nativewoodyvegcoveragerecruits <- ggplot(data = beltdatasummarynativeworks, aes(x=nativevegcoveragepercentage*100, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  geom_smooth (method = "glm", method.args = list(family = "poisson"),
               colour = "darkblue", se = TRUE) +
  labs(x = 'Remnant vegetation %', y = "Native woody recruits (stems/ha)") +
  annotate("text", x = 60, y = 3000, label = "italic(R) ^ 2 == 0.01",
           parse = TRUE) +
  theme_classic()
nativewoodyvegcoveragerecruits

#ggsave(nativewoodyvegcoveragerecruits, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativewoodyvegcoveragerecruits.tiff", width = 16, height = 12, units = "cm", dpi = 600)

nativevegcovervrecruits <- glm(norecruitsperha ~ nativevegcoveragepercentage, data = beltdatasummarynativeworks,   family = poisson(link = "log"))
summary(nativevegcovervrecruits)
r2(nativevegcovervrecruits)
with(summary(nativevegcovervrecruits), 1 - deviance/null.deviance)
model_performance(nativevegcovervrecruits)

#### Q3. BROWSER POO FREQUENCY ----

# load data
poodata <- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Pellet counts")

# bring in site type into poodata dataframe
poodata$sitetype <- sitedata$'Site type'[match(poodata$'Site Code', sitedata$'Site code')]

# bring in site pair into beltdata dataframe
poodata$pair <- sitedata$'Paired site code'[match(poodata$'Site Code', sitedata$'Site code')]

# rename 'messy' factor names
poodata <- poodata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         belt = `Belt point`)


# identify factors
poodata$sitetype <- as.factor(poodata$sitetype)
poodata$pair <- as.factor(poodata$pair)
poodata$site <- as.factor(poodata$site)
poodata$transect <- as.factor(poodata$transect)
poodata$belt <- as.factor(poodata$belt)


# replace P's and A's with 1s and 0s
poodata <- poodata[!(poodata$Deer == "NA"),] # remove one row with NAs
poodata[poodata  == "P"] <- "1"
poodata[poodata  == "A"] <- "0"
poodata$Kangaroo <- as.numeric(poodata$Kangaroo)
poodata$Wallaby <- as.numeric(poodata$Wallaby)
poodata$Rabbit <- as.numeric(poodata$Rabbit)
poodata$Deer <- as.numeric(poodata$Deer)

# create count able of browser poo counts
poodatacount <- poodata %>%
  group_by(sitetype, pair, site) %>%          
  summarise(
    count = n(),
    kanga = sum(Kangaroo),
    wall = sum(Wallaby),
    rabbit = sum(Rabbit),
    deer = sum(Deer))

# create table of browser poo frequency # NOTE DF NAME SAME AS ABOVE 
poodatacount <- poodata %>%
  group_by(sitetype, pair, site) %>%          
  summarise(
    count = n(),
    kanga = sum(Kangaroo)/count,
    wall = sum(Wallaby)/count,
    rabbit = sum(Rabbit)/count,
    deer = sum(Deer)/count)


# sum across all browser types
poodatacount$total <- poodatacount$kanga + poodatacount$wall + poodatacount$rabbit + poodatacount$deer 
describe(poodatacount)

# graph browser frequency 
browsertotal <- ggplot(
  data = poodatacount, aes(x=reorder(site,-total), y=total, fill = sitetype)) +
  geom_bar(stat="identity", colour = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Frequency of browser pellets") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
browsertotal

#ggsave(browsertotal, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/browsertotale.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# merge poodata and beltdata
beltdatasummarynative$browsers <- poodatacount$total[match(beltdatasummarynative$site, poodatacount$site)]  

beltdatasummarynativeworks <- filter(beltdatasummarynative, sitetype == "Works")

browservrecruit <- ggplot(data = beltdatasummarynativeworks, aes(x=browsers, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  geom_smooth (method = "glm",method.args = list(family = "poisson"),
               colour = "darkblue") +
  labs(x = 'Browser pellet frequency ', y = "Native woody recruits stems/ha)") +
  annotate("text", x = 0.2, y = 3000, label = "italic(R) ^ 2 == 0.25",
           parse = TRUE) +
  theme_classic()
browservrecruit

#ggsave(browservrecruit, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/browservrecruit.tiff", width = 16, height = 12, units = "cm", dpi = 600)


browservrecruitmodel <- glm(norecruitsperha ~ browsers, data = beltdatasummarynativeworks, family = poisson(link = "log"))
summary(browservrecruitmodel)
r2(browservrecruitmodel)
with(summary(browservrecruitmodel), 1 - deviance/null.deviance)
model_performance(browservrecruitmodel)

count_belt2 <- count_belt %>% group_by(WptID, Belt_ID)%>% count()
beltdata %>% group_by(WptID, visit, Belt_ID) %>% count(WptID)
], Belt_ID, visit) %>%count(WptID, sort = TRUE)

# bring in site pair into beltdata dataframe
#beltdata$visit <- visitdata$'Visit_number'[match(beltdata$'VisitID', visitdata$'VisitID')]

# BRING IN TRANSECT LENGTH DATA ??? 

#### Calculate the number native tree and shrub species per site
beltdata2 <- beltdata %>%
  rename(site = 'WptID',
         species = 'SCI_NAME',
         origin = 'ORIGIN',
         count = 'Count',
         recruits = 'Type',
         alive = 'Count_alive',
         visit = 'Visit_number')
#        nativevegcoveragecounts = 'native veg coverage counts',
#       nativevegcoveragepercentage = 'Woody_perc') 

# identify factors
beltdata2$site <- as.factor(beltdata2$site)
beltdata2$species <- as.factor(beltdata2$species)
beltdata2$origin <- as.factor(beltdata2$origin)
beltdata2$sitetype <- as.factor(beltdata2$sitetype)
beltdata2$pair <- as.factor(beltdata2$pair)
beltdata2$visit <- as.factor(beltdata2$visit)

# identify variables
beltdata2$count <- as.numeric(beltdata2$count)
#beltdata$recruits <-as.numeric(beltdata$recruits)

# data exploration - belts
glimpse(beltdata2)
describe(beltdata2)
describeBy(beltdata2, beltdata2$sitetype)

# species level summaries
# summary totals - propagules - taxa
sort(unique(beltdata2$species)) # SOME OF THESE ARE NOT WOODY SPECIES!!!
beltdata2 <- filter(beltdata2, species != "Atriplex spp.")
beltdata2 <- filter(beltdata2, species != "Senecio spp.")
beltdata2 <- filter(beltdata2, species != "Atriplex semibaccata")
beltdata2 <- filter(beltdata2, species != "Atriplex cinerea")
beltdata2 <- filter(beltdata2, species != "Carpha spp.")
beltdata2 <- filter(beltdata2, species != "Mentha australis")


# CHECK FOR OTHERS

with(beltdata2, tapply(species, list(sitetype), FUN = function(x) length(unique(x))))
with(beltdata2, tapply(species, list(site), FUN = function(x) length(unique(x))))
with(beltdata2, tapply(species, list(sitetype, origin), FUN = function(x) length(unique(x))))
with(beltdata2, tapply(species, list(origin), FUN = function(x) length(unique(x))))

#survival
survival <- filter(beltdata2, visit %in% c(2,3), recruits %in% c("reveg"))
survival1 <- survival %>% group_by (site) %>% filter (visit==3)
survival2 <- filter(survival, site %in% c(2, 3, 4, 6, 7, 8, 9, 10, 12))
survival2 <- filter(survival2, origin != 'Introduced')

sp_abun1 <- survival2 %>% 
  group_by(site, VisitID, Belt_ID, visit) %>%
  summarise_at(vars(count, alive), tibble::lst(mean, median))
summaryabun <- sp_abun1 %>% mutate(percent=alive_mean/count_mean)
summaryabun

sp_abun1 <- summaryabun %>% 
  group_by(visit) %>%
  mutate(Count = n()) %>%
  summarise_at(vars(percent), funs(mean, sd, se=sd(.)/sqrt(n())))
sp_abun1
periodlabs <- c("Planting", "Post-planting")

fig1 <- ggplot(sp_abun1, aes(x = reorder(visit, -mean), y = mean)) +
  geom_bar(position=position_dodge(), stat="identity", fill="lightblue") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.2,
                position=position_dodge(.9)) +
  labs(y = "Percent survival per site", x = "Survey Period") +
  theme_classic() + theme(legend.position = c(0.8, 0.8)) +
  scale_x_discrete(labels = periodlabs)+
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 11), axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=11), legend.title=element_text(size=11), 
        legend.text=element_text(size=11))+
  ylim(0,1)
fig1


# create summary of species per sitetype by no of sites for paper
remreveg <- filter(beltdata2, origin != 'Introduced', visit %in% c(2))
speciessummarybysite <- remreveg %>%
  group_by(species, origin, sitetype, site, visit) %>%
  summarise(
    noofstems = sum(count))

speciessummary <- speciessummarybysite %>%
  group_by(species, origin) %>%
  summarise(
    noofremsites = sum(sitetype == "Remnant"),
    noofworkssites = sum(sitetype == "Revegetation"),
    meannoofstemsrem = mean(noofstems[sitetype == "Remnant"], na.rm = TRUE),
    meannoofstemsworks = mean(noofstems[sitetype == "Revegetation"], na.rm = TRUE))

speciessummary <- speciessummary %>%
  mutate(sitedifremminusworks = noofremsites - noofworkssites)

speciessummary[is.na(speciessummary)] <- 0
# write.csv(speciessummary, "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/specieslistforTableS1.csv", row.names=FALSE)

# summarise
beltdatasummary <- remreveg %>%
  group_by(sitetype, pair, site, visit, origin) %>%
  summarise(
    richness = n_distinct(species),
    nostems = sum(count))
#   norecruits = sum(recruits))
# summarise

beltdatasurvive <- remreveg %>%
  group_by(sitetype, pair, site, visit, origin) %>%
  summarise(
    richness = n_distinct(species),
    surv = (count/alive))
#   norecruits = sum(recruits))
beltsummarycomplete <-  complete(beltdatasummary, origin, fill = list(richness= 0, nostems = 0))

# filter to consider only native
beltdatasummarynative <- filter(beltdatasummary, origin == "Native")
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)

# graph native species richness
nativetreerichness <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site, -richness), y=richness, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Richness of native trees and shrubs") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreerichness

# model it
nativerichnessbysitetypemodel <- glmmTMB(richness ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativerichnessbysitetypemodel)
r2(nativerichnessbysitetypemodel)

# BOX PLOT OF NATIVE RICHNESS

nativetreerichnessBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=richness, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Native woody richness") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativetreerichnessBOX

############################ Individual Species Abun ##########################################
Indiv_data <- remreveg %>%
  group_by(site, species) %>%
  summarise_at(vars(count, alive), tibble::lst(mean, median))

#calculate the mean and CI for alive and total
spdata_sp <- Indiv_data %>%
  group_by(species) %>%
  summarise(n = n(), 
            meantotal = mean(count_mean),
            setotal = sd(count_mean)/sqrt(n()),
            meanalive = mean(alive_mean),
            sealive = sd(alive_mean)/sqrt(n()))
Rich <- spdata_sp %>% mutate(surv = meanalive / meantotal)
#spdata_sp <- Rich %>% arrange(desc(-surv))
Species_data <- spdata_sp %>% arrange(desc(n))

# cut the first 20 rows
Data_sp16 <-Species_data[1:10,]
tail(Data_sp16)
#Data_sp16 <- dplyr::select (Data_sp16,-c(n))
Data_sp16 <- dplyr::select (Data_sp16,-c(surv))
#which(is.na(Sp_reshapeses))

# reshaping of data (this is pretty messy sorry)
Sp_reshape <- gather(Data_sp16, key = "variable", value = "mean", -c(species, n, setotal, sealive))
Sp_reshape <- subset(Sp_reshape, select = -c(setotal, sealive))
Sp_reshapeses <- gather(Data_sp16, key = "variable", value = "error", -c(species, n, meantotal, meanalive))
Sp_reshape$error <- Sp_reshapeses$error
Sp_reshape <- as.data.frame(Sp_reshape)

## Question - how do I include CI in this plot
Sp_reshape$variable <- factor(Sp_reshape$variable, levels = c("meantotal", "meanalive"))
#Sp_reshape$lifeform <- as.factor(Sp_reshape$lifeform)

#reorder by lifeform, variable and mean (heighest to lowest)
Sp_reshape <- Sp_reshape[with(Sp_reshape, order(variable, -mean)),]
Sp_reshape$species <- factor(Sp_reshape$species, levels = unique(Sp_reshape$species))
#ggsave(path = path, width = width, height = height, device='jpeg', dpi=300)

fig3 <- ggplot(Sp_reshape, aes(x = species, y = mean, fill = variable)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-error, ymax=mean+error), 
                width=.2,
                position=position_dodge(.9)) +
  #facet_grid(~ lifeform, scales = "free_x", space = "free") +
  scale_fill_grey(start = 0.8, end = 0.6, name="Individual plant counts", breaks=c("meantotal", "meanalive"),
                  labels=c("Total plants", "Alive plants")) + 
  #scale_x_discrete(labels=c("Ozothamnus ferrugineus" = "Ozo. ferrugineus", "Leptospermum continentale" = "Lepto. continentale",
  #                          "Leptospermum lanigerum" = "Lepto. lanigerum", "Cassinia aculeata subsp. aculeata" = "Cassinia aculeata", 
  #                          "Allocasuarina verticillata" = "Allo. verticillata", 
  #                          "Acacia verniciflua s.l." = "Acacia verniciflua")) +
  labs(y = "Mean abundance per site", x = "Species") +
  theme_classic() + theme(legend.position = c(0.7, 0.8)) +
  theme(axis.text.x = element_text(size = 9, angle = 55, hjust = 1),
        axis.text.y = element_text(size = 9), axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9), legend.title=element_text(size=9), 
        legend.text=element_text(size=9))
fig3

new <- beltdata[which(beltdata$visit_num >= 2),]

new #had 2 and 3ss only
new$site_num <- as.numeric(as.character(new$site))

data_int <- NULL
for(i in 1:NROW(levels(as.factor(new$site_num)))) {
  rows_site_i <- which(new$site_num == as.numeric(levels(as.factor(new$site_num)))[i])
  if(sum(unique(new$visit_num[rows_site_i])) >= 5) { p <- new[rows_site_i,] } else {}
  data_int <- rbind(data_int, p) }


data_int2 <-data_int %>% 
  group_by (site, species, VisitID, Belt_ID, recruits, counts, Avg_Height, visit_num)









# create summary of species per sitetype by no of sites for paper
remreveg <- filter(beltdata, origin != 'Introduced')
survival2 <- filter(survival, site %in% c(2, 3, 4, 6, 7, 8, 9, 10, 12))
speciessummarybysite <- remreveg %>%
  group_by(species, origin, sitetype, site, visit) %>%
  summarise(
    noofstems = sum(counts))

speciessummary <- speciessummarybysite %>%
  group_by(species, origin) %>%
  summarise(
    noofremsites = sum(sitetype == "Remnant"),
    noofworkssites = sum(sitetype == "Revegetation"),
    meannoofstemsrem = mean(noofstems[sitetype == "Remnant"], na.rm = TRUE),
    meannoofstemsworks = mean(noofstems[sitetype == "Revegetation"], na.rm = TRUE))

speciessummary <- speciessummary %>%
  mutate(sitedifremminusworks = noofremsites - noofworkssites)

speciessummary[is.na(speciessummary)] <- 0
# write.csv(speciessummary, "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/specieslistforTableS1.csv", row.names=FALSE)

# summarise
remreveg <- filter(beltdata, origin != 'Introduced', recruits %in% c("reveg"))
beltdatasummary <- remreveg %>%
  group_by(sitetype, pair, site, visit, origin) %>%
  summarise(
    richness = n_distinct(species),
    nostems = sum(counts))
#   norecruits = sum(recruits))
# summarise

beltdatasurvive <- beltdata %>%
  group_by(sitetype, pair, site, visit, origin) %>%
  summarise(
    richness = n_distinct(species),
    surv = (counts/Count_alive))
#   norecruits = sum(recruits))
beltsummarycomplete <-  complete(beltdatasummary, origin, fill = list(richness= 0, nostems = 0))

# filter to consider only native
beltdatasummarynative <- filter(beltdatasummary, origin == "Native")
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)

# graph native species richness
nativetreerichness <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site, -richness), y=richness, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Richness of native trees and shrubs") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreerichness


#nativetreerichness <- nativetreerichness + theme(legend.position = "none")

#boxplot species richness, survival, plant density accross sites/surveys

#Comparing reveg to remnant at year 3
# Comparing survival between visit 2 and visit 3 and how env variables effect that
# Composition between remnant and reveg areas after planting and after first year
#ggsave(nativetreerichness, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativetreerichness.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
nativerichnessbysitetypemodel <- glmmTMB(richness ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativerichnessbysitetypemodel)
r2(nativerichnessbysitetypemodel)

# BOX PLOT OF NATIVE RICHNESS

nativetreerichnessBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=nostems, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Native woody richness") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativetreerichnessBOX
# nativetreerichnessBOX <- nativetreerichnessBOX + theme(legend.position = "none")


#### Q1/2a. NATIVE TREE AND SHRUB ABUNDANCE ----

nativetreeabundance <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site,-nostems), y=nostems, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Abundance of native trees and shrubs") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreeabundance

# but need to consider TRANSECT length!!!

# rename belt transect length
sitedata2 <- sitedata %>%
  rename(beltlength = 'Plot Length (m)',
         site = 'Site Name') 

# calculate average transect length by site
beltlengthdata <- sitedata2 %>%
  group_by(site) %>%
  summarise(
    nbelts = n(),
    sumbeltlength = sum(beltlength),
    totalbeltarea = sumbeltlength*nbelts)


# incorporate belt area into count data
beltdatasummarynative$totalbeltarea <- beltlengthdata$totalbeltarea[match(beltdatasummarynative$site, beltlengthdata$site)]
beltdatasummarynative$stemsperha <- round(beltdatasummarynative$nostems/beltdatasummarynative$totalbeltarea*10000)
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)



# plot with new data
nativetreeabundance <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site,-stemsperha), y=stemsperha, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Native tree and shrub density (stems/ha)") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
nativetreeabundance


# model it - old analyses
nativeabundancebysitetypemodel <- glm(nostems ~ sitetype, data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativeabundancebysitetypemodel)
r2(nativeabundancebysitetypemodel)


# model it - new analyses (considers total surveyed area)
nativeabundanceperhabysitetypemodel <- glmmTMB(stemsperha ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(nativeabundanceperhabysitetypemodel)
r2(nativeabundanceperhabysitetypemodel)
1/exp(-0.566653)

# BOX PLOT OF NATIVE WOODY PLANT ABUNDANCE

nativetreeabundanceBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=stemsperha, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  labs(x = 'Site', y = "Native woody plant density") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  expand_limits(y = 0) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativetreeabundanceBOX
nativetreeabundanceBOX <- nativetreeabundanceBOX + theme(legend.position = "none")



#### Q1/2d. RECRUIT ABUNDANCE ----

# calculate recruits per belt area into count data
beltdatasummarynative$norecruitsperha <- round(beltdatasummarynative$norecruits/beltdatasummarynative$totalbeltarea*10000)
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)


# graph native and exotic species richness
recruitnativetreerichness <- ggplot(
  data = beltdatasummarynative, aes(x=reorder(site,-norecruits), y=norecruits, fill = sitetype)) +
  geom_bar(stat="identity", color = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Number of native woody recruits") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

recruitnativetreerichness

#ggsave(recruitnativetreerichness, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/recruitnativetreerichness.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
recruitrichnessbysitetypemodel <- glmmTMB(norecruitsperha ~ sitetype + (1|pair), data = beltdatasummarynative, family = poisson(link = "log"))
summary(recruitrichnessbysitetypemodel)
r2(recruitrichnessbysitetypemodel)
1/exp(-0.88515)

# NATIVE RECRUITS BOXPLOT
nativerecruitsBOX <- ggplot(
  data = beltdatasummarynative, aes(x=sitetype, y=norecruitsperha, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  expand_limits(y = 0) +
  labs(x = 'Site', y = "Native woody recruits (stems/ha)") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
nativerecruitsBOX
nativerecruitsBOX <- nativerecruitsBOX + theme(legend.position = "none")



#### Q1/3a. EXOTIC TREE AND SHRUB Richness----

# fill in missing cases for exotic taxa
# beltdatasummary <- beltdatasummary[complete.cases(beltdatasummary),]
beltdatasummary <-  complete(beltdatasummary, origin, fill = list(richness= 0, nostems = 0, norecruits = 0))

# filter to only consider exotic
beltdatasummaryexotic <- filter(beltdatasummary, origin == "Exotic")
describeBy(beltdatasummaryexotic, beltdatasummaryexotic$sitetype)

# old figure of raw data
exotictreerichness <- ggplot(data = beltdatasummaryexotic, aes(x=reorder(site,-richness), y=richness, fill = sitetype)) +
  geom_bar(stat="identity", color ="black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Richness of exotic trees and shrubs") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

exotictreerichness

#ggsave(exotictreerichness, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/exotictreerichness.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
exoticrichnessbysitetypemodel <- glmmTMB(richness ~ sitetype +(1|pair), data = beltdatasummaryexotic, family = poisson(link = "log"))
summary(exoticrichnessbysitetypemodel)
r2(exoticrichnessbysitetypemodel)


# EXOTIC RICHNESS BOXPLOT
exotictreerichnessBOX <- ggplot(
  data = beltdatasummaryexotic, aes(x=sitetype, y=richness, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  expand_limits(y = 0) +
  labs(x = 'Site', y = "Exotic woody richness") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
exotictreerichnessBOX
exotictreerichnessBOX <- exotictreerichnessBOX + theme(legend.position = "none")

#### Q1/3a. EXOTIC TREE AND SHRUB ABUNDANCE ----

# incorporate belt area into exotic count data
beltdatasummaryexotic$totalbeltarea <- beltlengthdata$totalbeltarea[match(beltdatasummaryexotic$site, beltlengthdata$site)]
beltdatasummaryexotic$stemsperha <- round(beltdatasummaryexotic$nostems/beltdatasummarynative$totalbeltarea*10000)
describeBy(beltdatasummaryexotic, beltdatasummaryexotic$sitetype)


# graph of raw data
exotictreeabundance <- ggplot(
  data =  beltdatasummaryexotic, aes(x=reorder(site,-nostems), y=nostems, fill = sitetype)) +
  geom_bar(stat="identity", color="black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Abundance of exotic trees and shrubs") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

exotictreeabundance


# ggsave(exotictreeabundance, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/exotictreeabundance.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model it
exoticabundancebysitetypemodel <- glmmTMB(stemsperha ~ sitetype + (1|pair), data = beltdatasummaryexotic, family = poisson(link = "log"))
summary(exoticabundancebysitetypemodel)
r2(exoticabundancebysitetypemodel)
exp(0.7529)

# EXOTIC Woody plant abundance BOXPLOT
exotictreeabundanceBOX <- ggplot(
  data = beltdatasummaryexotic, aes(x=sitetype, y=stemsperha, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  expand_limits(y = 0) +
  labs(x = 'Site', y = "Exotic woody density") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
exotictreeabundanceBOX
legend <- get_legend(exotictreeabundanceBOX)
exotictreeabundanceBOX <- exotictreeabundanceBOX + theme(legend.position = "none")


# combine belt transect data boxplots
boxmedley <- ggarrange(nativetreerichnessBOX, exotictreerichnessBOX, legend,  nativetreeabundanceBOX, exotictreeabundanceBOX, nativerecruitsBOX, align = "v", ncol = 3, nrow = 2, labels = c("A", "B", "", "C", "D", "E"))
boxmedley <- ggarrange(nativetreerichnessBOX, nativetreeabundanceBOX, nativerecruitsBOX, exotictreerichnessBOX, exotictreeabundanceBOX, legend, align = "v", ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "")) # different order
boxmedley

# ggsave(boxmedley, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/boxplotmedley.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)
# ggsave(boxmedley, filename = "boxplotmedley.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)






#### Q1/2b. 1m QUADRAT - Native Ground Cover species & other ground cover ----
# bring in site type into quadratdata dataframe
quadratdata$sitetype <- sitedata$'Site type'[match(quadratdata$'Site Name', sitedata$'Site Name')]

# bring in site pair into beltdata dataframe
quadratdata$pair <- sitedata$'Paired site code'[match(quadratdata$'Site Name', sitedata$'Site Name')]


# Rename the different coverage types
quadratdata <- quadratdata %>%
  rename(site = 'Site Name',
         transectno = 'Transect Number',
         beltpoint = 'Belt point',
         exotic = 'Exotic',
         native = 'Native plants',
         coarsewd = 'Coarse WD',
         finewd = 'Fine WD & Litter',
         bareground = 'Bare Ground',
         rockcover = 'Rock Cover',
         water = 'Water') 


# summarise data
describe(quadratdata)
describeBy(quadratdata, quadratdata$sitetype)


# identify factors
quadratdata$site <- as.factor(quadratdata$site)
quadratdata$beltpoint <- as.factor(quadratdata$beltpoint)
quadratdata$sitetype <- as.factor(quadratdata$sitetype)


# identify variables
quadratdata$transectno <- as.numeric(quadratdata$transectno)
quadratdata$exotic <-as.numeric(quadratdata$exotic)
quadratdata$native <-as.numeric(quadratdata$native)
quadratdata$coarsewd <-as.numeric(quadratdata$coarsewd)
quadratdata$finewd <-as.numeric(quadratdata$finewd)
quadratdata$bareground <-as.numeric(quadratdata$bareground)
quadratdata$rockcover <-as.numeric(quadratdata$rockcover)
quadratdata$water <-as.numeric(quadratdata$water)


# summarise
quadratdatasummary <- quadratdata %>%
  group_by(sitetype, pair, site) %>%
  summarise(
    exotic = mean(exotic, na.rm = TRUE), 
    native = mean(native, na.rm = TRUE),
    coarsewd = mean(coarsewd,na.rm = TRUE ),
    finewd = mean(finewd,na.rm = TRUE ),
    bareground = mean(bareground, na.rm = TRUE ),
    rockcover = mean(rockcover, na.rm = TRUE ),
    water = mean(water, na.rm = TRUE)
  )
quadratdatasummary


# summarise data
describe(quadratdatasummary)
describeBy(quadratdatasummary, quadratdatasummary$sitetype)



# graph native ground cover percentages
nativegroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-native), y=native, fill = sitetype)) +
  geom_bar(stat="identity", color= "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage native of ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

nativegroundcoverpercentage


#ggsave(nativegroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativegroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model native ground cover using glmm
quadratdata$nativebeta <- (quadratdata$native+1)/102
nativegroundcoverbysitetypeglmm <- glmmTMB(nativebeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(nativegroundcoverbysitetypeglmm)
r2(nativegroundcoverbysitetypeglmm)


# graph exotic ground cover percentages
exoticgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-exotic), y=exotic, fill = sitetype)) +
  geom_bar(stat="identity",colour = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Exotic groundcover (%)") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  theme(axis.ticks.x = element_blank())

exoticgroundcoverpercentage

# ggsave(exoticgroundcoverpercentage, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/exoticgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 300)


# model exotic ground cover using glmm
quadratdata$exoticbeta <- (quadratdata$exotic+1)/102
exoticgroundcoverbysitetypeglmm <- glmmTMB(exoticbeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(exoticgroundcoverbysitetypeglmm)
r2(exoticgroundcoverbysitetypeglmm)


# graph finewd ground cover percentages
finewdgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-finewd), y=finewd, fill = sitetype)) +
  geom_bar(stat="identity",colour = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage of fine woody debris ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
finewdgroundcoverpercentage

#ggsave(finewdgroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/finewdgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model finewd ground cover using glmm
quadratdata$finewdbeta <- (quadratdata$finewd+1)/102
finewdgroundcoverbysitetypeglmm <- glmmTMB(finewdbeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(finewdgroundcoverbysitetypeglmm)


# graph coarsewd ground cover percentages
coarsewdgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-coarsewd), y=coarsewd, fill = sitetype)) +
  geom_bar(stat="identity",color = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage of coarse woody debris ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

coarsewdgroundcoverpercentage

# ggsave(coarsewdgroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/coarsewdgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model coarsewd ground cover using glmm 
quadratdata$coarsewdbeta <- (quadratdata$coarsewd+1)/102
coarsewdgroundcoverbysitetypeglmm <- glmmTMB(coarsewdbeta ~ sitetype + (1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(coarsewdgroundcoverbysitetypeglmm)


# graph  bareground ground cover percentages
baregroundgroundcoverpercentage <- ggplot(
  data = quadratdatasummary,  aes(x=reorder(site,-bareground), y=bareground, fill = sitetype)) +
  geom_bar(stat="identity",color = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Percentage of bare ground cover") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
baregroundgroundcoverpercentage

# ggsave(baregroundgroundcoverpercentage, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/baregroundgroundcoverpercentage.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model  bareground ground cover using glmm 
quadratdata$baregroundbeta <- (quadratdata$ bareground+1)/102
baregroundgroundcoverbysitetypeglmm <- glmmTMB(baregroundbeta ~ sitetype+(1|site/pair), data = quadratdata, family=beta_family(link="logit"))
summary(baregroundgroundcoverbysitetypeglmm)



#### Q1/2b 1m QUADRAT - COMPAIRING AVERAGE COVERAGE TYPE ----

#COMBINED GRAPH
# pivot longer
quadratdatasummarylong <- quadratdatasummary %>%
  pivot_longer(cols = exotic:bareground, names_to = "covertype", values_to = "cover")
quadratdatasummarylong$covertype <- factor(quadratdatasummarylong$covertype, levels = c("native", "exotic", "finewd", "coarsewd", "bareground"))


#summarise
quadratdatalongsummary <- quadratdatasummarylong %>%
  group_by(sitetype, covertype) %>%
  summarise(
    meancover = mean(cover),
    coverse = sd(cover)/sqrt(n()))



# graph it
covertypecomparisonfigure <- ggplot(
  data = quadratdatalongsummary, aes(x=covertype, y= meancover, fill = sitetype)) +
  geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin=meancover-coverse, ymax=meancover+coverse), position = position_dodge(0.9, preserve = "single"), width=0.1, size=0.5, color="black") +
  ylim(0, 100) +
  labs(x= NULL, y = "Understorey cover %") +
  scale_x_discrete(breaks=c("native","exotic","finewd", "coarsewd", "bareground"),
                   labels=c("Native plants","Exotic plants","Leaf litter", "Woody debris", "Bare ground")) +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site type", labels = c("Remnant", "Revegetated")))+
  guides(fill=guide_legend(title="Site type", labels = c("Remnant", "Revegetated")))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.ticks.x = element_blank())

covertypecomparisonfigure

# ggsave(covertypecomparisonfigure, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/understoreycovers.tiff", width = 210, height = 120, units = c("mm"), bg = "white", dpi = 300)




#### Q1/2c STRATA HEIGHTS - CANOPY COVER VEGETATION HEIGHT DIVERSITY ----

#  CANOPY COVER 
# lets look at the canopy first
# bring in site type into canopydata dataframe
canopydata$sitetype <- sitedata$'Site type'[match(canopydata$'Site Code', sitedata$'Site code')]

# bring in site pair into beltdata dataframe
canopydata$pair <- sitedata$'Paired site code'[match(canopydata$'Site Name', sitedata$'Site Name')]

# rename 'messy' factor names
canopydata <- canopydata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         beltm = `Belt intercept`,
         canopy = 'Canopy')


# identify factors
canopydata$sitetype <- as.factor(canopydata$sitetype)
canopydata$pair <- as.factor(canopydata$pair)
canopydata$site <- as.factor(canopydata$site)
canopydata$transect <- as.factor(canopydata$transect)
canopydata$beltm <- as.factor(canopydata$beltm)
canopydata$canopy <- as.factor(canopydata$canopy)
canopydata$canopy <- factor(canopydata$canopy,levels=c("P","A"),labels=c('presence','absence'))

# exclude transects rows with NAs
canopydata <- canopydata[complete.cases(canopydata), ]


# create frequency table of canopy
dffreq <- canopydata %>%
  count(sitetype, pair, site, canopy) %>%
  group_by(site) %>%          
  mutate(prop = prop.table(n))
dffreq <- as.data.frame(dffreq)


# consider only presence
dffreqpres <- filter(dffreq, canopy == "presence")


# freq figure for presence of canopy by siteypte NB. Definitely probably worth separating out East and West sites.
frequencyofcanopy <- ggplot(
  data = dffreqpres, aes(x=reorder(site,-prop), y= prop*100, fill = sitetype)) +   geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  ylim(0, 100) +
  labs(x="Site", y = "Canopy cover %")+
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())

frequencyofcanopy

# ggsave(frequencyofcanopy, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/frequencyofcanopy.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# model it
canopybysitetypeglmm <- glmmTMB(prop ~ sitetype + (1|pair), data = dffreqpres, family=beta_family(link="logit"))
summary(canopybysitetypeglmm)

describe(dffreqpres)
describeBy(dffreqpres, dffreqpres$sitetype)


# Canopy presence BOXPLOT
canopyBOX <- ggplot(
  data = dffreqpres, aes(x=sitetype, y=prop, fill = sitetype)) +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_jitter(width = 0.1) +
  ylim(0,1) +
  labs(x = 'Site type', y = "Canopy presence %") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())
canopyBOX
# legend <- get_legend(canopyBOX)
# canopyBOX <- canopyBOX + theme(legend.position = "none")

# combine with boxmedley to create Figure 3 for paper
Fig3 <- ggarrange(nativetreerichnessBOX, nativetreeabundanceBOX,nativerecruitsBOX, exotictreerichnessBOX, exotictreeabundanceBOX, canopyBOX, common.legend = TRUE, legend = "bottom", align = "v", ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F")) # different order
Fig3

# ggsave(Fig3, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/Fig3.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)


### ASSESS VEGETATION STRUCTURE - TREE HEIGHT AND SIZE ----

# bring in site type
heightdata$sitetype <- sitedata$'Site type'[match(heightdata$'Site Code', sitedata$'Site code')]

# bring in site pair into beltdata dataframe
heightdata$pair <- sitedata$'Paired site code'[match(heightdata$'Site Name', sitedata$'Site Name')]


# rename 'messy' factor names
heightdata <- heightdata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         species = `Species name`,
         height = 'Heights (first 5 of each Sp to 0.5m)')

head(heightdata)
species <- levels(as.factor(heightdata$species))
species


# identify factors
heightdata$sitetype <- as.factor(heightdata$sitetype)
heightdata$pair <- as.factor(heightdata$pair)
heightdata$site <- as.factor(heightdata$site)
heightdata$transect <- as.factor(heightdata$transect)
heightdata$species <- as.factor(heightdata$species)


# identify variables
heightdata$height <- as.numeric(heightdata$height)

head(heightdata)


# graph raw heights using boxplots?
heightfig <- ggplot(
  data = heightdata, aes(x=species, y= height, fill = sitetype)) +
  geom_boxplot() +
  labs(x= NULL, y = "Height (m)") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfig

# there are LOTS of different species... subset?

# graph raw heights using boxplots at sitetype level only?
heightfig2 <- ggplot(
  data = heightdata, aes(x=sitetype, y= height, fill = sitetype)) +
  geom_boxplot() +
  labs(x= NULL, y = "Height (m)") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfig2


# summarise heightdata
heightdatasummary <- heightdata %>%
  group_by(sitetype, pair, site, species) %>%
  summarise(
    count = n(),
    meanheight = mean(height))

describeBy(heightdatasummary, heightdatasummary$sitetype)
# new <- heightdatasummary[which(heightdatasummary$count >= 3),]


# filter out species for which <3 plants were surveyed at a site?
heightdatasummarytrunc <- filter(heightdatasummary, count > 2)

# at how many sites were those species found? 
speciesbysite <- heightdatasummarytrunc %>%
  group_by(sitetype, species) %>%
  summarise(
    count = n(),
    height = mean(meanheight))


# Matt's fancy code to select only species which occur at at least three site of a given sitetype
new <- speciesbysite[which(speciesbysite$count > 2),]

want <- unique(as.character(new$species))

test <- NULL
for(i in 1:NROW(want)) {
  
  get <- heightdata[which(as.character(heightdata$species) == want[i]),]
  test <- rbind(test, get) }

test <- filter(test, species != 'Eucalyptus ovata')


# graph raw heights using boxplots using truncated data?
heightfigtrunc <- ggplot(
  data = test, aes(x=species, y= height, fill = sitetype)) +
  geom_boxplot(position = position_dodge(preserve = 'single')) +
  labs(x= "Species", y = "Plant height (m)") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(fill=guide_legend(title="Site type"))+
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfigtrunc

# ggsave(heightfigtrunc, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/Fig S1.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)




# graph raw heights using boxplots at sitetype level only?
heightfigtrunc2 <- ggplot(
  data = test, aes(x=sitetype, y= height, fill = sitetype)) +
  geom_boxplot() +
  labs(x= NULL, y = "Height (m)") +
  theme_classic() +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))
heightfigtrunc2


# Model of plant heights by species for truncated data
heightslmm <- glmmTMB(height ~ sitetype*species + (1|site/pair), data = test, family=gaussian())
summary(heightslmm)
plot_model(heightslmm)




### STRATA DATA - ASSESS VEGETATION STRUCTURE - VEGETATION STRUCTURE ----
# load data
sitedata <- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Site info")
stratadata<- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Strata Heights")

# lets look at the canopy first
# bring in site type into strata dataframe
stratadata$sitetype <- sitedata$'Site type'[match(stratadata$'Site Code', sitedata$'Site code')]

# bring in site pair into strata dataframe
stratadata$pair <- sitedata$'Paired site code'[match(stratadata$'Site Name', sitedata$'Site Name')]


# rename 'messy' factor names
stratadata <- stratadata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         beltm = `Belt intercept`,
         height = 'Poll height')


# identify factors
stratadata$sitetype <- as.factor(stratadata$sitetype)
stratadata$pair <- as.factor(stratadata$pair)
stratadata$site <- as.factor(stratadata$site)
stratadata$transect <- as.factor(stratadata$transect)
stratadata$beltm <- as.factor(stratadata$beltm)
stratadata$height <- as.factor(stratadata$height)
canopydata$canopy <- factor(canopydata$canopy,levels=c("P","A"),labels=c('presence','absemce'))

# identify variables
stratadata$Tree <- as.numeric(stratadata$Tree)
stratadata$Shrub <- as.numeric(stratadata$Shrub)
stratadata$Fern <- as.numeric(stratadata$Fern)
stratadata$Herbaceous <- as.numeric(stratadata$Herbaceous)

# exclude transects rows with NAs
stratadata <- stratadata[complete.cases(stratadata), ] # OR
stratadata <- na.omit(stratadata)

# create count able of different plant types
stratacount <- stratadata %>%
  group_by(sitetype, pair, site, height) %>%          
  summarise(
    count = n(),
    tree = sum(Tree),
    shrub = sum(Shrub),
    fern = sum(Fern),
    herb = sum(Herbaceous))

stratacountworks <- filter(stratacount, sitetype == "Works")


describe(stratacountworks)

# create frequency table of different plant types
stratasum <- stratadata %>%
  group_by(sitetype, height) %>%          
  summarise(
    count = n(),
    tree = sum(Tree)/n(),
    shrub = sum(Shrub)/n(),
    fern = sum(Fern)/n(),
    herb = sum(Herbaceous)/n())


# convert to long form for graphing
stratasumlong <- stratasum %>%
  pivot_longer(cols = tree:herb, names_to = "planttype", values_to = "frequency")


# Figure for frequency of different plant types at different heights by sitetype


# OPTION 1
stratavegtypecomparison1 <- ggplot(
  data = stratasumlong, aes(x=height, y= frequency, fill = sitetype)) +
  geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  ylim(0, 1) +
  labs(x= NULL, y = "Frequency") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~planttype)

stratavegtypecomparison1

# ggsave(stratavegtypecomparison1, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/stratavegtypecomparison1.tiff", width = 16, height = 12, units = "cm", dpi = 600)



# OPTION 2 #### THE ONE I"M USING # NB. NO ERROR BARS! Need to create df from site summary (which is below) to create errorbars
stratasitesypecomparison <- ggplot(
  data = stratasumlong, aes(x=height, y= frequency, fill = planttype)) +
  geom_bar(stat = "identity", aes(fill = planttype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  ylim(0, 1) +
  labs(x= NULL, y = "Frequency") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Plant form"))+
  guides(fill=guide_legend(title="Plant form"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~sitetype, ncol = 1)

stratasitesypecomparison

#ggsave(stratasitesypecomparison, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/stratasitesypecomparison.tiff", width = 16, height = 12, units = "cm", dpi = 600)

stratacomparesitetypefig <- ggplot(data =stratasumlong, aes(y=height, x=planttype, fill=frequency, size=frequency, colour = frequency)) +
  geom_point(alpha=0.5, shape=21) +
  scale_size(range = c(1, 30), name="Frequency of hits") +
  facet_wrap(.~sitetype) +
  theme_classic()+
  labs(x= "Plant type", y="Pole Heights")
stratacomparesitetypefig


stratacomparesitetypefigurecircle <- ggplot(data =stratasumlong, aes(y=height, x=planttype, fill=frequency, size=frequency)) +
  geom_point(alpha=1, shape=21, color="black") +
  facet_wrap(.~sitetype) +
  labs(x= "Plant type", y="Pole Heights") +
  scale_fill_continuous(name = "Frequency",  limits=c(0, 1), breaks=seq(0, 1, by=0.2)) + 
  scale_size_continuous(name = "Frequency", range = c(1, 15), limits=c(0, 1), breaks=seq(0, 1, by=0.2)) +
  scale_fill_steps2(
    low = "darkorange", 
    mid = "purple", 
    high = "darkgreen", 
    midpoint = 0.6)+
  guides(fill=guide_colourbar(reverse=TRUE)) + 
  guides(fill=guide_legend(), size = guide_legend()) +
  theme_classic()

stratacomparesitetypefigurecircle

#ggsave(stratacomparesitetypefigurecircle, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/stratacomparesitetypefigurecircle.tiff", width = 16, height = 12, units = "cm", dpi = 600)

#summarise 
stratasumbysite <- stratadata %>%
  group_by(sitetype, pair, site, height) %>%    #  NEED TO ADD SITE AS GROUPING FACTOR HERE     
  summarise(
    count = n(),
    tree = sum(Tree)/n(),
    shrub = sum(Shrub)/n(),
    fern = sum(Fern)/n(),
    herb = sum(Herbaceous)/n())


# convert to long form for graphing
stratasumbysitelong <- stratasumbysite %>%
  pivot_longer(cols = tree:herb, names_to = "planttype", values_to = "frequency")


#summarise by sitetype and planttype
stratasumbysitelongsum <- stratasumbysitelong %>%
  group_by(sitetype, height, planttype) %>%    #  NEED TO ADD SITE AS GROUPING FACTOR HERE     
  summarise(
    meanfreq = mean(frequency),
    sefreq = sd(frequency)/sqrt(n()))


# FIGURE FOR PAPER
stratasumbysitelongsum$planttype <- factor(stratasumbysitelongsum$planttype, levels=c("tree","shrub","herb","fern"),labels=c("Tree","Shrub","Herb","Fern"))

stratavegtypesitetypecomparison <- ggplot(
  data = stratasumbysitelongsum, aes(x=height, y= meanfreq, fill = sitetype)) +
  geom_bar(stat = "identity", aes(fill = sitetype), colour = "black", position = position_dodge2(width = 0.9, preserve = "single")) +
  geom_errorbar(aes(ymin=meanfreq-sefreq, ymax=meanfreq+sefreq), position = position_dodge(0.9, preserve = "single"), width=0.1, size=0.5, color="black") +
  ylim(0, 1) +
  labs(x= "Vegetation height (m)", y = "Frequency of hits") +
  scale_fill_brewer(palette="Set2", labels = c("Remnant", "Revegetated"))+
  guides(colour=guide_legend(title="Site type"))+
  guides(fill=guide_legend(title="Site type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  facet_wrap(~planttype) +
  theme(strip.text.x = element_text(size = 14))
stratavegtypesitetypecomparison 
stratavegtypesitetypecomparison <- stratavegtypesitetypecomparison +coord_flip()
stratavegtypesitetypecomparison

# ggsave(stratavegtypesitetypecomparison, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/stratavegtypesitetypecomparison.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)


## STRATA MODELS
# if you need summary stats adapt this
stratasumbysitelongsumtemp <- filter(stratasumbysitelongsum, planttype == "Fern")
describeBy(stratasumbysitelongsumtemp, stratasumbysitelongsumtemp$sitetype)


# model it - tree frequency by height and site type
stratasumbysite$treebeta <- (stratasumbysite$tree+.01)
stratatreehitfrequencyglmm <- glmmTMB(treebeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(stratatreehitfrequencyglmm)
Anova(stratatreehitfrequencyglmm)

emmip(stratatreehitfrequencyglmm, height ~ sitetype)

#model it - fern frequency by height and site type
stratasumbysite$fernbeta <- (stratasumbysite$fern+.01)
stratafernhitfrequencyglmm <- glmmTMB(fernbeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(stratafernhitfrequencyglmm)
Anova(stratafernhitfrequencyglmm)

emmip(stratafernhitfrequencyglmm, height ~ sitetype)



stratasumbysite$shrubbeta <- (stratasumbysite$shrub+.01)
stratashrubhitfrequencyglmm <- glmmTMB(shrubbeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(stratashrubhitfrequencyglmm)
Anova(stratashrubhitfrequencyglmm)

emmip(stratashrubhitfrequencyglmm, height ~ sitetype)

stratasumbysite$herbbeta <- (stratasumbysite$herb+.01)/1.02
strataherbhitfrequencyglmm <- glmmTMB(herbbeta ~ sitetype * height + (1|pair/site), data = stratasumbysite, family=beta_family(link="logit"))
summary(strataherbhitfrequencyglmm)
Anova(strataherbhitfrequencyglmm)

emmip(strataherbhitfrequencyglmm, height ~ sitetype)






#### Q3 LANDSCAPE CONTEXT - NATIVE WOODY VEGETATION COVERAGE - RECRUIT ABUNDANCE ---- 

# consider only works sites
beltdatasummarynativeworks <- filter(beltdatasummarynative, sitetype == "Works")

# graph native and exotic species richness
nativewoodyvegcoveragerecruits <- ggplot(data = beltdatasummarynativeworks, aes(x=nativevegcoveragepercentage*100, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  geom_smooth (method = "glm", method.args = list(family = "poisson"),
               colour = "darkblue", se = TRUE) +
  labs(x = 'Remnant vegetation %', y = "Native woody recruits (stems/ha)") +
  annotate("text", x = 60, y = 3000, label = "italic(R) ^ 2 == 0.01",
           parse = TRUE) +
  theme_classic()
nativewoodyvegcoveragerecruits

#ggsave(nativewoodyvegcoveragerecruits, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativewoodyvegcoveragerecruits.tiff", width = 16, height = 12, units = "cm", dpi = 600)

nativevegcovervrecruits <- glm(norecruitsperha ~ nativevegcoveragepercentage, data = beltdatasummarynativeworks,   family = poisson(link = "log"))
summary(nativevegcovervrecruits)
r2(nativevegcovervrecruits)
with(summary(nativevegcovervrecruits), 1 - deviance/null.deviance)
model_performance(nativevegcovervrecruits)


# #COMPARE TO EXOTIC COVER
# 
# groundcovervsnativecoverage <- full_join(quadratdatasummary, beltdatasummary, by='site')
# groundcovervsnativecoverageworks <- filter(groundcovervsnativecoverage, sitetype.x == "Works")
# describe(groundcovervsnativecoverageworks)
# 
# #graph it 
# exoticcovervsnativecoveragegraph <- ggplot(data = groundcovervsnativecoverageworks, aes(x=nativevegcoveragepercentage, y=exotic)) +
#   geom_jitter(size = 2, color = "darkorange3") +
#   
#   labs(x = 'Percentage of Native Woody Vegetation', y = "Exotic ground cover (%)") +
#   theme_classic()
# exoticcovervsnativecoveragegraph

# #model it - need to do beta
# groundcovervsnativecoverage$exoticbeta <- groundcovervsnativecoverage$exotic/100
# 
# exoticgroundcovervsnativecoverageglm <- glm(exoticbeta ~ nativevegcoveragepercentage, data = groundcovervsnativecoverage,   family = beta_family(link="logit"))
# summary(exoticgroundcovervsnativecoverageglm)
# summ(exoticgroundcovervsnativecoverageglm)


# #COMPARE TO NATIVE COVER
# #graph it 
# nativecovervsnativecoveragegraph <- ggplot(data = groundcovervsnativecoverageworks, aes(x=nativevegcoveragepercentage, y=native)) +
#   geom_jitter(size = 2, color = "darkorange3") +
#   
#   labs(x = 'Percentage of Native Woody Vegetation', y = "Native ground cover(%)") +
#   theme_classic()
# nativecovervsnativecoveragegraph
# 
# #ggsave(nativecovervsnativecoveragegraph, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nativecovervsnativecoveragegraph.tiff", width = 16, height = 12, units = "cm", dpi = 600)
# 
# 
# #model it - need to do beta
# groundcovervsnativecoverage$nativebeta <- groundcovervsnativecoverage$native/100
# 
# nativecovervsnativecoverageglm <- glm(nativebeta ~ nativevegcoveragepercentage, data = groundcovervsnativecoverage,   family = beta_family())
# summary(nativecovervsnativecoverageglm)
# summ(nativecovervsnativecoverageglm)



#### Q3. BROWSER POO FREQUENCY ----

# load data
poodata <- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Pellet counts")

# bring in site type into poodata dataframe
poodata$sitetype <- sitedata$'Site type'[match(poodata$'Site Code', sitedata$'Site code')]

# bring in site pair into beltdata dataframe
poodata$pair <- sitedata$'Paired site code'[match(poodata$'Site Code', sitedata$'Site code')]

# rename 'messy' factor names
poodata <- poodata %>%
  rename(site = 'Site Name',
         transect = 'Transect Number',
         belt = `Belt point`)


# identify factors
poodata$sitetype <- as.factor(poodata$sitetype)
poodata$pair <- as.factor(poodata$pair)
poodata$site <- as.factor(poodata$site)
poodata$transect <- as.factor(poodata$transect)
poodata$belt <- as.factor(poodata$belt)


# replace P's and A's with 1s and 0s
poodata <- poodata[!(poodata$Deer == "NA"),] # remove one row with NAs
poodata[poodata  == "P"] <- "1"
poodata[poodata  == "A"] <- "0"
poodata$Kangaroo <- as.numeric(poodata$Kangaroo)
poodata$Wallaby <- as.numeric(poodata$Wallaby)
poodata$Rabbit <- as.numeric(poodata$Rabbit)
poodata$Deer <- as.numeric(poodata$Deer)

# create count able of browser poo counts
poodatacount <- poodata %>%
  group_by(sitetype, pair, site) %>%          
  summarise(
    count = n(),
    kanga = sum(Kangaroo),
    wall = sum(Wallaby),
    rabbit = sum(Rabbit),
    deer = sum(Deer))

# create table of browser poo frequency # NOTE DF NAME SAME AS ABOVE 
poodatacount <- poodata %>%
  group_by(sitetype, pair, site) %>%          
  summarise(
    count = n(),
    kanga = sum(Kangaroo)/count,
    wall = sum(Wallaby)/count,
    rabbit = sum(Rabbit)/count,
    deer = sum(Deer)/count)


# sum across all browser types
poodatacount$total <- poodatacount$kanga + poodatacount$wall + poodatacount$rabbit + poodatacount$deer 
describe(poodatacount)

# graph browser frequency 
browsertotal <- ggplot(
  data = poodatacount, aes(x=reorder(site,-total), y=total, fill = sitetype)) +
  geom_bar(stat="identity", colour = "black", position = position_dodge(preserve = "single")) +
  labs(x = 'Site', y = "Frequency of browser pellets") +
  scale_fill_brewer(palette="Set2")+
  guides(colour=guide_legend(title="Site Type"))+
  guides(fill=guide_legend(title="Site Type"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  theme(axis.ticks.x = element_blank())
browsertotal

#ggsave(browsertotal, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/browsertotale.tiff", width = 16, height = 12, units = "cm", dpi = 600)

# merge poodata and beltdata
beltdatasummarynative$browsers <- poodatacount$total[match(beltdatasummarynative$site, poodatacount$site)]  

beltdatasummarynativeworks <- filter(beltdatasummarynative, sitetype == "Works")

browservrecruit <- ggplot(data = beltdatasummarynativeworks, aes(x=browsers, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  geom_smooth (method = "glm",method.args = list(family = "poisson"),
               colour = "darkblue") +
  labs(x = 'Browser pellet frequency ', y = "Native woody recruits stems/ha)") +
  annotate("text", x = 0.2, y = 3000, label = "italic(R) ^ 2 == 0.25",
           parse = TRUE) +
  theme_classic()
browservrecruit

#ggsave(browservrecruit, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/browservrecruit.tiff", width = 16, height = 12, units = "cm", dpi = 600)


browservrecruitmodel <- glm(norecruitsperha ~ browsers, data = beltdatasummarynativeworks, family = poisson(link = "log"))
summary(browservrecruitmodel)
r2(browservrecruitmodel)
with(summary(browservrecruitmodel), 1 - deviance/null.deviance)
model_performance(browservrecruitmodel)



#### Q3. SOIL ANALYSIS bulk density ----

beltdatasummarynative <- filter(beltdatasummarynative, sitetype == "Works")
describeBy(beltdatasummarynative, beltdatasummarynative$sitetype)

#Graph bulk density to recruits 
bulkdensityvsrecruits <- ggplot(data = beltdatasummarynative, aes(x=bulkdensityaverage, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  labs(x = 'Bulk density (g/cm3)', y = "Native woody recruits (stems/ha)")+
  # annotate("text", x = 1, y = 3000, label = "italic(R) ^ 2 == 0.28",
  #         parse = TRUE) +
  theme_classic()
bulkdensityvsrecruits

#ggsave(bulkdensityvsrecruits, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/bulkdensityvsrecruitss.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model it 
bulkdensityvsrecruitsglm <- glm(norecruitsperha ~ bulkdensityaverage, data = beltdatasummarynative,   family = poisson(link = "log"))
summary(bulkdensityvsrecruitsglm)
r2(bulkdensityvsrecruitsglm)
with(summary(bulkdensityvsrecruitsglm), 1 - deviance/null.deviance)
model_performance(bulkdensityvsrecruitsglm)


#### Q3. SOIL ANALYSIS nutrients ----

# Nitrogen affect on recruits - Graph it
nitrogenvsrecruits <- ggplot(data = beltdatasummarynative, aes(x=totalnitrogen, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  geom_smooth (method = "glm",method.args = list(family = "poisson"), colour = "darkblue") +
  annotate("text", x = 6000, y = 3000, label = "italic(R) ^ 2 == 0.14",
           parse = TRUE) +
  labs(x = 'Total Nitrogen mg/kg ', y = "Native woody recruits (stems/ha)") +
  theme_classic()
nitrogenvsrecruits

#ggsave(nitrogenvsrecruits, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nitrogenvsrecruits.tiff", width = 16, height = 12, units = "cm", dpi = 600)


# model it
nitrogenvsrecruitsglm <- glm(norecruitsperha ~ totalnitrogen, data = beltdatasummarynativeworks,   family = poisson(link = "log"))
summary(nitrogenvsrecruitsglm)
r2(nitrogenvsrecruitsglm)
with(summary(nitrogenvsrecruitsglm), 1 - deviance/null.deviance)
model_performance(nitrogenvsrecruitsglm)


# Phosphorus affect on recruits - Graph it 
phosphorusvsrecruits <- ggplot(data = beltdatasummarynative, aes(x=totalphosporus, y=norecruitsperha)) +
  geom_jitter(size = 2, color = "#FC8D62") +
  geom_smooth (method = "glm",method.args = list(family = "poisson"), colour = "darkblue") +
  annotate("text", x = 1000, y = 3000, label = "italic(R) ^ 2 == 0.02",
           parse = TRUE) +
  labs(x = 'Total Phosphorus mg/kg ', y = "Native woody recruits (stems/ha)") +
  theme_classic()
phosphorusvsrecruits


#ggsave(phosphorusvsrecruits, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/phosphorusvsrecruits.tiff", width = 16, height = 12, units = "cm", dpi = 600)

#model it
phosphorusvsrecruitsglm <- glm(norecruitsperha ~ totalphosporus, data = beltdatasummarynative,   family = poisson(link = "log"))
summary(phosphorusvsrecruitsglm)
r2(phosphorusvsrecruitsglm)
with(summary(phosphorusvsrecruitsglm), 1 - deviance/null.deviance)
model_performance(phosphorusvsrecruitsglm)


# MEDLY PANELS OF DRIVERS OF RECRUITMENT
linearpanels <- ggarrange(nativewoodyvegcoveragerecruits, browservrecruit, bulkdensityvsrecruits, nitrogenvsrecruits, phosphorusvsrecruits, align = "hv", ncol = 3, nrow = 2, labels = c("A", "B", "C", "D", "E", "F"), hjust = -17) # different order
linearpanels

# ggsave(linearpanels, filename = "~/uomShare/wergProj/Eliza_Thesis_Nov22/figures/linearpanels.tiff", width = 210, height = 140, units = c("mm"), bg = "white", dpi = 300)


# #### Soil analysis - ground coverage data  ----
# soildata <- read_xlsx("Riparian works monitoring data_20231013.xlsx", sheet = "Soil samples")
# soildata$sitetype <- sitedata$'Site type'[match(soildata$'Site Name', sitedata$'Site Name')]
# 
# 
# soildata <- soildata %>%
#   rename(site = 'Site Name',
#          bulkdensityraverage = 'Bulk density average',
#          totalnitrogen = 'total nitrogen mg/kg',
#          totalphosporus = 'total phosporus mg/kg') 
# 
# # identify factors
# soildata$site <- as.factor(soildata$site)
# 
# # identify variables
# soildata$bulkdensityraverage <- as.numeric(soildata$bulkdensityraverage)
# soildata$totalnitrogen <- as.numeric(soildata$totalnitrogen)
# soildata$totalphosporus <- as.numeric(soildata$totalphosporus)
# 
# soildatasummary <- soildata %>%
#   group_by(sitetype, site) %>%
#   summarise(
#     bulkdensityraverage = mean(bulkdensityraverage),
#     totalnitrogen = mean(totalnitrogen),
#     totalphosporus = mean(totalphosporus))
# 
# 
# 
# #merge quadratdata summary and soildata 
# covervssoil <- full_join(quadratdatasummary, soildata, by='site')
# 
# #compare remant vs work sites
# boxplot(totalnitrogen ~ sitetype.x, data =covervssoil)
# boxplot(totalphosporus ~ sitetype.x, data =covervssoil)
# boxplot(bulkdensityraverage ~ sitetype.x, data =covervssoil)
# 
# #restrict analysis to worksites 
# covervssoilworks <- filter(covervssoil, sitetype.x == "Works")
# 
# 
# # Nitrogen affect on recruits - Graph it
# nitrogenvsecoticcover <- ggplot(data = covervssoilworks, aes(x=totalnitrogen, y=exotic)) +
#   geom_jitter(size = 2, color = "darkorange2") +
#   labs(x = 'Total nitrogen (mg/kg) ', y = "Exotic ground cover (%)") +
#   theme_classic()
# nitrogenvsecoticcover
# 
# #ggsave(nitrogenvsecoticcover, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/nitrogenvsecoticcovers.tiff", width = 16, height = 12, units = "cm", dpi = 600)
# 
# 
# # model it - need to do beta
# covervssoilworks$exoticbeta <- covervssoilworks$exotic/100
# nitrogenvsecoticcoverglm <- glm(exoticbeta ~ totalnitrogen, data = covervssoilworks,   family = beta_family())
# summary(nitrogenvsecoticcoverglm)
# 
# 
# # Phosphorus affect on recruits 
# phosphorusvsexoticcover <- ggplot(data = covervssoilworks, aes(x=totalphosporus, y=exotic)) +
#   geom_jitter(size = 2, color = "darkorange2") +
#   labs(x = 'Total phosphorus (mg/kg) ', y = "Exotic ground cover (%)") +
#   theme_classic()
# 
# phosphorusvsexoticcover
# 
# #ggsave(phosphorusvsexoticcover, filename = "C:/Users/Eliza.Foley-Congdon/OneDrive - Water Technology Pty Ltd/Desktop/Eliza Uni/My thesis/phosphorusvsexoticcover.tiff", width = 16, height = 12, units = "cm", dpi = 600)
# 
# 
# 
# #model it - need to do beta
# covervssoilworks$exoticbeta <- covervssoilworks$exotic/100
# phosphorusvsecoticcoverglm <- glm(exoticbeta ~ totalphosporus, data = covervssoilworks,   family = beta_family())
# summary(phosphorusvsecoticcoverglm)
# r2(phosphorusvsecoticcoverglm)
# 