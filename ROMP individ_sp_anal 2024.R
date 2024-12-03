
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
beltdata24$lifeform <- floravic$'NVIS_GF'[match(beltdata24$'Sp_name', floravic$'SCI_NAME')]
beltdata24$Month_Yr <- format(as.Date(beltdata24$Date), "%Y-%m")
beltdata24 <- mutate(beltdata24, 
                     origin = ifelse(origin == "Native but some stands may be alien", "Native", origin)) 
beltdata24 <- beltdata24 %>%
  rename(site = 'WptID',
         species = 'Sp_name',
         type = 'Type (reveg, existing, recruit)',
         alive = 'Count_alive',
         dead = 'Count_dead',
         visit = 'Visit_number', 
         belt = 'Belt_ID')

# identify factors
beltdata24$site <- as.factor(beltdata24$site)
beltdata24$species <- as.factor(beltdata24$species)
beltdata24$origin <- as.factor(beltdata24$origin)
beltdata24$sitetype <- as.factor(beltdata24$sitetype)
beltdata24$pair <- as.factor(beltdata24$pair)
describeBy(beltdata24, beltdata24$sitetype)

# Check species without any records
missing_val <- beltdata24 %>% filter(alive==0, dead==0)
#writexl::write_xlsx(missing_val, '~/uomShare/wergProj/W12 - Revegetation/ROMP surveys/2024_25 outputs/missing_val.xlsx')
# Check which columns have missing values
sapply(beltdata24, function(x) sum(is.na(x)))

#### Q1/2a. NATIVE TREE AND SHRUB Richness ----

beltdata24 <- beltdata24 %>% filter(!is.na(site))
beltdata24 <- beltdata24 %>% replace_na(list(dead = 0))
beltdata24 <- beltdata24 %>% replace_na(list(alive = 0))
beltdata24 <- beltdata24 %>% mutate(count= alive + dead)

beltdata24 <- filter(beltdata24, species != 'Unknown')
beltdata24 <- filter(beltdata24, species != 'Nil')

############################ Individual Species ##########################################

# Analysing revegetated sites that only have native species undertaken on the 13 sites visited 3 times
beltdata24_no_weeds <- filter(beltdata24, origin != 'Introduced')
beltdata24_reveg <- filter(beltdata24_no_weeds, type == c('reveg'))
#writexl::write_xlsx(beltdata24_reveg, '~/uomShare/wergProj/W12 - Revegetation/ROMP surveys/2024_25 outputs/beltdata24_reveg.xlsx')
beltdata24_reveg_2 <- filter(beltdata24_reveg, visit == '2')
beltdata24_reveg_3 <- filter(beltdata24_reveg, visit == '3')
Indiv_data_3 <- beltdata24_reveg_3 %>%
  group_by(site, catchment, species, lifeform) %>%
  summarise_at(vars(count, alive), tibble::lst(mean, median))
Indiv_data_3a <- Indiv_data_3 %>% group_by(lifeform) %>% count (species)
Indiv_data_2 <- beltdata24_reveg_2 %>%
  group_by(site, catchment, species, lifeform) %>%
  summarise_at(vars(count, alive), tibble::lst(mean, median))
Indiv_data_2a <- Indiv_data_2 %>% group_by(lifeform) %>% count (species)
#Indiv_sp_diff <- full_join(Indiv_data_3a, Indiv_data_2a, by = "species")
#setdiff(Indiv_data_2a,Indiv_data_3a)
Indiv_data_2 %>% group_by(site) %>%count ()
Indiv_data_3 %>% group_by(site) %>%count ()
cma_sp <- beltdata24_reveg_3 %>% 
  group_by(species, catchment) %>%
  summarise(number = n()) %>% tally(sort = T)

#calculate the mean and CI for alive and total
spdata_sp <- Indiv_data_3 %>%
  group_by(species, lifeform) %>%
  summarise(n = n(), 
            meantotal = mean(count_mean),
            setotal = sd(count_mean)/sqrt(n()),
            meanalive = mean(alive_mean),
            sealive = sd(alive_mean)/sqrt(n()))
Rich <- spdata_sp %>%
  mutate(surv = meanalive / meantotal)
#remove all species recorded at only 1 CMA
#Rich <- filter(Rich, n != '1')
#Rich <- filter(Rich, n != '2')
# order by abundance
spdata_sp <- Rich %>%
  arrange(desc(-surv))

#combine the 2 tables by species and arrange by CMA
Species_data <- full_join(spdata_sp, cma_sp, by = "species") %>%
  arrange(desc(n.y))

# cut the first 20 rows
Data_sp16 <-Species_data[1:20,]
tail(Data_sp16)
Data_sp16 <- dplyr::select (Data_sp16,-c(n.y))
Data_sp16 <- dplyr::select (Data_sp16,-c(surv))
#which(is.na(Sp_reshapeses))

# reshaping of data (this is pretty messy sorry)
Sp_reshape <- gather(Data_sp16, key = "variable", value = "mean", -c(species, lifeform, n.x, setotal, sealive))
Sp_reshape <- subset(Sp_reshape, select = -c(setotal, sealive))
Sp_reshapeses <- gather(Data_sp16, key = "variable", value = "error", -c(species, lifeform, n.x, meantotal, meanalive))
Sp_reshape$error <- Sp_reshapeses$error
Sp_reshape <- as.data.frame(Sp_reshape)

## Question - how do I include CI in this plot
Sp_reshape$variable <- factor(Sp_reshape$variable, levels = c("meantotal", "meanalive"))
Sp_reshape$lifeform <- as.factor(Sp_reshape$lifeform)

#reorder by lifeform, variable and mean (heighest to lowest)
Sp_reshape <- Sp_reshape[with(Sp_reshape, order(lifeform, variable, -mean)),]
Sp_reshape$species <- factor(Sp_reshape$species, levels = unique(Sp_reshape$species))
#ggsave(path = path, width = width, height = height, device='jpeg', dpi=300)

fig1 <- ggplot(Sp_reshape, aes(x = species, y = mean, fill = variable)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-error, ymax=mean+error), 
                width=.2,
                position=position_dodge(.9)) +
  facet_grid(~ lifeform, scales = "free_x", space = "free") +
  scale_fill_grey(start = 0.8, end = 0.6, name="Individual plant counts", breaks=c("meantotal", "meanalive"),
                  labels=c("Total plants", "Alive plants")) + 
  #scale_x_discrete(labels=c("Ozothamnus ferrugineus" = "Ozo. ferrugineus", "Leptospermum continentale" = "Lepto. continentale",
  #                          "Leptospermum lanigerum" = "Lepto. lanigerum", "Cassinia aculeata subsp. aculeata" = "Cassinia aculeata", 
  #                          "Allocasuarina verticillata" = "Allo. verticillata", 
  #                          "Acacia verniciflua s.l." = "Acacia verniciflua")) +
  labs(y = "Mean abundance per site", x = "Species") +
  theme_classic() + theme(legend.position = c(0.9, 0.8)) +
  theme(axis.text.x = element_text(size = 9, angle = 55, hjust = 1),
        axis.text.y = element_text(size = 9), axis.title.x = element_text(size=9),
        axis.title.y = element_text(size=9), legend.title=element_text(size=9), 
        legend.text=element_text(size=9))
fig1

###SURVIVAL ANALYSIS
survdat_3 <- beltdata24_reveg_3 %>% 
  group_by(site, catchment, species, lifeform) %>% summarise_at(vars(count, alive), tibble::lst(sum))
survdat_2 <- beltdata24_reveg_2 %>% 
  group_by(site, catchment, species, lifeform) %>% summarise_at(vars(count, alive), tibble::lst(sum))
survdat_3 <- survdat_3 %>% mutate(percent3=alive_sum/count_sum)
survdat_2 <- survdat_2 %>% mutate(percent2=alive_sum/count_sum)
survdat_3$catchment <- as.factor(survdat_3$catchment)
survdat_3$site <- as.factor(survdat_3$site)
#Indiv_data_3<- Indiv_data_3 %>% mutate(percent3=alive_mean/count_mean)
surv1_mod <- glmmTMB((percent2 + 0.01)/1.02 ~ catchment + (1|site), data = survdat_2, family = 'beta_family')
summary(surv1_mod)

site_sp_data <- left_join(survdat_3, survdat_2, join_by(site, species, catchment, lifeform))
site_sp_data<- site_sp_data %>% mutate(percent.x = alive_sum.x/count_sum.y)
site_sp_data <- site_sp_data %>% na.omit()
site_sp_data$catchment <- as.factor(site_sp_data$catchment)
site_sp_data$site <- as.factor(site_sp_data$site)
surv2_mod <- glmmTMB((percent3 + 0.01)/1.02 ~ catchment + (1|site), data = site_sp_data, family = 'binomial')
summary(surv2_mod)

##########Species Richness survival
site_sprich_visit3 <- beltdata24_reveg_3 %>%
  group_by(site, belt, catchment) %>%
  summarise(
    richness = n_distinct(species))
site_sprich_visit2 <- beltdata24_reveg_2 %>%
  group_by(site, belt, catchment) %>%
  summarise(
    richness = n_distinct(species))
site_sprich_dat <- left_join(site_sprich_visit3, site_sprich_visit2, join_by(site, belt, catchment))
site_sprich_dat <- site_sprich_dat %>% na.omit()
site_sprich_dat <- site_sprich_dat %>% mutate(perc_rich=richness.x/richness.y)
site_sprich_dat$perc_rich[site_sprich_dat$perc_rich >= 1] <- 1
surv2_mod <- glmmTMB(perc_rich ~ catchment + (1|site), data = site_sprich_dat, family = 'binomial')
summary(surv2_mod)
r.squaredGLMM(surv2_mod)
plot_model(surv2_mod)
coef(surv2_mod)

mean_survrich <- site_sprich_dat%>% 
  group_by(site, belt, catchment) %>% summarise_at(vars(richness.x, richness.y), tibble::lst(mean))
