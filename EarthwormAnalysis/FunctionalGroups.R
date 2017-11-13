
library(dplyr)
########################################################
# 0. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\Users\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}



##########################################################
## Load in data
##########################################################

data_in <-"0_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
date <- max(file_dates, na.rm = TRUE)
loadin <- files[grep(date, files)]
loadinspecies <- loadin[grep("species_", loadin)]

species <- read.csv(file.path(data_in, loadinspecies))



data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates, na.rm = TRUE)
date <- max(file_dates)
loadin <- files[grep(date, files)]

sites <- read.csv(file.path(data_in, loadin))

##########################################################
## Calculate site level measures of the three main functional groups
##########################################################


Summary.bio <- species %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_site) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    Epi_biomass = sum(WetBiomass[Functional_Type == "Epigeic"]),
    Endo_biomass = sum(WetBiomass[Functional_Type == "Endogeic"]),
    Ane_biomass = sum(WetBiomass[Functional_Type == "Anecic"])
  )

summary.bio <- as.data.frame(Summary.bio)
summary.bio <- summary.bio[complete.cases(summary.bio),]


Summary.abund <- species %>% # Start by defining the original dataframe, AND THEN...
  group_by(Study_site) %>% # Define the grouping variable, AND THEN...
  summarise( # Now you define your summary variables with a name and a function...
    Epi_abundance = sum(Abundance[Functional_Type == "Epigeic"]),
    Endo_abundance = sum(Abundance[Functional_Type == "Endogeic"]),
    Ane_abundance = sum(Abundance[Functional_Type == "Anecic"])
  )

summary.abund <- as.data.frame(Summary.abund)
summary.abund <- summary.abund[complete.cases(summary.abund),]

##########################################################
## Match with site level dataset
##########################################################

sites_fg <- (merge(sites, summary.bio, by = "Study_site", all.x = TRUE))
sites_fg <- (merge(sites_fg, summary.abund, by = "Study_site", all.x = TRUE))


##########################################################
## simple model
##########################################################
## Abundance has the most data

fg_abundance <- sites_fg[!(is.na(sites_fg$Epi_abundance)),]

hist(log(fg_abundance$Epi_abundance + 1))
fg_abundance$logEpiAbundance <- log(fg_abundance$Epi_abundance + 1)

fg_abundance$logEndoAbundance <- log(fg_abundance$Endo_abundance + 1)
fg_abundance$logAneAbundance <- log(fg_abundance$Ane_abundance + 1)


fg_abundance <- fg_abundance[complete.cases(fg_abundance$scalePH),]


table(fg_abundance$ESA)

Enough <- c("Broadleaf deciduous forest", "Herbaceous", "Production - Herbaceous",
            "Production - Plantation")
fg_abundance <- droplevels(fg_abundance[which(fg_abundance$ESA %in% Enough),])


fg_abundance$bio10_5_scaled <- scale(fg_abundance$bio10_5)
fg_abundance$bio10_13_scaled <- scale(fg_abundance$bio10_13)
fg_abundance$bio10_14_scaled <- scale(fg_abundance$bio10_14)

abund1epi <- lmer(logEpiAbundance ~ scalePH * ESA + 
                    (bio10_5_scaled + bio10_13_scaled + bio10_14_scaled)^2 + 
                 (1|file/Study_Name), data = fg_abundance)

a2epia <- update(abund1epi, .~. -bio10_13_scaled:bio10_14_scaled)
anova(abund1epi, a2epia) ## Not significant
a2epib <- update(abund1epi, .~. -bio10_5_scaled:bio10_14_scaled )
anova(abund1epi, a2epib) ## Not significant
a2epic <- update(abund1epi, .~. -bio10_5_scaled:bio10_13_scaled)
anova(abund1epi, a2epic) ## Not significant
a2epid <- update(abund1epi, .~. -scalePH:ESA)
anova(abund1epi, a2epid) ## significant

a3epia <- update(a2epia, .~. -bio10_5_scaled:bio10_14_scaled)
anova(a2epia, a3epia) #Not significant
a3epib <- update(a2epia, .~. -bio10_5_scaled:bio10_13_scaled)
anova(a2epia, a3epib) # significant
a3epic <- update(a2epia, .~. -scalePH:ESA)
anova(a2epia, a3epic) # significant


a4epia <- update(a3epia, .~. -bio10_5_scaled:bio10_13_scaled)
anova(a3epia, a4epia) ## Significant
a4epib <- update(a3epia, .~. -scalePH:ESA)
anova(a3epia, a4epib) ## Significant


a5epia <- update(a3epia, .~. -bio10_14_scaled)
anova(a3epia, a5epia) ## Not significant

####### 

abund1endo <- lmer(logEndoAbundance ~ scalePH * ESA + 
                     (bio10_5_scaled + bio10_13_scaled + bio10_14_scaled)^2 +
                    (1|file/Study_Name), data = fg_abundance)

a2endoa <- update(abund1endo, .~. -bio10_13_scaled:bio10_14_scaled)
anova(abund1endo, a2endoa) # Not sign
a2endob <- update(abund1endo, .~. -bio10_5_scaled:bio10_14_scaled)
anova(abund1endo, a2endob) # Sign
a2endoc <- update(abund1endo, .~. -bio10_5_scaled:bio10_13_scaled)
anova(abund1endo, a2endoc) # Not sign
a2endod <- update(abund1endo, .~. -scalePH:ESA)
anova(abund1endo, a2endod) # Sign


a3endoa <- update(a2endoc, .~. -bio10_13_scaled:bio10_14_scaled)
anova(a2endoc, a3endoa) ## Not sign
a3endob <- update(a2endoc, .~. -bio10_5_scaled:bio10_14_scaled)
anova(a2endoc, a3endob) # Sign
a3endoc <- update(a2endoc, .~. -scalePH:ESA)
anova(a2endoc, a3endoc) # Sign


a4endoa <- update(a3endoa, .~. - bio10_5_scaled:bio10_14_scaled)
anova(a3endoa, a4endoa) # Sign
a4endob <- update(a3endoa, .~. - scalePH:ESA)
anova(a3endoa, a4endob) # Sign

a5endoa <- update(a3endoa, .~. - bio10_13_scaled)
anova(a3endoa, a5endoa) # Sign

###############

abund1ane <- lmer(logAneAbundance ~ scalePH * ESA + 
                    (bio10_5_scaled + bio10_13_scaled + bio10_14_scaled)^2 + 
                     (1|file/Study_Name), data = fg_abundance)



a2anea <- update(abund1ane, .~. -bio10_13_scaled:bio10_14_scaled)
anova(abund1ane, a2anea) # Not sign.
a2aneb <- update(abund1ane, .~. -bio10_5_scaled:bio10_14_scaled)
anova(abund1ane, a2aneb) #Not sign.
a2anec <- update(abund1ane, .~. -bio10_5_scaled:bio10_13_scaled)
anova(abund1ane, a2anec) # Not sign.
a2aned <- update(abund1ane, .~. -scalePH:ESA)
anova(abund1ane, a2aned) #Sign


a3anea <- update(a2aneb, .~. -bio10_13_scaled:bio10_14_scaled)
anova(a2aneb, a3anea) # Not sign.
a3aneb <- update(a2aneb, .~. -bio10_5_scaled:bio10_13_scaled)
anova(a2aneb, a3aneb) # Not sign.
a3anec <- update(a2aneb, .~. -scalePH:ESA)
anova(a2aneb, a3anec) #Sign

a4anea <- update(a3anea, .~. -bio10_5_scaled:bio10_13_scaled)
anova(a3anea, a4anea) #  sign.
a4aneb <- update(a3anea, .~. -scalePH:ESA)
anova(a3anea, a4aneb) #Sign

a5anea <- update(a3anea, .~. -bio10_14_scaled)
anova(a3anea, a5anea) ## Sign

####################################################################
## A PLOT
####################################################################
ref <- 0 ## Getting mean value of pH, and other variables?


epi <- createNewdata(model= a5epia, modelFixedEffects = c("scalePH", "ESA", "bio10_5_scaled", "bio10_13_scaled","bio10_14_scaled"),
                     data = fg_abundance, mainEffect = "ESA")
epi <- predictValues(model = a5epia, newdata = epi, responseVar = "logEpiAbundance", re.form = NA, seMultiplier = 1)

endo <- createNewdata(model= a3endoa, modelFixedEffects = c("scalePH", "ESA", "bio10_5_scaled", "bio10_13_scaled","bio10_14_scaled"),
                      data = fg_abundance, mainEffect = "ESA")
endo <- predictValues(model = a3endoa, newdata = endo, responseVar = "logEndoAbundance", re.form = NA, seMultiplier = 1)

ane <- createNewdata(model= a3anea,  modelFixedEffects = c("scalePH", "ESA", "bio10_5_scaled", "bio10_13_scaled","bio10_14_scaled"),
                     data = fg_abundance, mainEffect = "ESA")
ane <- predictValues(model = a3anea, newdata = ane, responseVar = "logAneAbundance", re.form = NA, seMultiplier = 1)

ymin <- min(c(epi$lower, endo$lower, ane$lower))
ymax <- max(c(epi$upper, endo$upper, ane$upper))

fg_cols <- c("#FFCCCC", "#FF0000", "#990033")

## Plot
# pdf(file = file.path(figures, "Biomass_PastureIntensity.pdf"), height = 4)
jpeg(file = file.path(figures, "Abundance_functionalgroups.jpg"), quality = 100, res = 200, height = 1000, width = 1500)

par(mar=c(12, 4, 1, 2))
plot(-1e+05, -1e+05, ylim = c(ymin, ymax),
     xlim = c(0, (nrow(epi))),  ylab = "log(Abundance)", xlab = " ",  xaxt='n', axes = FALSE)
Axis(side = 2 )
errbar(0:(nrow(epi)-1), epi$logEpiAbundance, epi$upper, epi$lower,
       add = TRUE, col = "black", errbar.col = fg_cols[1], cex = 1.5)

errbar(0:(nrow(epi)-1), epi$logEpiAbundance, epi$upper, epi$lower,
       add = TRUE, col = fg_cols[1], errbar.col = fg_cols[1], cex = 1.3)


errbar(0:(nrow(endo)-1), endo$logEndoAbundance, endo$upper, endo$lower,
       add = TRUE, col = fg_cols[2], errbar.col = fg_cols[2], cex = 1.5)

errbar(0:(nrow(ane)-1), ane$logAneAbundance, ane$upper, ane$lower,
       add = TRUE, col = fg_cols[3], errbar.col = fg_cols[3], cex = 1.5)


legend(3.05, y = 4, legend = c("Epigeics", "Endogeics", "Anecics"), pch = 19, col = fg_cols, bty = "n", cex = 1.5)
       
       
       
axis(side=1, at = 0:(nrow(epi)-1), labels = epi$ESA, las = 2)
dev.off()


