
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
date <- max(file_dates)
loadin <- files[grep(date, files)]
loadinspecies <- loadin[grep("species_", loadin)]

species <- read.csv(file.path(data_in, loadinspecies))



data_in <-"3_Data"

files <- list.files(file.path(data_in))
file_dates <- sapply(strsplit(files, "_"), "[", 2) ## Split the string by date, which produces a list, then take second element of each list i.e. the date
file_dates <- sapply(strsplit(file_dates, "\\."), "[", 1) ## Split the string by date, which produces a list, then take first element of each list i.e. the date

file_dates <- as.Date(file_dates)
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


table(fg_abundance$LU_Mgmt)

notEnough <- c("Integrated systems", "Tree plantations", "Unknown", "Urban")
fg_abundance <- droplevels(fg_abundance[!(fg_abundance$LU_Mgmt %in% notEnough),])


abund1epi <- lmer(logEpiAbundance ~ scalePH * LU_Mgmt + 
                 (1|file/Study_Name), data = fg_abundance)
abund2epi <- update(abund1epi, .~. -scalePH:LU_Mgmt)
anova(abund1epi, abund2epi) # Significant

abund1endo <- lmer(logEndoAbundance ~ scalePH * LU_Mgmt + 
                    (1|file/Study_Name), data = fg_abundance)
abund2endo <- update(abund1endo, .~. -scalePH:LU_Mgmt)
anova(abund1endo, abund2endo) ## Significant

abund1ane <- lmer(logAneAbundance ~ scalePH * LU_Mgmt + 
                     (1|file/Study_Name), data = fg_abundance)
abund2ane <- update(abund1ane, .~. -scalePH:LU_Mgmt)
anova(abund1ane, abund2ane) ## Not Significant
abund3ane <- update(abund2ane, .~. -scalePH)
abund4ane <- update(abund2ane, .~. -LU_Mgmt)
anova(abund2ane, abund3ane) # Significant
anova(abund2ane, abund4ane) # Significant

####################################################################
## A PLOT
####################################################################
ref <- 0 ## Getting mean value of pH


epi <- createNewdata(model= abund1epi, modelFixedEffects = c("LU_Mgmt", "scalePH"), data = fg_abundance)
epi <- predictValues(model = abund1epi, newdata = epi, responseVar = "logEpiAbundance", re.form = NA, seMultiplier = 1)
epi <- epi[which(abs(epi$scalePH-ref)==min(abs(epi$scalePH-ref))),]
epi <- epi[c(4, 5, 2, 1, 3),]
levels(epi$LU_Mgmt)[levels(epi$LU_Mgmt) == "Annual crop"] <- "Annual crops"

endo <- createNewdata(model= abund1endo, modelFixedEffects = c("LU_Mgmt", "scalePH"), data = fg_abundance)
endo <- predictValues(model = abund1endo, newdata = endo, responseVar = "logEndoAbundance", re.form = NA, seMultiplier = 1)
endo <- endo[which(abs(endo$scalePH-ref)==min(abs(endo$scalePH-ref))),]
endo <- endo[c(4, 5, 2, 1, 3),]

ane <- createNewdata(model= abund2ane, modelFixedEffects = c("LU_Mgmt", "scalePH"), data = fg_abundance)
ane <- predictValues(model = abund2ane, newdata = ane, responseVar = "logAneAbundance", re.form = NA, seMultiplier = 1)
ane <- ane[which(abs(ane$scalePH-ref)==min(abs(ane$scalePH-ref))),]
ane <- ane[c(4, 5, 2, 1, 3),]

ymin <- min(c(epi$lower, endo$lower, ane$lower))
ymax <- max(c(epi$upper, endo$upper, ane$upper))

fg_cols <- c("#FFCCCC", "#FF0000", "#990033")

## Plot
# pdf(file = file.path(figures, "Biomass_PastureIntensity.pdf"), height = 4)
jpeg(file = file.path(figures, "Abundance_functionalgroups.jpg"), quality = 100, res = 200, height = 1000, width = 1500)

par(mar=c(10, 4, 1, 2))
plot(-1e+05, -1e+05, ylim = c(ymin, ymax),
     xlim = c(0, (nrow(epi)-1)),  ylab = "log(Abundance)", xlab = " ",  xaxt='n', axes = FALSE)
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
       
       
       
axis(side=1, at = 0:(nrow(epi)-1), labels = epi$LU_Mgmt, las = 2)
dev.off()


