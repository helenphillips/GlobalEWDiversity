########################################################
# 1. Set Working Directory
########################################################

if(Sys.info()["nodename"] == "IDIVNB193"){
  setwd("C:\\restore2\\hp39wasi\\sWorm\\EarthwormAnalysis\\")
}

if(Sys.info()["nodename"] == "IDIVNB179"){
  setwd("C:\\USers\\hp39wasi\\WORK\\sWorm\\EarthwormAnalysis\\")
}


########################################################
# 2. Load data
########################################################
source("Functions/loadMostRecent.R")



data_in_new <- "0_Data"
newdata <- loadMostRecent("sites_", data_in_new)
newdata <- read.csv(file.path(data_in_new, newdata))



data_in_old <- "8_Data"
## THIS WORKS FOR NOW, AS I HAVEN'T RE-RUN THE SCRIPT.
richness <- loadMostRecent("sitesRichness", data_in_old)
richness <- read.csv(file.path(data_in_old, richness))
abundance <- loadMostRecent("sitesAbundance_", data_in_old)
abundance <- read.csv(file.path(data_in_old, abundance))
biomass <- loadMostRecent("sitesBiomass_", data_in_old)
biomass <- read.csv(file.path(data_in_old, biomass))

########################################################
# 3. Richness
########################################################

richnessStudies <- unique(richness$Study_Name)
new_richness <- newdata[newdata$Study_Name %in% richnessStudies,]

nrow(new_richness) # 8133
length(unique(new_richness$Study_Name)) # 168

nrow(richness) # 5416
length(unique(richness$Study_Name)) # 172 ???

richnessStudies[!(richnessStudies %in% new_richness$Study_Name)]
# Moscow_1983   Satino_1      Kamenn_Step_1 coors2016a
## These are the studies with 1 site which I have now removed

old <- data.frame(table(richness$Study_Name))
names(old)[2] <- "oldN"
new <- data.frame(table(new_richness$Study_Name))
names(new)[2] <- "newN"
diff <- merge(old, new, by = "Var1", all.x = TRUE)
diff$diff <- diff$newN - diff$oldN
diff <- diff[complete.cases(diff),]

nrow(diff[diff$diff > 200,]) # 3
nrow(diff[diff$diff > 50,]) # 15
nrow(diff[diff$diff > 10,]) # 41
nrow(diff[diff$diff  == 0,]) # 80
########################################################
# 4. Abundance
########################################################

abundanceStudies <- unique(abundance$Study_Name)
new_abundance <- droplevels(newdata[newdata$Study_Name %in% abundanceStudies,])

nrow(new_abundance) # 9249
length(unique(new_abundance$Study_Name)) #  207

nrow(abundance) # 6388
length(unique(abundance$Study_Name)) # 211 ??? 1 NA

abundanceStudies[!(abundanceStudies %in% new_abundance$Study_Name)]
# Moscow_1983   Satino_1      Kamenn_Step_1 coors2016a
## These are the studies with 1 site which I have now removed

new_abundance$Study_Name[!(new_abundance$Study_Name %in% abundanceStudies)]


old <- data.frame(table(abundance$Study_Name))
names(old)[2] <- "oldN"
new <- data.frame(table(new_abundance$Study_Name))
names(new)[2] <- "newN"
diff <- merge(old, new, by = "Var1", all = TRUE)
diff$diff <- diff$newN - diff$oldN
diff <- diff[complete.cases(diff),]

nrow(diff[diff$diff > 200,]) # 3
nrow(diff[diff$diff > 50,]) # 17
nrow(diff[diff$diff > 10,]) # 44
nrow(diff[diff$diff  == 0,]) # 116

########################################################
# 5. Biomass
########################################################

biomassStudies <- unique(biomass$Study_Name)
new_biomass <- droplevels(newdata[newdata$Study_Name %in% biomassStudies,])

nrow(new_biomass) #  4333
length(unique(new_biomass$Study_Name)) #  125

nrow(biomass) # 3326
length(unique(biomass$Study_Name)) # 125  1 NA

# biomassStudies[!(biomassStudies %in% new_biomass$Study_Name)]


new_biomass$Study_Name[!(new_biomass$Study_Name %in% biomassStudies)]


old <- data.frame(table(biomass$Study_Name))
names(old)[2] <- "oldN"
new <- data.frame(table(new_biomass$Study_Name))
names(new)[2] <- "newN"
diff <- merge(old, new, by = "Var1", all = TRUE)
diff$diff <- diff$newN - diff$oldN
diff <- diff[complete.cases(diff),]

nrow(diff[diff$diff > 200,]) # 0 
nrow(diff[diff$diff > 50,]) # 7
nrow(diff[diff$diff > 10,]) # 22 
nrow(diff[diff$diff  == 0,]) # 73

########################################################
# 6. Manual Check
########################################################

usefulCols <- c("file", "Study_Name", "Site_Name", "SpeciesRichness",             
               "SpeciesRichnessUnit","Site_WetBiomass","Site_WetBiomassUnits",        
               "Site_Abundance","Site_AbundanceUnits")

## Studies I have found problems with

## 4013_Fusilero2013

t1 <- newdata[which(newdata$file == "4013_Fusilero2013"),] # 24 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]

## 4130_Kuntz2013
t1 <- newdata[which(newdata$file == "4130_Kuntz2013"),] # 4 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]

## 3331_Watmough2014
t1 <- newdata[which(newdata$file == "3331_Watmough2014"),] # 4 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]

## 7543_Bedano2016
t1 <- newdata[which(newdata$file == "7543_Bedano2016"),] # 250 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]

## 000_Baker1998b
t1 <- newdata[which(newdata$file == "000_Baker1998b"),] # 280 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]

## 46_Kernecker2015
t1 <- newdata[which(newdata$file == "46_Kernecker2015"),] # 12 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]

## 7601_GutierrezLopez2016
t1 <- newdata[which(newdata$file == "7601_GutierrezLopez2016"),] # 327 rows
t1 <- t1[,which(names(t1) %in% usefulCols)]
########################################################
# 7. Random Check
########################################################

# r_samp <- sample(newdata$file, size = 10)

smaller <- data.frame(table(newdata$file))
smaller <- smaller[smaller$Freq < 10,]

r_samp <- sample(smaller$Var1, size = 10, replace = FALSE)


for(r in 1:length(r_samp)){
  
  print(r_samp[r])
  t2 <- newdata[which(newdata$file == r_samp[r]),]
  
  t2 <- t2[,which(names(t2) %in% usefulCols)]
  
  print(nrow(t2))
  
  if(nrow(t2) < 20){
    print(t2)
    readline(prompt="Press [enter] to continue")
  }
  
  if(nrow(t2) > 20){
    print(head(t2, 10))
    readline(prompt="Press [enter] to see bottom of dataframe")
    
    print(tail(t2, 10))
    readline(prompt="Press [enter] to continue")
    
  }
  
}

#########################

t3 <- data.frame(table(newdata$file, newdata$Site_WetBiomassUnits))

t3 <- t3[t3$Freq > 0,]

t4 <- data.frame(table(newdata$file))


t5 <- merge(t3, t4, by = "Var1", all = TRUE)

t5 <- t5[!(is.na(t5$Freq.x)),]

t5$diff <- t5$Freq.y - t5$Freq.x
t5 <- t5[t5$diff > 0,]


#[1] 000_Baker1998b          000_Beausejour2015      000_Bescansa2010       
#[4] 000_MulderUnpublished   113_Davalos2015         350_Iannone2015        
#[7] 3532_Arai2014           3644_Hishi2014          4421_Loss2012          
#[10] 4643_Nieminen2011       4716_Loss2011           5184_Li2010            
#[13] 5597_Hirth2009          5907_Joschko2009        7471_Moos2016          
#[16] 7601_GutierrezLopez2016



# 000_Bescansa2010

newdata[newdata$file == "000_Bescansa2010",]
# NAs are expected

# 000_Beausejour2015

ba <- newdata[newdata$file == "000_Beausejour2015",]

