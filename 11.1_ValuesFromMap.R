
if(Sys.info()["nodename"] == "TSGIS02"){
  setwd("C:/sWorm/EarthwormAnalysis")
}


library(data.table)


#results <- "D:\\sWorm\\Results\\Richness"
# results <- "D:\\sWorm\\Results\\Abundance"
results <- "D:\\sWorm\\Results\\Biomass"

# regions <- c("africa", "asia", "europe", "latin_america", "north_america", "west_asia")

resultCSV <- "predictedValues.csv"

asia <-  fread(file.path(results, "asia", resultCSV), sep = ",")
asia <- asia[complete.cases(asia),]
asia$V1 <- exp(asia$V1)
print(min(asia$V1, na.rm = TRUE))
print(max(asia$V1, na.rm = TRUE))
print(mean(asia$V1, na.rm = TRUE))
print(median(asia$V1, na.rm = TRUE))
print(sd(asia$V1, na.rm = TRUE))


africa <- fread(file.path(results, "africa", resultCSV), sep = ",")
africa <- africa[complete.cases(africa),]
africa$V1 <- exp(africa$V1)
print(min(africa$V1, na.rm = TRUE))
print(max(africa$V1, na.rm = TRUE))
print(mean(africa$V1, na.rm = TRUE))
print(median(africa$V1, na.rm = TRUE))
print(sd(africa$V1, na.rm = TRUE))

all <- rbind(africa, asia)
rm(africa); rm(asia)

europe <- fread(file.path(results, "europe", resultCSV), sep = ",")
europe <- europe[complete.cases(europe),]
europe$V1 <- exp(europe$V1) 
print(min(europe$V1, na.rm = TRUE))
print(max(europe$V1, na.rm = TRUE))
print(mean(europe$V1, na.rm = TRUE))
print(median(europe$V1, na.rm = TRUE))
print(sd(europe$V1, na.rm = TRUE))

all <- rbind(all, europe)
rm(europe)

gc()

latin_america <- fread(file.path(results, "latin_america", resultCSV), sep = ",")
latin_america <- latin_america[complete.cases(latin_america),]
latin_america$V1 <- exp(latin_america$V1) 
print(min(latin_america$V1, na.rm = TRUE))
print(max(latin_america$V1, na.rm = TRUE))
print(mean(latin_america$V1, na.rm = TRUE))
print(median(latin_america$V1, na.rm = TRUE))
print(sd(latin_america$V1, na.rm = TRUE))

all <- rbind(all, latin_america)
rm(latin_america)

north_america <- fread(file.path(results, "north_america", resultCSV), sep = ",")
north_america <- north_america[complete.cases(north_america),]
north_america$V1 <- exp(north_america$V1) 
print(min(north_america$V1, na.rm = TRUE))
print(max(north_america$V1, na.rm = TRUE))
print(mean(north_america$V1, na.rm = TRUE))
print(median(north_america$V1, na.rm = TRUE))
print(sd(north_america$V1, na.rm = TRUE))

all <- rbind(all, north_america)
rm(north_america)

gc()

west_asia <- fread(file.path(results, "west_asia", resultCSV), sep = ",")
west_asia <- west_asia[complete.cases(west_asia),]
west_asia$V1 <- exp(west_asia$V1)
print(min(west_asia$V1, na.rm = TRUE))
print(max(west_asia$V1, na.rm = TRUE))
print(mean(west_asia$V1, na.rm = TRUE))
print(median(west_asia$V1, na.rm = TRUE))
print(sd(west_asia$V1, na.rm = TRUE))

all <- rbind(all, west_asia)
rm(west_asia)

gc()

mean(all$V1)
sd(all$V1)
min(all$V1)
max(all$V1)
median(all$V1)

