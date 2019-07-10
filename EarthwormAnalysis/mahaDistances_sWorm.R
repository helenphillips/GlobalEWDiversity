## Computation of the distances

# The mahnalobis distance gives an indication of how far is a given point from the distribution.
# Create a data frame that only contains the predictor variables.
# This data frame is called predictorsData

setwd("C:/Research/Projects/sWorm/Files")

# Load you data here:
data = read.csv("richness.csv",header=TRUE)
summary(data)

# Removal of no data (this is only valuable if you populate your dataset directly from the maps)
data[data == -9999] <- NA
#data = na.omit(data)
summary(data)

# Create a dataset in a given order - Land Cover is not included here because we are working with averages and LC is a categorical variable
predictorsData = data.frame(ID = data$ID,
                            bio10_1 = scale(data$bio10_1),
bio10_4  = scale(data$bio10_4),
bio10_7  = scale(data$bio10_7),
bio10_12  = scale(data$bio10_12),
bio10_15  = scale(data$bio10_15),
cation = scale(data$cation),
elevation  = scale(data$elevation),
ai           = scale(data$ai),
ph           = scale(data$ph),
pet          = scale(data$pet),
pet_sd       = scale(data$pet_sd),
snow         = scale(data$snow),
clay         = scale(data$clay),
silt         = scale(data$silt),
carbon       = scale(data$carbon))

summary(predictorsData[,-1])


# These two lines compute the covariance matrix and the vector of means
# of the predictors. This are a 15x15 matrix and a vector with 15 entries.
sigma<-cov(predictorsData[,-1])
mu<-colMeans(predictorsData[,-1]) # in both cases take care what columns are you eliminating (IDs and Coordinates should not be here)


# Read the file with the new input variables - To be mapped - The shapefile that I am going to send you
distance=readOGR("Dataset_072019", layer = "point_grid") 
# distance <- read.csv("Dataset_072019\\data_globe_v2.txt",header=TRUE)
summary(distance)

# Create a new dataframe - keep the same order
distanceData = with(data, data.frame(distance$POINTID))
distanceData$bio10_1      <- distance$bio10_1
distanceData$bio10_4      <- distance$bio10_4
distanceData$bio10_7      <- distance$bio10_7
distanceData$bio10_12     <- distance$bio10_12
distanceData$bio10_15     <- distance$bio10_15
distanceData$cation       <- distance$cation
distanceData$elevation    <- distance$elevation
distanceData$ai           <- distance$ai
distanceData$ph           <- distance$ph
distanceData$pet          <- distance$pet
distanceData$pet_sd       <- distance$pet_sd
distanceData$snow         <- distance$snow
distanceData$clay         <- distance$clay
distanceData$silt         <- distance$silt
distanceData$carbon       <- distance$carbon

# Remove NAs if existing  (this method does not cope with NAs)
distanceData[distanceData == -9999] <- NA
distanceData = na.omit(distanceData)
summary(distanceData[,-1])
d = as.data.frame(distanceData[,-1]) # create a subset dataset to be easier to handle - not critical and you can remove it if you prefer

# Create a new column in the dataframe that saves the distance for each point.
distanceData$mahaDistance <- mahalanobis(d,mu,sigma) # Calculate the multidimentional distance
summary(distanceData)



length(which(distanceData$mahaDistance < qchisq(.975, df=15))) # How many points are in the 97.5% quantile
100 * (length(which(distanceData$mahaDistance < qchisq(.975, df=15))) / length(distanceData$distance.POINTID)) # estimating the proportion of pixels included (please note that the statistical definition of an outlier here is .975)

distanceData_final = with(data, data.frame(distanceData$distance.POINTID))
distanceData_final$mahaDist = distanceData$mahaDistance

write.csv(distanceData_final,"mahaDistances_richness_v1.csv") # Write the file with the distance.
