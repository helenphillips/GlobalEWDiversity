createSplits <- function(dat, kfold = 10){
  rows <- nrow(dat)
  rows <- sample(rows, size = length(1:rows))
  
  chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  
  splits <- chunk(rows, kfold)
  
  return(splits)
}

calculateMSE <- function(data)
{
  cat("Function only works if the predicted and observed values are in their original units - NOT LOGGED \n")
  obsMinuspredsqr <- (data$observed - data$predicted)^2
  meanMSE <- mean(obsMinuspredsqr)
  
  return(meanMSE)
  
}


calculateMSEofQuantiles <- function(data, quantProbs = c(0.33, 0.66)){
  
  cat("Function only works if the predicted and observed values are in their original units - NOT LOGGED \n")
  obsMinuspredsqr <- (data$observed - data$predicted)^2
  low <- quantile(data$observed, probs = c(quantProbs[1], quantProbs[2]))[[1]]
  high <- quantile(data$observed, probs = c(quantProbs[1], quantProbs[2]))[[2]]
  
  data$quant <- 1
  data$quant <- ifelse(data$observed > low, 2, data$quant)
  data$quant <- ifelse(data$observed > high, 3, data$quant)
  
  
  MSE_low <- mean(obsMinuspredsqr[data$quant == 1])
  MSE_med <- mean(obsMinuspredsqr[data$quant == 2])
  MSE_high <- mean(obsMinuspredsqr[data$quant == 3])
  
  allMSEs = c( MSE_low , MSE_med ,MSE_high)
  names( allMSEs ) = c( "low" , "medium", "high" )
  
  return(allMSEs)
}


df_variables_sensitivity <- function(data){
  toMatch <- c("^PHIHOX$",
               "^CLYPPT$", "^SLTPPT$",
               "^CECSOL$", "^ORCDRC$",
               
               "^bio10_1$", "^bio10_4$",
               "^bio10_7$", "^bio10_12$",
               "^bio10_15$",
               
               "^Aridity$", "^PETyr$", "^PET_SD$", "^elevation$")
  
  matches <- unique (grep(paste(toMatch,collapse="|"), 
                          names(data), value=FALSE))
  
  return(matches)
  
}