findVariables <- function(df, VIFThreshold = 3){
  for(var in 1:ncol(df)){
    t <- corvif(df)
    maxVIF <- which(t == max(t))
    if(t[maxVIF,] > VIFThreshold){
      # remove col from df
      print(t[maxVIF,])
      print(paste("removing", names(df)[maxVIF]))
      df <- df[,-maxVIF]
    } else {
      break
    }
  }
  return(names(df))  
}

