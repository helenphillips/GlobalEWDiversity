df_variables <- function(data){
toMatch <- c("phFinal",
"ClayFinal", "SiltFinal",
"^CECSOL$", "OCFinal",

"^bio10_1$", "^bio10_4$",
"^bio10_7$", "^bio10_12$",
"^bio10_15$",

"^Aridity$", "^PETyr$", "^PET_SD$")

matches <- unique (grep(paste(toMatch,collapse="|"), 
                        names(data), value=FALSE))

return(matches)

}

findVariables <- function(df, VIFThreshold = 3){
  for(v in 1:ncol(df)){
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

