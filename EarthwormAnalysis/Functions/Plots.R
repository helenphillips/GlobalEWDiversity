createNewdata <- function(model, modelFixedEffects, data){
  
  ## Create Values
  vars <- list()
  for(e in 1:length(modelFixedEffects)){
    if(is.numeric(data[,which(names(data) == modelFixedEffects[e])])){
      vars[[e]] <- seq(min(data[,which(names(data) == modelFixedEffects[e])], na.rm = TRUE), max(data[,which(names(data) == modelFixedEffects[e])], na.rm = TRUE), length.out = 100)
    } else { vars[[e]] <- levels(model@frame[,which(names(model@frame) == modelFixedEffects[e])])}
    
  }

  
  newdata <- expand.grid(vars)
  names(newdata) <- modelFixedEffects
  
  
  return(newdata)
}

predictValues <- function(model, newdata, responseVar, re.form = NA, seMultiplier){
  newdata$Response <- predict(model,newdata, re.form = NA)
  names(newdata)[which(names(newdata) == "Response")] <- responseVar
  mm <- model.matrix(terms(model), newdata)
  
  coef.names<-names(fixef(model))
  original.mm.names<-dimnames(mm)[[2]]
  if(length(coef.names)!=length(original.mm.names)){mm<-mm[,dimnames(mm)[[2]]%in%coef.names]} 
  
  
  pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)), mm))
  
  newdata$upper <- newdata[,which(names(newdata) == responseVar)] + seMultiplier * sqrt(pvar1)
  newdata$lower <- newdata[,which(names(newdata) == responseVar)] - seMultiplier * sqrt(pvar1)
  
  return(newdata)
}

#################################################################################################
##
#################################################################################################


plotSingle <- function(model, modelFixedEffs, Effect1, responseVar, seMultiplier = 1.96, data, cols = "000000", legend.position, 
                       ylabel = "", xlabel = "", otherContEffectsFun = "median"){
  
  newdata <- createNewdata(model = model, modelFixedEffects = modelFixedEffs, data = data)
  
  ### Get rid of any unnesecary variable levels
  ## From model, which columns are additional variables?
  othervars <- model@frame[,names(model@frame) %in% modelFixedEffs[!(modelFixedEffs %in% Effect1)]]
  for(col in 1:ncol(othervars))
  {
    a <- which(names(newdata) == names(othervars)[col])
    
    if(names(othervars)[col] == "scalePH"){
      ref <- 0
      newdata <- newdata[which(abs(newdata[,a]-ref)==min(abs(newdata[,a]-ref))),]
      
    }
    
    if(class(othervars[,col]) == "factor"){
      # Use the reference level
      ref <- levels(data[,names(data) == names(othervars)[col]])[1]
      
      newdata <- newdata[newdata[,a] == ref,]
      
    }
    if(class(othervars[,col]) == "numeric"){
      f <- get(otherContEffectsFun)
      ref <- f(othervars[,col])
     
      newdata <- newdata[which(abs(newdata[,a]-ref)==min(abs(newdata[,a]-ref))),]
     
    }
    
  }
  
  ## Predict response and se values
  newdata <- predictValues(model, newdata, responseVar, re.form = NA, seMultiplier = seMultiplier)
  
  
  r <- which(names(newdata) == responseVar)
  p <- which(names(newdata) == Effect1)
  
  
  #####COLOURS
  ####
  

  
  
  ## If predictor is numeric
  if(class(newdata[,p]) == "numeric"){
    pt_cols <- paste("#", cols, sep="")
    
    ci_cols <- adjustcolor(pt_cols, alpha.f = 0.3)
    
    plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
       xlim = c(min(newdata[,p],na.rm = TRUE), max(newdata[,p], na.rm = TRUE)),  ylab = ylabel, xlab = xlabel)
    
    X.Vec <- c(newdata[,p], max(newdata[,p]), 
               rev(newdata[,p]), min(newdata[,p]))
    Y.Vec <- c(newdata$lower, tail(newdata$upper, 1), rev(newdata$upper), newdata$lower[1])
    
    polygon(X.Vec, Y.Vec, col = ci_cols, border = NA)
    
    points(newdata[,p], newdata[,r], 
           col = pt_cols,  type = "l", lwd = 5)
    
    
    }
  
  ## If predictor is factor
  if(class(newdata[,p]) == "factor"){
    
    if(length(cols) != length(levels(newdata[,p]))){
      print("The length of colours does not match the number of factor levels. Some may be removed or duplicated")
    }
    cols <- rep(cols, length.out = length(levels(newdata[,p])))
    pt_cols <- paste("#", cols, sep="")
    
    
    par(mar=c(10, 3, 1, 1))
    plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
         xlim = c(0, (nrow(newdata)-1)),  ylab = ylabel, xlab = xlabel,  xaxt='n', axes = FALSE)
    Axis(side = 2 )
    errbar(0:(nrow(newdata)-1), newdata[,r], newdata$upper, newdata$lower,
           add = TRUE, col = pt_cols, errbar.col = pt_cols, cex = 1.5)
    axis(side=1, at = 0:(nrow(newdata)-1), labels = levels(newdata[,p]), las=2)
    
    }
  # legend(legend.position, legend = fac, col = pt_cols, lwd = 2, bty = "n")
}


#################################################################################################
##
#################################################################################################



plotInteraction <- function(model, modelFixedEffs, Effect1, Effect2, responseVar, seMultiplier = 1.96, data, 
                            otherContEffectsFun = "median", cols = "000000", 
                            legend.position = "topleft", ylabel = "", xlabel = ""){
  # ci <- confint(model, method="Wald")
  
  # Create newdata 
  newdata <- createNewdata(model = model, modelFixedEffects = modelFixedEffs, data = data)
  
  
  ## Which column on newdata is Effect1
  ndEff1 <- which(names(newdata) == Effect1)
  ## Which column of newdata is Effect2
  ndEff2 <- which(names(newdata) == Effect2)
  
  ## Remove continuous variables not represented by data
  if(is.numeric(data[,which(names(data) == Effect1)]) && !(is.numeric(data[,which(names(data) == Effect2)]))){
    for(l in levels(newdata[,which(names(newdata) == Effect2)])){
      tmp <- data[data[which(names(data) == Effect2)] == l,]
      remove <- c(intersect(which(newdata[,ndEff1] < min(tmp[,which(names(data) == Effect1)], na.rm = TRUE)-0.1), which(newdata[,ndEff2] == l)),
                  intersect(which(newdata[,ndEff1] > max(tmp[,which(names(data) == Effect1)], na.rm = TRUE)+0.1), which(newdata[,ndEff2] == l)))
      if(length(remove) > 0){ ## Otherwise entire dataframes are removed if none of the points meet the criteria
        newdata <- newdata[-remove,]
      }
    }
  }
  
  if(is.numeric(data[,which(names(data) == Effect2)]) && !(is.numeric(data[,which(names(data) == Effect1)]))){
    for(l in levels(newdata[,which(names(newdata) == Effect1)])){
      tmp <- data[data[which(names(data) == Effect1)] == l,]
      remove <- c(intersect(which(newdata[,ndEff2] < min(tmp$scalePH, na.rm = TRUE)), which(newdata[,ndEff1] == l)),
                  intersect(which(newdata[,ndEff2] > max(tmp$scalePH, na.rm = TRUE)), which(newdata[,ndEff1] == l)))
      if(length(remove) > 0){ ## Otherwise entire dataframes are removed if none of the points meet the criteria
        newdata <- newdata[-remove,]
      }
    }
  }
  
  ## Predict response and variance
  newdata <- predictValues(model, newdata, responseVar, seMultiplier = seMultiplier, re.form = NA)
  
  # if(length(which(c(class(newdata[,1]), class(newdata[,2])) == "factor")) != 1){
  #  stop("Plotting only works for one factor level and one continuous at the moment. Sorry")
  # }
  
  
  #####################################
  ## What other variables are there, where we want the median or the reference?
  ########################################
  otherContEffectsFun = "median"
  
  othervars <- as.data.frame(model@frame[,intersect(which(names(model@frame) %in% modelFixedEffs), which(!names(model@frame) %in% c(Effect1, Effect2)))]  )
  names(othervars) <- names(model@frame)[intersect(which(names(model@frame) %in% modelFixedEffs), which(!names(model@frame) %in% c(Effect1, Effect2)))]
  for(col in 1:ncol(othervars))
  {
    a <- which(names(newdata) == names(othervars)[col])
    
    if(class(othervars[,col]) == "factor"){
      # Use the reference level
      ref <- levels(data[,names(data) == names(othervars)[col]])[1]
      
      newdata <- newdata[newdata[,a] == ref,]
      
    }
    if(class(othervars[,col]) == "numeric"){
      f <- get(otherContEffectsFun)
      ref <- f(othervars[,col])
      
      newdata <- newdata[which(abs(newdata[,a]-ref)==min(abs(newdata[,a]-ref))),]
      
    }
    
  }
  
  

  r <- which(names(newdata) == responseVar)
  p1 <- which(names(newdata) == Effect1)
  p2 <- which(names(newdata) == Effect2)
  
  
  ###########################################################################
  ### With one categorical and one factorial effect
  ############################################################################
  
  if(any(c(class(newdata[,p1]) == "numeric", class(newdata[,p2]) == "numeric"))){
    
    fac <- ifelse(class(newdata[,p1]) == "factor", p1, p2)
    num <- ifelse(class(newdata[,p2]) == "factor", p1, p2)
    
    if(length(cols) != length(levels(newdata[,fac]))){
      print("The length of colours does not match the number of factor levels. Some may be removed or duplicated")
    }
    cols <- rep(cols, length.out = length(levels(newdata[,fac])))
    pt_cols <- paste("#", cols, sep="")
    
    ci_cols <- adjustcolor(pt_cols, alpha.f = 0.3)
    
    plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
         xlim = c(min(newdata[,num],na.rm = TRUE), max(newdata[,num], na.rm = TRUE)),  ylab = ylabel, xlab = xlabel)
    
    X.Vec <- c(newdata[,num][newdata[,fac] == levels(newdata[,fac])[1]], max(newdata[,num][newdata[,fac] == levels(newdata[,fac])[1]]), 
               rev(newdata[,num][newdata[,fac] == levels(newdata[,fac])[1]]), min(newdata[,num][newdata[,fac] == levels(newdata[,fac])[1]]))
    Y.Vec <- c(newdata$lower[newdata[,fac] == levels(newdata[,fac])[1]], tail(newdata$upper[newdata[,fac] == levels(newdata[,fac])[1]], 1), rev(newdata$upper[newdata[,fac] == levels(newdata[,fac])[1]]), newdata$lower[newdata[,fac] == levels(newdata[,fac])[1]][1])
    
    polygon(X.Vec, Y.Vec, col = ci_cols[1], border = NA)
    
    points(newdata[,num][newdata[,fac] == levels(newdata[,fac])[1]], newdata[,r][newdata[,fac] == levels(newdata[,fac])[1]], 
           col = pt_cols[1],  type = "l", lwd = 2)
    
    for(lev in 2:length(levels(newdata[,fac]))){
      
      
      X.Vec <- c(newdata[,num][newdata[,fac] == levels(newdata[,fac])[lev]], max(newdata[,num][newdata[,fac] == levels(newdata[,fac])[lev]]), 
                 rev(newdata[,num][newdata[,fac] == levels(newdata[,fac])[lev]]), min(newdata[,num][newdata[,fac] == levels(newdata[,fac])[lev]]))
      Y.Vec <- c(newdata$lower[newdata[,fac] == levels(newdata[,fac])[lev]], tail(newdata$upper[newdata[,fac] == levels(newdata[,fac])[lev]], 1), rev(newdata$upper[newdata[,fac] == levels(newdata[,fac])[lev]]), newdata$lower[newdata[,fac] == levels(newdata[,fac])[lev]][1])
      
      polygon(X.Vec, Y.Vec, col = ci_cols[lev], border = NA)
      
      points(newdata[,num][newdata[,fac] == levels(newdata[,fac])[lev]], newdata[,r][newdata[,fac] == levels(newdata[,fac])[lev]], 
             col = pt_cols[lev],  type = "l", lwd = 2)
      
      
      }
    legend(legend.position, legend = fac, col = pt_cols, lwd = 2, bty = "n")
    
  }
  
  ###########################################################################
  ### If both effects are categorical
  ############################################################################
  if(class(newdata[,p1]) == "factor" && class(newdata[,p1]) == "factor"){
    
    print("Plot will be grouped by Effect 1")
    
    if(length(cols) != length(levels(newdata[,ndEff1]))){
      print("The length of colours does not match the number of factor levels. Some may be removed or duplicated")
    }
    
    cols <- cols[1:length(levels(newdata[,ndEff1]))]
    
    cols <- rep(cols, each = length(levels(newdata[,ndEff2])))
    pt_cols <- paste("#", cols, sep="")
  
    newdata$Colour <- pt_cols
    
    #### Get rid of levels not represented by data
    
    tokeep <- c()
    for(l in levels(newdata[,ndEff2])){
      tmp <- data[data[which(names(data) == Effect2)] == l,]
      keep <- intersect(which(newdata[,ndEff2] == l), which(newdata[,ndEff1] %in% levels(droplevels(tmp$LU_Mgmt))))
      
      tokeep <- c(tokeep, keep)
      
    }
    newdata <- newdata[tokeep,]
    
  
    ##### Group by Effect1
    print("Grouping plot by Effect1")
    
    order <- c()
    for(i in (levels(data[,which(names(data) == Effect1)]))){
      order <- c(order, which(newdata[,ndEff1] == i))
    }
    
    newdata <- newdata[order,]
    
    plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
         xlim = c(0, nrow(newdata)),  ylab = ylabel, xlab = xlabel)
    
    errbar(1:nrow(newdata), newdata[,r], newdata$upper, newdata$lower,
           add = TRUE, col = newdata$Colour, errbar.col = newdata$Colour, cex = 1.5)
    text(1:nrow(newdata), rep(0, nrow(newdata)), labels = newdata$intensity) 
    # legend(legend.position, legend = levels(newdata[,ndEff1]), col = unique(newdata$Colour[order]), lwd = 2, bty = "n")
  }
  
  
  ###########################################################################
  ### If both effects are continuous
  ############################################################################
  
   if(class(newdata[,p1]) == "numeric" && class(newdata[,p2]) == "numeric"){
    stop("Plotting does not work for two continuous variables at the moment. Sorry")
   }
  
  
}
