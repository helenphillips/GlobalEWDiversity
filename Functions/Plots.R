plotInteraction <- function(model, Effect1, Effect2, responseVar, seMultiplier = 1.96, data, cols = "000000", legend.position = "topleft", ylabel = "", xlabel = ""){
  # ci <- confint(model, method="Wald")
  
  ## Create values to predict over
  if(is.numeric(data[,which(names(data) == Effect1)])){
    eff1 <- seq(min(data[,which(names(data) == Effect1)], na.rm = TRUE), max(data[,which(names(data) == Effect1)], na.rm = TRUE), length.out = 100)
  } else { eff1 <- levels(model@frame[,which(names(model@frame) == Effect1)])}
  
  if(is.numeric(data[,which(names(data) == Effect2)])){
    eff2 <- seq(min(data[,which(names(data) == Effect2)], na.rm = TRUE), max(data[,which(names(data) == Effect2)], na.rm = TRUE), length.out = 100)
  } else { eff2 <- levels(model@frame[,which(names(model@frame) == Effect2)])}
  
  newdata <- with(data, expand.grid(E1 = eff1, E2 = eff2))
  names(newdata) <- c(Effect1, Effect2)
  
  ## Remove continuous variables not represented by data
  if(is.numeric(data[,which(names(data) == Effect1)]) && !(is.numeric(data[,which(names(data) == Effect2)]))){
    for(l in levels(newdata[,which(names(newdata) == Effect2)])){
      tmp <- data[data[which(names(data) == Effect2)] == l,]
      remove <- c(intersect(which(newdata[,1] < min(tmp[,which(names(data) == Effect1)], na.rm = TRUE)-0.1), which(newdata[,2] == l)),
                  intersect(which(newdata[,1] > max(tmp[,which(names(data) == Effect1)], na.rm = TRUE)+0.1), which(newdata[,2] == l)))
      if(length(remove) > 0){ ## Otherwise entire dataframes are removed if none of the points meet the criteria
        newdata <- newdata[-remove,]
      }
    }
  }
  
  if(is.numeric(data[,which(names(data) == Effect2)]) && !(is.numeric(data[,which(names(data) == Effect1)]))){
    for(l in levels(newdata[,which(names(newdata) == Effect1)])){
      tmp <- data[data[which(names(data) == Effect1)] == l,]
      remove <- c(intersect(which(newdata[,1] < min(tmp$scalePH, na.rm = TRUE)), which(newdata[,2] == l)),
                  intersect(which(newdata[,1] > max(tmp$scalePH, na.rm = TRUE)), which(newdata[,2] == l)))
      if(length(remove) > 0){ ## Otherwise entire dataframes are removed if none of the points meet the criteria
        newdata <- newdata[-remove,]
      }
    }
  }
  
  newdata$Response <- predict(model,newdata, re.form = NA)
  names(newdata)[which(names(newdata) == "Response")] <- responseVar
  mm <- model.matrix(terms(model), newdata)
  
  coef.names<-names(fixef(model))
  original.mm.names<-dimnames(mm)[[2]]
  if(length(coef.names)!=length(original.mm.names)){mm<-mm[,dimnames(mm)[[2]]%in%coef.names]} 
  
  
  pvar1 <- diag(mm %*% base::tcrossprod(as.matrix(vcov(model)), mm))
  
  newdata$upper <- newdata[,which(names(newdata) == responseVar)] + seMultiplier * sqrt(pvar1)
  newdata$lower <- newdata[,which(names(newdata) == responseVar)] - seMultiplier * sqrt(pvar1)
  
  
  if(length(which(c(class(newdata[,1]), class(newdata[,2])) == "factor")) != 1){
    stop("Plotting only works for one factor level and one continuous at the moment. Sorry")
  }
  
  
  #######################################
  ## One factor level and one continuous
  ########################################
  fac <- levels(newdata[,which(c(class(newdata[,1]), class(newdata[,2])) == "factor")])
  f <- which(c(class(newdata[,1]), class(newdata[,2])) == "factor")
  n <- which(c(class(newdata[,1]), class(newdata[,2])) == "numeric")
  r <- which(names(newdata) == responseVar)
  
  if(length(cols) != length(fac)){
    print("The length of colours does not match the number of factor levels. Some may be removed or duplicated")
  }
  cols <- rep(cols, length.out = length(fac))
  pt_cols <- paste("#", cols, sep="")
  
  ci_cols <- adjustcolor(pt_cols, alpha.f = 0.3)
  
  
  plot(-1e+05, -1e+05, ylim = c(min(newdata$lower,na.rm = TRUE), max(newdata$upper, na.rm = TRUE)),
       xlim = c(min(newdata[,n],na.rm = TRUE), max(newdata[,n], na.rm = TRUE)),  ylab = ylabel, xlab = xlabel)

  X.Vec <- c(newdata[,n][newdata[,f] == fac[1]], max(newdata[,n][newdata[,f] == fac[1]]), 
             rev(newdata[,n][newdata[,f] == fac[1]]), min(newdata[,n][newdata[,f] == fac[1]]))
  Y.Vec <- c(newdata$lower[newdata[,f] == fac[1]], tail(newdata$upper[newdata[,f] == fac[1]], 1), rev(newdata$upper[newdata[,f] == fac[1]]), newdata$lower[newdata[,f] == fac[1]][1])
  
  polygon(X.Vec, Y.Vec, col = ci_cols[1], border = NA)
  
  points(newdata[,n][newdata[,f] == fac[1]], newdata[,r][newdata[,f] == fac[1]], 
              col = pt_cols[1],  type = "l", lwd = 2)
  
  for(lev in 2:length(fac)){
    
    X.Vec <- c(newdata[,n][newdata[,f] == fac[lev]], max(newdata[,n][newdata[,f] == fac[lev]]), 
               rev(newdata[,n][newdata[,f] == fac[lev]]), min(newdata[,n][newdata[,f] == fac[lev]]))
    Y.Vec <- c(newdata$lower[newdata[,f] == fac[lev]], tail(newdata$upper[newdata[,f] == fac[lev]], 1), rev(newdata$upper[newdata[,f] == fac[lev]]), newdata$lower[newdata[,f] == fac[lev]][1])
    
    polygon(X.Vec, Y.Vec, col = ci_cols[lev], border = NA)
    
    points(newdata[,n][newdata[,f] == fac[lev]], newdata[,r][newdata[,f] == fac[lev]], 
           type = "l", lwd = 2, col = pt_cols[lev])
    
  }
  legend(legend.position, legend = fac, col = pt_cols, lwd = 2, bty = "n")
}
