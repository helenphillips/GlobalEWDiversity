library(lme4)


modelSimplification <- function(model = model, optimizer = "bobyqa", Iters = 2e5,  data, alpha = 0.05){
  
  
  #cars$type <- as.factor(rep(letters[1:4], length = 50))
  #cars$rand <- as.factor(rep(letters[5:7], length = 50))
  
  #
  
  
  form <- as.character(model@call$formula)
  form <- form[-grep("~", form)]
  response <- form[-grep("\\+", form)]
  
  randEffect <- form[grep("\\+", form)]
  randEffect <- gsub("^.*\\(", "", randEffect) ## This also removes the first bracket
  randEffect <- paste0("+ (", randEffect)
  
  fam <- as.character(model@call$family)
  ## lmers don't have a family, so....
  if(length(fam) == 0) {fam <- "notpoisson"}
  
  
  all.terms <- rownames(anova(model))
  interactions <- all.terms[grep(":", all.terms)]
  main <- all.terms[!(all.terms %in% interactions)]
  
  res <- data.frame(term = NA, pVal = NA)
  
  if(fam == "poisson"){
    
    refModel <- glmer(formula = model@call$formula, data = data, family = fam, 
                      control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
    
  }else{
    refModel <- lmer(formula = model@call$formula, data = data, 
                     control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
  } ## How do I change this 
  
  
  if(length(interactions) > 0){
    while(length(interactions) > 0){
      ## break loop if remaining interactions have pValues greater than alpha
      # (and none of the pvalues are NA)
      if(all(!(is.na(res$pVal))) && max(res$pVal) < alpha)
      {
        break
      }
      
      
      res <- data.frame(term = NA, pVal = NA)
      
      for( inter in 1:length(interactions)){
        cat(paste("Testing interaction: ", interactions[inter], "\n", sep=""))
        used <- c(main, interactions[-inter])
        res[inter, 'term'] <- interactions[inter]
        
        fixedeffs <- paste(used,collapse="+")
        
        call <- paste(response, "~", fixedeffs, randEffect)
        
        if(fam == "poisson"){
          model2 <- glmer(formula = call, data = data, family = fam, 
                          control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
        }else{
          model2 <- lmer(formula = call, data = data,
                         control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
          
        }
        res[inter, 'pVal'] <- anova(refModel, model2)[2,'Pr(>Chisq)']
        print(res[inter, 'pVal'])
      }
      
      if( max(res$pVal, na.rm = TRUE) > alpha){ 
        toRemove <- res$term[which(res$pVal == max(res$pVal, na.rm = TRUE))]
        cat(paste("Removing", toRemove, "\n"))
        interactions <- interactions[-which(toRemove == interactions)]
      }
      # recreate a ref model
      used <- c(main, interactions)
      fixedeffs <- paste(used,collapse="+")
      call <- paste(response, "~", fixedeffs, randEffect)
      if(fam == "poisson"){
        refModel <- glmer(formula = call, data = data, family = fam, 
                          control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
      }else{
        refModel <- lmer(formula = call, data = data, 
                         control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
      }
    }
  }
  
  
  ## Deal with any main effects now
  
  ## TODO: What is no interactions to begin with
  
  
  ## What main effects are still in interactions
  matches <- unique (grep(paste(main,collapse="|"), 
                          interactions, value=TRUE))
  
  
  keep <- unique(unlist(strsplit(matches, "[:]")))
  
  left <- main[!(main %in% keep)]
  
  cat(paste("Testing removal of main effects\n"))
  
  inInteractions <- main[main %in% keep]
  
  res <- data.frame(term = NA, pVal = NA)
  
  while(length(left) > 0){
    
    ## Remove any remaining main effects
    if(all(!(is.na(res$pVal))) && max(res$pVal) < alpha)
    {
      break
    }
    
    res <- data.frame(term = NA, pVal = NA)
    
    for(effect in 1:length(left)){
      
      used <- c(inInteractions, interactions, left[-effect])
      cat(paste("Testing main effect: ", left[effect], "\n", sep=""))
      res[effect, 'term'] <- left[effect]
      
      fixedeffs <- paste(used,collapse="+")
      
      call <- paste(response, "~", fixedeffs, randEffect)
      
      if(fam == "poisson"){
        model2 <- glmer(formula = call, data = data, family = fam, 
                        control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
      }else{
        model2 <- lmer(formula = call, data = data,
                       control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
        
      }
      res[effect, 'pVal'] <- anova(refModel, model2)[2,'Pr(>Chisq)']
      print(res[effect, 'pVal'])
    }
    if( max(res$pVal, na.rm = TRUE) > alpha){ 
      toRemove <- res$term[which(res$pVal == max(res$pVal, na.rm = TRUE))]
      cat(paste("Removing", toRemove, "\n"))
      left <- left[-which(left == toRemove)]
    }
    
    used <- c(inInteractions, interactions, left)
    fixedeffs <- paste(used,collapse="+")
    call <- paste(response, "~", fixedeffs, randEffect)
    if(fam == "poisson"){
      refModel <- glmer(formula = call, data = data, family = fam, 
                        control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
    }else{
      refModel <- lmer(formula = call, data = data,
                       control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
      
    }
    
  }
  
  cat(paste("Returning final model\n"))
  return(refModel)
}




modelSimplificationAIC <- function(model = model, optimizer = "bobyqa", Iters = 2e5,  data){
  
  form <- as.character(model@call$formula)
  form <- form[-grep("~", form)]
  response <- form[-grep("\\+", form)]
  
  randEffect <- form[grep("\\+", form)]
  randEffect <- gsub("^.*\\(", "", randEffect) ## This also removes the first bracket
  randEffect <- paste0("+ (", randEffect)
  
  fam <- as.character(model@call$family)
  ## lmers don't have a family, so....
  if(length(fam) == 0) {fam <- "notpoisson"}
  
  
  all.terms <- rownames(anova(model))
  interactions <- all.terms[grep(":", all.terms)]
  main <- all.terms[!(all.terms %in% interactions)]
  
  res <- data.frame(term = NA, AIC = 0)
  
  
  if(fam == "poisson"){
    
    refModel <- glmer(formula = model@call$formula, data = data, family = fam, 
                      control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
    
  }else{
    refModel <- lmer(formula = model@call$formula, data = data, 
                     control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
  } ## How do I change this 
  
  AICRefModel <- 10000000 # an imaginary number, so that it doesn't fail on the first loop through. This gets replaced
  
    if(length(interactions) > 0){
      while(length(interactions) > 0){
        ## break loop if remaining interactions have AIC values lower than the ref model
        # (and none of the pvalues are NA)
        if(min(res$AIC, na.rm = TRUE) > AICRefModel)
        {
          break
        }
        
        
        res <- data.frame(term = NA, AIC = NA)
        
        for( inter in 1:length(interactions)){
          cat(paste("Testing interaction: ", interactions[inter], "\n", sep=""))
          used <- c(main, interactions[-inter])
          res[inter, 'term'] <- interactions[inter]
          
          fixedeffs <- paste(used,collapse="+")
          
          call <- paste(response, "~", fixedeffs, randEffect)
          
          if(fam == "poisson"){
            model2 <- glmer(formula = call, data = data, family = fam, 
                            control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
          }else{
            model2 <- lmer(formula = call, data = data,
                           control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
            
          }
          
          anov <- anova(refModel, model2)
          
          AnovaRow <- grep("model2", rownames(anov))
          res[inter, 'AIC'] <- anov[AnovaRow,'AIC']
          AICRefModel <- anov[grep("refModel", rownames(anov)),'AIC']
          
          print(AICRefModel - res[inter, 'AIC'])
          
          
          
        }
        
        if(min(res$AIC, na.rm = TRUE) < AICRefModel){
          toRemove <- res$term[which(res$AIC == min(res$AIC, na.rm = TRUE))]
          cat(paste("Removing", toRemove, "\n"))
          interactions <- interactions[-which(toRemove == interactions)]
        }
        
        # recreate a ref model
        used <- c(main, interactions)
        fixedeffs <- paste(used,collapse="+")
        call <- paste(response, "~", fixedeffs, randEffect)
        if(fam == "poisson"){
          refModel <- glmer(formula = call, data = data, family = fam, 
                            control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
        }else{
          refModel <- lmer(formula = call, data = data, 
                           control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
        }
      }
    }
    
    
    ## Deal with any main effects now
    
    ## TODO: What is no interactions to begin with
    
    
    ## What main effects are still in interactions
    matches <- unique (grep(paste(main,collapse="|"), 
                            interactions, value=TRUE))
    
    
    keep <- unique(unlist(strsplit(matches, "[:]")))
    
    left <- main[!(main %in% keep)]
    
    cat(paste("Testing removal of main effects\n"))
    
    inInteractions <- main[main %in% keep]
    
    res <- data.frame(term = NA, AIC = 0) ## unrealistic number, so doesn't fail first time
    
    while(length(left) > 0){
      
      ## Remove any remaining main effects
      if(min(res$AIC, na.rm = TRUE) > AICRefModel)
      {
        break
      }
      
      res <- data.frame(term = NA, AIC = NA)
      
      for(effect in 1:length(left)){
        
        used <- c(inInteractions, interactions, left[-effect])
        cat(paste("Testing main effect: ", left[effect], "\n", sep=""))
        res[effect, 'term'] <- left[effect]
        
        if(length(used) == 0){ used <- "1"
          cat("Reference model is intercept only model. Therefore, a negative AIC is good \n")}
        fixedeffs <- paste(used,collapse="+")
        
        call <- paste(response, "~", fixedeffs, randEffect)
        
        if(fam == "poisson"){
          model2 <- glmer(formula = call, data = data, family = fam, 
                          control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
        }else{
          model2 <- lmer(formula = call, data = data,
                         control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
          
        }
        
        anov <- anova(refModel, model2)
        AnovaRow <- grep("model2", rownames(anov))
        res[effect, 'AIC'] <- anov[AnovaRow,'AIC']
        AICRefModel <- anov[grep("refModel", rownames(anov)),'AIC']
        
        print(AICRefModel - res[effect, 'AIC'])
        
        
        
      }
      # print(paste('AICRefModel', AICRefModel))
      if( min(res$AIC, na.rm = TRUE) < AICRefModel){  
        toRemove <- res$term[which(res$AIC == min(res$AIC, na.rm = TRUE))]
        cat(paste("Removing", toRemove, "\n"))
        left <- left[-which(left == toRemove)]
      }
      
      used <- c(inInteractions, interactions, left)
      fixedeffs <- paste(used,collapse="+")
      call <- paste(response, "~", fixedeffs, randEffect)
      if(fam == "poisson"){
        refModel <- glmer(formula = call, data = data, family = fam, 
                          control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
      }else{
        refModel <- lmer(formula = call, data = data,
                         control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
        
      }
      
    }
  
  
  cat(paste("Returning final model\n"))
  return(refModel)
}
