
HypothesisTesting <- function(model = model, data, TestingGroups){
  
  optimizer <- "bobyqa"
  Iters <- 2e5
  
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
  
  res <- data.frame(term = NA, AIC = 0)
  
  
  if(fam == "poisson"){
    
    refModel <- glmer(formula = model@call$formula, data = data, family = fam, 
                      control = glmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters))) 
    
  }else{
    refModel <- lmer(formula = model@call$formula, data = data, 
                     control = lmerControl(optimizer = optimizer,optCtrl=list(maxfun=Iters)))
  } ## How do I change this 
  
  
  # AICRefModel <- 10000000 # an imaginary number, so that it doesn't fail on the first loop through. This gets replaced
  
  for(group in 1:length(TestingGroups)){
    cat(paste("Testing group: ", TestingGroups[group], "\n", sep = ""))
    
    res[group, 'term'] <- TestingGroups[group]
    
    used <- all.terms[-(all.terms %in% TestingGroups[group])]
    
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
    res[group, 'AIC'] <- anov[AnovaRow,'AIC']
    AICRefModel <- anov[grep("refModel", rownames(anov)),'AIC']
    
    # print(AICRefModel - res[inter, 'AIC'])
    
  }
  res[TestingGroups + 1, 'term'] <- "RefModel"
  res[TestingGroups + 1, 'term'] <- AICRefModel
  
  return(res)
}
        
        
        
        