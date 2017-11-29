library(lme4)

modelSimplification <- function(model = model, data, alpha = 0.05){


#cars$type <- as.factor(rep(letters[1:4], length = 50))
#cars$rand <- as.factor(rep(letters[5:7], length = 50))

#model <- lmer(speed ~ type * dist + (1|rand), data = cars)


  form <- as.character(model@call$formula)
  form <- form[-grep("~", form)]
  response <- form[-grep("\\+", form)]
  
  randEffect <- form[grep("\\+", form)]
  randEffect <- gsub("^.*\\(", "", randEffect) ## This also removes the first bracket
  randEffect <- paste0("+ (", randEffect)
  
  all.terms <- rownames(anova(model))
  interactions <- all.terms[grep(":", all.terms)]
  main <- all.terms[!(all.terms %in% interactions)]
  
  res <- data.frame(term = NA, pVal = NA)
  
  
  refModel <- model ## How do I change this 
  
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
        
          used <- c(main, interactions[-inter])
          res[inter, 'term'] <- interactions[inter]
        
          fixedeffs <- paste(used,collapse="+")
        
          call <- paste(response, "~", used, randEffect)
        
          model2 <- lmer(formula = call, data = data) ## What about any other model specifics
        
          res[inter, 'pVal'] <- anova(refModel, model2)[2,'Pr(>Chisq)']
        }
      
        if( max(res$pVal, na.rm = TRUE) > alpha){ 
          toRemove <- res$term[which(res$pVal == max(res$pVal, na.rm = TRUE))]
        
          interactions <- interactions[-which(toRemove == interactions)]
        }
      # recreate a ref model
      used <- c(main, interactions)
      fixedeffs <- paste(used,collapse="+")
      call <- paste(response, "~", used, randEffect)
      refModel <- lmer(formula = call, data = data) ## What about any other model specifics
      }
    }
  
  
  ## Deal with any main effects now
  
  ## TODO: What is no interactions to begin with
  
  
  ## What main effects are still in interactions
    matches <- unique (grep(paste(main,collapse="|"), 
                          interactions, value=TRUE))
  
  
    keep <- unique(unlist(strsplit(matches, "[:]")))
  
    left <- main[!(main %in% keep)]
    inInteractions <- main[main %in% keep]
  
    res <- data.frame(term = NA, pVal = NA)
  
    while(length(left) > 0){
    
      ## Remove any remaining main effects
      if(all(!(is.na(res$pVal))) && max(res$pVal) < alpha)
      {
        break
      }
    
      for(effect in 1:length(left)){
      
        used <- c(inInteraction, interactions, left[-effect])
      
        res[inter, 'term'] <- left[effect]
      
        fixedeffs <- paste(used,collapse="+")
      
        call <- paste(response, "~", used, randEffect)
      
        model2 <- lmer(formula = call, data = cars) ## What about any other model specifics
      
        res[inter, 'pVal'] <- anova(refModel, model2)[2,'Pr(>Chisq)']
      
      }
      if( max(res$pVal, na.rm = TRUE) > alpha){ 
        toRemove <- res$term[which(res$pVal == max(res$pVal, na.rm = TRUE))]
      
        left <- left[-which(toRemove == left)]
      }
    
      used <- c(inInteraction, interactions, left)
      fixedeffs <- paste(used,collapse="+")
      call <- paste(response, "~", used, randEffect)
      refModel <- lmer(formula = call, data = cars)
    
    }
  
  return(refModel)
}
