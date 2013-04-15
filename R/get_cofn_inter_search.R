
###########################################################################################################
get.cofn.int.search<-function(model){
  # get coef for a model from fit object after running interval search 
  f1 <- model$cvreg
  res<-f1$glmnet.fit
  cof <- as.vector(coef(res, s=model$lambda))
  names.cof <- rownames(res$beta)
  cofn <- cof[which(cof != 0)]
  names(cofn) <- names.cof[which(cof != 0)] 
  bet <- res$beta[match(names(cofn), rownames(res$beta)),]
  
  return(cofn)
}  

