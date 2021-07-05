#Library
library(tidyr)
library(dplyr) 
library(diffusr) #RWR
library(caret) #CV
library(ggplot2)
library(hrbrthemes)
library(ranger) #RANGER
library(viridis)
library(viridisLite)
library(forcats)
library(pROC)

#Parameters
#' @param p number of genes.
#' @param n number of sample size.
#' @param y a vector with one category for each sample.
#' @param data dataset (nx(p+1)) of expression levels from RNA-seq.
#' @param weights a vector with modified selection probability for each p.
#' @param ntrees_g vector with number of trees to grow in each RF.
#' @param mtry_g vector with number of variables available for splitting at each tree node.

#######################################################
########Knowledge-Slanted Random Forest (RF)###########
#######################################################

slanted_rf=function(data,weights, ntrees_g, mtry_g){
  
  
  set.seed(23)
  
  myControl <- trainControl(
    method = "LOOCV", 
    verboseIter = TRUE, 
  )
  
  hyper_grid <- expand.grid(
    mtry       = mtry_g,
    min.node.size= 1,
    splitrule="gini"
  )
  
  grid_ntree <- expand.grid(
    num.tree     = ntrees_g
  )
  
  #results
  results_final_bias=list()
  
  #slanted-random forest
  
  for (i in 1:nrow(grid_ntree)){
    
    model_rf_bias <- train(factor(y) ~., 
                           data = data, 
                           method='ranger',
                           trControl = myControl,  
                           tuneGrid=hyper_grid, 
                           num.trees=grid_ntree$num.tree[i],
                           metric = "Accuracy",
                           split.select.weights=weights, #Pesos de Selección
    )
    
    results_final_bias[[i]]=model_rf_bias$results
  }
  
  
 return(data.frame(ntree=rep(ntrees_g,each=length(mtry_g)),
                                        do.call(rbind, results_final_bias))
  )
  

  
}


#######################################################
########## Conventional Random Forest (RF)#############
#######################################################

conv_rf=function(data, ntrees_g, mtry_g){
  
  
  set.seed(23)
  
  myControl <- trainControl(
    method = "LOOCV", 
    verboseIter = TRUE, 
  )
  
  hyper_grid <- expand.grid(
    mtry       = mtry_g,
    min.node.size= 1,
    splitrule="gini"
  )
  
  grid_ntree <- expand.grid(
    num.tree     = ntrees_g
  )
  
  #results
  results_final_bias=list()
  
  #slanted-random forest
  
  for (i in 1:nrow(grid_ntree)){
    
    model_rf_bias <- train(factor(y) ~., 
                           data = data, 
                           method='ranger',
                           trControl = myControl,  
                           tuneGrid=hyper_grid, 
                           num.trees=grid_ntree$num.tree[i],
                           metric = "Accuracy",
                           #split.select.weights=weights,
    )
    
    results_final_bias[[i]]=model_rf_bias$results
  }
  
  
  return(data.frame(ntree=rep(ntrees_g,each=length(mtry_g)),
                    do.call(rbind, results_final_bias))
  )
  
  
  
}

####Example
#y=rep(1:0, each=5)
#x=matrix(rnorm(10*1000,mean=3,sd=3.0), 10,1000) 
#data=data.frame(y,x)
#ntrees_g=c(10,50,100,200,300,400,500,1000)
#mtry_g= c(10, 50, 125,200,300,400,500,600,700,800, 1000) 
#weights=runif(1000,0,1)
#slanted=slanted_rf(data,weights, ntrees_g, mtry_g)
#conventional=conv_rf(data, ntrees_g, mtry_g)


#########Graphs

# #methods=rep(c("slanted","RF")
# global_results= rbind(slanted, conventional)
# global_results=cbind(methods=rep(c("slanted","RF"), each=length(mtry_g)*length(ntrees_g)),
#                  global_results)
# 
# global_results %>% 
#   ggplot(aes(x=factor(ntree), y=Accuracy, fill=methods)) +
#   geom_boxplot() + 
#   theme(legend.title = element_blank(), 
#         legend.position = c(0.9, 0.2)) +
#   scale_x_discrete(name ="ntree") +
#   scale_y_continuous( breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1)) 
# 
# 
# global_results %>%
#   ggplot(aes(x=factor(mtry), y=Accuracy, fill=methods)) +
#   geom_boxplot() + 
#   theme(legend.title = element_blank(), 
#         legend.position = c(0.9, 0.2)) +
#   scale_x_discrete(name ="mtry") +
#   scale_y_continuous( breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1)) 
# 


