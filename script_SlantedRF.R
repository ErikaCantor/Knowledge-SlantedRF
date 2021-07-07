#######################################################################
### All rights are reserved by the authors.
### Title: Knowledge-Slanted Random Forest
### Author: Erika Cantor
### Created on: 2021/07/21
#######################################################################

#Functions
source("RWR_PriorKnowledge.R")
source("RF_methods.R")
source("Predictions_validation.R")
source("TopGenes.R")

#Parameters and Data
r=0.3 #back probability in RWR
seednodes=read.table("seed", sep = "\t") #a vector with prior known genes.
scores= read.table("scores", header = TRUE, sep = ",",row.names="col_0")  #scores of protein-protein interactions (PPI) network 
p=ncol(scores) #number of genes
data=read.table("data", header = TRUE, sep="\t") 
names(data)[1] <- "y"
n=nrow(data)

#'Parameters for RF
#'ntrees_g: grid of trees
#'mtry_g:  grid of mtry
#'ntree: optimal ntree
#'#mtry: optimal mtry

#weights for knowledge-slanted random forest
rwr_results=rwr_prior(p,scores,r,seednodes$V1)
weights=rwr_results$weights

#Knowledge_Slanted Random Forest
slanted=slanted_rf(data,weights, ntrees_g, mtry_g)

#Conventional Random Forest
conventional=conv_rf(data, ntrees_g, mtry_g)

#Performance visualization to choose optimal parameters.

methods=rep(c("slanted","RF"), each=length(mtry_g)*length(ntrees_g))
global_results= rbind(slanted, conventional)
global_results=cbind(methods=rep(c("slanted","RF"), each=length(mtry_g)*length(ntrees_g)),
                  global_results)
global_results %>% 
   ggplot(aes(x=factor(ntree), y=Accuracy, fill=methods)) +
   geom_boxplot() + 
   theme(legend.title = element_blank(), 
         legend.position = c(0.9, 0.2)) +
   scale_x_discrete(name ="ntree") +
  scale_y_continuous( breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1)) 
 
 
global_results %>%
   ggplot(aes(x=factor(mtry), y=Accuracy, fill=methods)) +
   geom_boxplot() + 
   theme(legend.title = element_blank(), 
         legend.position = c(0.9, 0.2)) +
   scale_x_discrete(name ="mtry") +
   scale_y_continuous( breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1)) 


#Predictions 
predicts=predictions_RF(n,weights,ntree, mtry, data)

#predictions for slanted-RF
predict_slanted=unlist(lapply(predicts$slanted_pre, which.max))

#predictions for conventional-RF
predict_conventional=unlist(lapply(predicts$conventional_pre, which.max))

#Confusion Matrix

confusionMatrix=list(slantedRF=table(prediction=predict_slanted, observed=data$y), 
                     conventionalRF=table(prediction=predict_conventional, observed=data$y))

#Top20genes Slanted Random Forest

model=ranger(factor(y) ~ ., 
                   num.trees=ntree, 
                   mtry=mtry, 
                   dependent.variable.name=y, 
                   data = data, 
                   split.select.weights=weights,
                   verbose=TRUE,
                   seed=23, 
                   classification=TRUE,
                   probability=TRUE,
)

topSlanted_RF=top_variables(model,ntree)



