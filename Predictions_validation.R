##############################################################
############Optimal Performance  and predictions##############
############           using LOOCV              ##############
##############################################################


#Parameters for optimal RF
#' @param p number of genes.
#' @param n number of sample size.
#' @param y a vector with one category for each sample.
#' @param data dataset (nx(p+1)) of expression levels from RNA-seq.
#' @param weights a vector with modified selection probability for each p.
#' @param ntree number of trees to grow in each RF.
#' @param mtry  number of variables available for splitting at each tree node.


predictions_RF=function(n, weights, ntree, mtry, data){
  y=data[,1]
  set.seed(23) 
  predictions_slanted=list()
  predictions_conv=list()
  for(i in 1:n){
    
    slanted_randomF=ranger(factor(y) ~ ., 
                        num.trees=ntree, 
                        mtry=mtry, 
                        dependent.variable.name=y, 
                        data = data[-i,],  
                        split.select.weights=weights,
                        verbose=TRUE,
                        seed=23, 
                        classification=TRUE,
                        probability=TRUE,
                        keep.inbag=TRUE
                        
    )
    
    
    predictions_slanted[[i]] <- predict(slanted_randomF, data=data[i,], type="response")$predictions 
    
    conv_randomF=ranger(factor(y) ~ ., 
                        num.trees=ntree, 
                        mtry=mtry, 
                        dependent.variable.name=y, 
                        data = data[-i,],  
                        #split.select.weights=weights,
                        verbose=TRUE,
                        seed=23, 
                        classification=TRUE,
                        probability=TRUE,
                        keep.inbag=TRUE
                        
    )
    
    predictions_conv[[i]] <- predict(conv_randomF, data=data[i,], type="response")$predictions 
    
    
  }
#predictions vector for each method
  return(list(slanted_pre=predictions_slanted,conventional_pre=predictions_conv))
  
}


#Examples
#y=rep(1:0, each=5)
#x=matrix(rnorm(10*1000,mean=3,sd=3.0), 10,1000) 
#data=data.frame(y,x)
#ntree=500
#mtry= 50
#n=nrow(data)
#weights=runif(1000,0,1)
#slanted=slanted_rf(data,weights, ntrees_g, mtry_g)
#conventional=conv_rf(data, ntrees_g, mtry_g)
#predicts=predictions_RF(n,weights,ntree, mtry, data)




