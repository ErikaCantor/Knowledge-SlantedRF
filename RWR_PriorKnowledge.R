##########################################################
################ Random Walk with restart ################
##########################################################

#'Packages
library(diffusr) 

#' @param p number of genes or nodes.
#' @param scores matrix (pxp) with scores from string for each gene.
#' @param r back-probability.
#' @param p0 the initial probability of being at node i.
#' @param seednodes a vector with prior known genes.
#' @return pt$p.inf the probability of each node being a candidate-gene.

rwr.prior<-function(p,scores,r,seednodes,seed1){
  
  #Set initial probability in all nodes
  p0 <- rep(0.00001, p) 
  
  #Set initial probability equal to 1 for seed nodes
  p0[which(rownames(scores)%in%seednodes)]<-1 
  
  #adjacency matrix 
  graph <- as.matrix(scores) 
  
  #computation of stationary distribution
  set.seed(seed1)
  pt  <- random.walk(p0, graph,r=r, niter = 10000, thresh = 1e-06)
  
  #Probability of each node is a candidate gene.
  return(pt$p.inf)
  
  
}

