#######################################################################
### All rights are reserved by the authors.
### Title: Knowledge-Slanted Random Forest
### Author: Erika Cantor
### Created on: 2021/06/21
#######################################################################

#'Packages
library(diffusr) 
library(caret)
library(ggplot2)
library(ranger) 

#######################################################
#Step 1: Random Walk with restart using prior knowledge
#######################################################

#' @param p number of genes or nodes.
#' @param scores matrix (pxp) with scores from string for each gene.
#' @param r back-probability.
#' @param p0 the initial probability of being at node i.
#' @param seednodes a vector with prior known genes.
#' @return pt$p.inf the probability of each node being a candidate-gene.


#Set initial probability in all nodes
p0 <- rep(0.00001, p) 

#Set initial probability equal to 1 for seed nodes
p0[which(rownames(scores)%in%seednodes)]<-1 

#adjacency matrix 
graph <- as.matrix(scores) 

#computation of stationary distribution
set.seed(23)
pt  <- random.walk(p0, graph,r=0.3, niter = 10000, thresh = 1e-06)

#Final probability of each node being a candidate-gene.
pt$p.inf

#######################################################
#Step 2: Knowledge-Slanted Random Forest (RF)
#######################################################

#' @param p number of genes.
#' @param n number of sample size.
#' @param y a vector with one category for each sample.
#' @param data dataset (nx(p+1)) of expression levels from RNA-seq.
#' @param split.select.weights a vector with modified selection probability for each p.
#' @param num.trees number of trees to grow in each RF.
#' @param mtry number of variables available for splitting at each tree node.


slanted_RF <- ranger(y ~ ., 
       num.trees=500, 
       mtry=500, 
       dependent.variable.name=y, 
       data = data,  
       split.select.weights=pt$p.inf, 
       seed=23,
       classification=TRUE,
       probability=TRUE,
       )
