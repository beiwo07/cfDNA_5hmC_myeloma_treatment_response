
##################################

#the script is meant to run on hpc 

##################################

library(dplyr)
library(tidyverse)
#library(mlr3proba)
library(mlr3tuning)
library(mlr3learners)
library(glmnet)
library(DESeq2)
library(conflicted)

sessionInfo()

args = commandArgs(trailingOnly=TRUE)
sbatch_index = args[1]

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#for test locally
#load("./data/candi_genes.RData")
#load("./data/df_train.RData")
#load("./data/df_test.RData")

load("/gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/data/candi_genes.RData")
load("/gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/data/df_train.RData")
load("/gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/data/df_test.RData")
#-----------------------------100 iterations of the following---------------------------------------------------------------------------------------------------------------------#
#2. Feature selection using el  
#a. Tune alph with the best c-index (training)
candi_gene<- df_train[, c(rownames(res_filtered), "response_abstracted_cr")]
task<- as_task_classif(candi_gene, target = "response_abstracted_cr")
measure<-  msr("classif.auc")
#cv<- rsmp("cv", folds=3)
cv<- rsmp("cv", folds=5)
#cv<- rsmp("loo")
#alpha_vals<- seq(0,1,by=0.2)
alpha_vals<- seq(0,1,by=0.05)
search_space<- tnr("grid_search")
penalty.factor<- unlist(lapply(res_filtered$log2FoldChange, function(x) 1/abs(x)))
glmnet<- lrn("classif.cv_glmnet", 
             alpha= to_tune(c(alpha_vals)), 
             #alignment= "lambda",
             penalty.factor=penalty.factor, 
             predict_type="prob")
#find best alpha
set.seed(sbatch_index)
instance_el<- mlr3tuning::tune(
  tuner = search_space, 
  task = task, 
  learner = glmnet, 
  resampling = cv, 
  measure= measure
)

#b. Use best alpha to refit model to find best lambda (training)
Y<- data.matrix(candi_gene$response_abstracted_cr)
X<- data.matrix(candi_gene[, colnames(candi_gene)!= "response_abstracted_cr"])
cv_mod<- cv.glmnet(x= X, y= Y, family="binomial", alpha=instance_el$result$alpha, 
                   penalty.factor=penalty.factor)
best_lambda<- cv_mod$lambda.min

#c. Additional feature selection: use best alpha and lambda to fit the model, remove genes with coef less of a threshold, if any, and generate list of final genes (training)
mod<- glmnet(x=X, y= Y, family = "binomial", alpha = instance_el$result$alpha, lambda = best_lambda, 
             penalty.factor=penalty.factor)
final_gene_ls<- as.data.frame(as.matrix(coef(mod))) %>% 
  dplyr::slice(-1) %>% #remove intercept 
  dplyr::filter(abs(s0)>=0.1) %>% #restrict the absolute value of coef to be greater than 0.2
  rownames_to_column(var="gene")
#--------------------------------------------------------------------------------------------------------------------------------------------------#
#save 
write.csv(final_gene_ls, file= paste0("/gpfs/data/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM/res/final_genes_", sbatch_index, ".csv"))