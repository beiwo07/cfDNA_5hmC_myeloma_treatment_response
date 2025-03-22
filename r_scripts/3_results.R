
library(dplyr)
library(glmnet)
library(tidyverse)
library(broom)

#read results of 100 iterations from hpc 
setwd("/Volumes/chiu-lab/bw_folder/dissertation/5hmc_tx/interation_100_allrace_newMM")
path<- "./res"

files<- list.files(path, pattern = "\\.csv$", full.names = T)

genes<- NULL

for (file in files){
  df<- read.csv(file)
  genes<- append(genes, df$gene)
}

gene_counts<- as.data.frame(table(genes)) %>% 
  arrange(Freq) %>% 
  dplyr::filter(Freq==100)

# use this genes in testing set
load("./data/df_train.RData")
load("./data/df_test.RData")

#1. Use final genes to fit in regular el model, to get age- and sex-adjusted coefs to calculate wp_scores
X<- df_train %>% 
  dplyr::select(c(gene_counts$genes, "age_diag_calculated", "sex_emr")) %>% 
  data.matrix()
Y<- data.matrix(df_train$response_abstracted_cr)

#find best lambda
cv_mod<- cv.glmnet(x= X, y= Y, family="binomial", alpha=0)
best_lambda<- cv_mod$lambda.min
mod<- glmnet(x=X, y= Y, family = "binomial", alpha = 0,  lambda = best_lambda)

#get weight from the 70% sample 
wt<- as.list(coef(mod))
#remove the coef for age and sex
wt<- wt[2:(length(wt)-2)]

#get scores 
final_gene_train<- df_train %>% 
  dplyr::select(gene_counts$genes)
final_gene_test<- df_test %>% 
  dplyr::select(gene_counts$genes)

score_train<- as.data.frame(final_gene_train*wt) %>% 
  mutate(wp_score= rowSums(across(everything()))) %>% 
  dplyr::select(wp_score) %>% 
  rownames_to_column(var= "errcid") %>% 
  left_join(df_train %>% 
              rownames_to_column(var="errcid"), by= "errcid") 

score_test<- as.data.frame(final_gene_test*wt) %>% 
  mutate(wp_score=rowSums(across(everything()))) %>% 
  dplyr::select(wp_score) %>% 
  rownames_to_column(var= "errcid") %>% 
  left_join(df_test %>% 
              rownames_to_column(var="errcid"), by= "errcid")  

#2. bootstrapping to calculate CIs for the coefs in the last model 
set.seed(1)
n_bootstraps <- 1000  # Number of bootstrap iterations
n <- nrow(X)  # Number of rows in the dataset
boot_mod_sample <- glmnet(x = X, y = Y, family = "binomial", alpha = 0, lambda = best_lambda)
n_coefs<- length(coef(boot_mod_sample))
boot_coefs <- matrix(NA, ncol = n_coefs, nrow = n_bootstraps)

# Bootstrapping loop
for (i in 1:n_bootstraps) {
  # Resample the data
  indices <- sample(1:n, size = n, replace = TRUE)
  X_boot <- X[indices, ]
  Y_boot <- Y[indices, ]
  
  # Fit the elastic net model to the bootstrapped sample
  boot_mod <- glmnet(x = X_boot, y = Y_boot, family = "binomial", alpha = 0, lambda = best_lambda)
  
  # Store the coefficients for the final genes
  boot_coefs[i, ] <- as.vector(coef(boot_mod))
}

# Calculate 95% confidence intervals
lower_ci <- apply(boot_coefs, 2, quantile, 0.025)
upper_ci <- apply(boot_coefs, 2, quantile, 0.975)

# Combine the coefficients and confidence intervals into a dataframe
final_gene_names <- rownames(coef(mod))  # Get gene names
ci_df <- data.frame(
  Gene = final_gene_names,
  Coefficient = as.vector(coef(mod)),
  Lower_CI = lower_ci,
  Upper_CI = upper_ci
)

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#assess HR 
#in training set 
mod_unadj_train<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated, 
                      data= score_train, family = "binomial")
mod_adj_train<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated +sex_emr+iss_bi, 
                    data= score_train, family = binomial("logit"))

#in testing set  
mod_unadj_test<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated, 
                     data= score_test, family = binomial("logit"))
mod_adj_test<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated + sex_emr+iss_bi, 
                   data= score_test, family = binomial("logit"))

wp_score_p<- data.frame(
  train_unadj_coef= summary(mod_unadj_train)$coefficients[2],
  train_unadj_p= summary(mod_unadj_train)$coefficients[11],
  train_adj_coef= summary(mod_adj_train)$coefficients[2],
  train_adj_p= summary(mod_adj_train)$coefficients[23],
  test_unadj_coef= summary(mod_unadj_test)$coefficients[2],
  test_unadj_p= summary(mod_unadj_test)$coefficients[11],
  test_adj_coef= summary(mod_adj_test)$coefficients[2],
  test_adj_p= summary(mod_adj_test)$coefficients[23]
)

tidy(mod_unadj_train, exponentiate= T, conf.int = T)
tidy(mod_adj_train, exponentiate= T, conf.int = T)

tidy(mod_unadj_test, exponentiate= T, conf.int = T)
tidy(mod_adj_test, exponentiate= T, conf.int = T)

#adjusting for different variables 
df_response<- readRDS("Y:/bw_folder/dissertation/5hmc_tx/analysis_drug/treatment_labeling/df_response.RData")
#only examine in testing set for now
score_test<- score_test %>% 
  mutate(mrn= as.character(mrn)) %>% 
  left_join(df_response %>% 
              select(mrn, dtq), by="mrn") %>% 
  mutate(dtq_bi= ifelse(is.na(dtq)|dtq%in% c("quad","doub"), 0, 1))
#kappa-lambda 
mod<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated + sex_emr+ kappa_lambda_clean_cat, 
          data= score_test, family = binomial("logit"))
tidy(mod, exponentiate= T, conf.int = T)

#treatment type
mod<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated + sex_emr+ dtq_bi, 
          data= score_test, family = binomial("logit"))
tidy(mod, exponentiate= T, conf.int = T)

#all
mod<- glm(response_abstracted_cr~ wp_score+ age_diag_calculated + sex_emr+ dtq_bi+ ldh_bi+ iss_bi, 
          data= score_test, family = binomial("logit"))
tidy(mod, exponentiate= T, conf.int = T)

#--------------------------------------------------------------------------------------------------------------------------------------------------#
write.csv(ci_df, file= "./output/final_gene.csv")
write.csv(score_train, file= "./output/score_train.csv")
write.csv(score_test, file= "./output/score_test.csv")