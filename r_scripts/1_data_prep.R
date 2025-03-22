
library(dplyr)
library(tidyverse)
library(glmnet)
library(DESeq2)

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#read raw count data 
df<-read.csv("./data/df_allrace_new.csv", row.names=1, check.names=FALSE)
df<- df[!is.na(df$response_abstracted_cr), ]
df<- df %>% 
  mutate(kappa_lambda_clean_cat= relevel(factor(kappa_lambda_clean_cat), ref = "normal"), 
         response_abstracted_cr = case_when(response_abstracted_cr == "<CR" ~ 0,
                                            response_abstracted_cr == ">=CR" ~ 1,
                                            TRUE ~ NA_real_),
         response_abstracted_cr= factor(response_abstracted_cr, levels= c(0,1)), 
         iss_bi= ifelse(iss_derived_m %in% c("2", "3"), "2/3", iss_derived_m)
         )
#format clinical data and raw count 
col_order<- c(colnames(df[1:88]),"iss_bi")
coldata<- df[, colnames(df) %in% col_order] 
count<- df[, !colnames(df) %in% colnames(coldata)]
count<- t(count)
#--------------------------------------------------------------------------------------------------------------------------------------------------#
## split training and testing sets and normalization 
set.seed(18)
#stratified sampling
coldata<- coldata %>% rownames_to_column(var="errcid") 
train_col<- coldata %>% 
  group_by(response_abstracted_cr, age_diag_cat, sex_emr, race_composite, iss_derived_m) %>% 
  sample_frac(size=0.6, replace = F) %>% 
  ungroup() %>% 
  column_to_rownames(var="errcid")
coldata<- coldata %>% column_to_rownames(var="errcid")
test_col<- coldata[!rownames(coldata)%in%rownames(train_col), ]

count_train<- count[ ,colnames(count) %in% rownames(train_col)]
count_test<- count[ ,colnames(count) %in% rownames(test_col)]

## normalization after the split 
#train set
#match dfs
order<- match(rownames(train_col), colnames(count_train))
count_train<- count_train[,order]
#normalization
dds<- DESeqDataSetFromMatrix(countData = count_train, colData = train_col, design = ~ 1)
vsd<- vst(dds, blind = FALSE)
norm_train<- assay(vsd)
norm_train<- t(norm_train)
df_train<- cbind(train_col, norm_train)
#testFea
#match dfs
order<- match(rownames(test_col), colnames(count_test))
count_test<- count_test[,order]
#normalization
dds<- DESeqDataSetFromMatrix(countData = count_test, colData = test_col, design = ~ 1)
vsd<- vst(dds, blind = FALSE)
norm_test<- assay(vsd)
norm_test<- t(norm_test)
df_test<- cbind(test_col, norm_test)

## save split ids
df_train$group<- "train"
df_test$group<- "test"
split<- rbind(df_train, df_test) %>% 
  dplyr::select(group) %>% 
  rownames_to_column(var="errcid") 
split<- as.data.frame(split)
df_train<- df_train %>% 
  dplyr::select(-group)
df_test<- df_test %>% 
  dplyr::select(-group)

#--------------------------------------------------------------------------------------------------------------------------------------------------#
#DE analysis in train set to select candidate genes with p<0.05
dds<- DESeqDataSetFromMatrix(countData = count_train, 
                             colData = train_col, 
                             design = ~ age_diag_cat+ response_abstracted_cr)
dds<- DESeq(dds)
res<- results(dds, 
              contrast = c("response_abstracted_cr",0, 1), 
              pAdjustMethod= "BH")

#save a res for plotting
save(res, file= "./data/res.RData")

#keep p<0.1 
res_filtered<- res[res$pvalue<0.05 & abs(res$log2FoldChange)>0.1,]
#--------------------------------------------------------------------------------------------------------------------------------------------------#
save(df_train, file="./data/df_train.RData")
save(df_test, file= "./data/df_test.RData")
save(res_filtered, file= "./data/candi_genes.RData")