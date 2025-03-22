
library(tidyverse)
library(gprofiler2)
library(msigdbr)
library("org.Hs.eg.db")
require(DOSE)
library(DESeq2)
library(survival)
#---------------------among new MM----------------------------
df<-read.csv("./data/df_allrace_new.csv", row.names=1, check.names=FALSE) 
#get coldata 
df_raw<- df[, 89:length(df)]
coldata<- df[, !colnames(df) %in% colnames(df_raw)]
#normalization
dds<- DESeqDataSetFromMatrix(countData = t(df_raw), colData = coldata, design = ~ 1)
vsd<- vst(dds, blind = FALSE)
count_norm<- assay(vsd)
count_norm<- t(count_norm)
df_norm<- cbind(coldata, count_norm)

#-------------------with candidate genes---------------------
#load data 
#selected genes names 
load("./data/candi_genes.RData")

#rank gene list by log2foldchange
coef<- unlist(res_filtered$log2FoldChange, recursive = T, use.names = F)
names(coef)<- rownames(res_filtered)
genelist<- sort(coef, decreasing=TRUE)
gene_names<- names(genelist)
write.table(gene_names, file = "./output/candi_genelist.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

#---------g:profiler------------------------------------------
glfe<- gost(query = names(genelist),
            organism = "hsapiens", 
            ordered_query = T,
            correction_method = "fdr",
            user_threshold = 0.05,
            sources = c("GO", "KEGG"),
            significant = T, 
            domain_scope = "annotated", 
            evcodes = T
) 

#filter to generate gmt for enrichment map 
filtered_path<- glfe$result %>% 
  filter(term_size>=5 & term_size<=350, intersection_size>=3) 

gem <- filtered_path[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
head(gem, 3)
#save 
write.table(gem, file = "./output/candi_gene_gem.txt", sep = "\t", quote = F, row.names = F)

#pick the top 5 
highlight_terms<- filtered_path %>% 
  arrange(source, p_value) %>% 
  group_by(source) %>% 
  dplyr::slice(1:5)

plot<- gostplot(glfe, interactive = F)+ scale_size_continuous(range=c(1,1.5))
pdf("./plots/gprofiler.pdf", h= 8, w=12)
publish_gostplot(plot, highlight_terms =highlight_terms$term_id)
dev.off()

#-------------------with final genes---------------------
#load data 
#selected genes names 
final_genes<- read.csv("./output/final_gene.csv") %>% 
  filter(!Gene %in% c("age_diag_calculated", "sex_emr", "(Intercept)")) %>% 
  column_to_rownames(var="Gene")

#rank gene list by coefs
coef<- unlist(final_genes$Coefficient, recursive = T, use.names = F)
names(coef)<- rownames(final_genes)
genelist<- sort(coef, decreasing=TRUE)
gene_names<- names(genelist)
write.table(gene_names, file = "./output/final_genelist.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\n")

#---------g:profiler------------------------------------------
glfe<- gost(query = names(genelist),
            organism = "hsapiens", 
            ordered_query = T,
            correction_method = "fdr",
            user_threshold = 0.05,
            sources = c("GO", "KEGG"),
            significant = T, 
            domain_scope = "annotated", 
            evcodes = T
) 

#filter to generate gmt for enrichment map 
filtered_path<- glfe$result %>% 
  filter(term_size>=5 & term_size<=350, intersection_size>=3) 

gem <- filtered_path[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
head(gem, 3)
#save 
write.table(gem, file = "./output/final_gene_gem.txt", sep = "\t", quote = F, row.names = F)

#pick the top 5 
highlight_terms<- filtered_path %>% 
  arrange(source, p_value) %>% 
  group_by(source) %>% 
  dplyr::slice(1:5)

plot<- gostplot(glfe, interactive = F)+ scale_size_continuous(range=c(1,1.5))
pdf("./plots/gprofiler_finalgene.pdf", h= 8, w=12)
publish_gostplot(plot, highlight_terms =highlight_terms$term_id)
dev.off()
