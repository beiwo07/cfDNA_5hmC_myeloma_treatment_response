
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)


#load deseq results done within training set 
load("./data/res.RData")

pdf("./plots/volcano.pdf", h= 8, w=12)
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-0.3,0.3),
                ylim= c(0,5),
                pCutoff = 0.05,
                FCcutoff = 0.1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 3.0,
                labSize = 5.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','Log (base 2) FC>0.1','p-value<0.05',
                               'p-value<0.05 & Log (base 2) FC>0.1'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0)
dev.off()

