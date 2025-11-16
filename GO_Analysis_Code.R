#Gene ontolgy analysis based on differntial gene expression 

gn<-read.delim('4 vs. 2.txt')

write.table(gn[,1],'gn.txt',row.names = FALSE,col.names = FALSE)

q <-read.table("fastB")
head(q)

library(AnnotationDbi)
library(org.At.tair.db)
library(BiocManager)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(BiocGenerics)
library(KEGGREST)
library(biomartr)
library(tidyverse)
library(enrichplot)
library(ggridges)
library(ggrepel)



columns(org.At.tair.db)

h<- read.csv('gene_name.csv',row.names = 1)
head(h)


q$gene_name <- mapIds(org.At.tair.db,
                      keys=q$V2,
                      keytype="TAIR",
                      column="GENENAME",
                      multiVals="first")

q$symbol <- mapIds(org.At.tair.db,
                   keys=q$V2,
                   keytype="TAIR",
                   column="SYMBOL",
                   multiVals="first")
q$go <- mapIds(org.At.tair.db,
               keys=q$V2,
               keytype="TAIR",
               column="GO",
               multiVals="first")
q$go <- mapIds(h,
               keys=q$V1,
               keytype="prelim",
               column="Feature.ID",
               multiVals="first")

j<- q[c(1,2,13,14)]
head(j)

write.csv(j, "allgenes.csv",row.names = FALSE)


#GO analysis


GO<- function(x,y,z) {
  rr <- read.csv(paste('allgenes',z,'.csv',sep=''))
  f<- rr %>%
    select('arabidopsis_homolog', paste('Logfc',x,sep=''),paste('pvalue',x,sep=''))
  
  colnames(f)[3] <- 'pvalue'
  colnames(f)[2] <- 'logfc'
  
  f<- f %>% filter(pvalue <0.05)
  
  fu<- f %>% filter(logfc >0.5)
  fd<- f %>% filter(logfc < (-0.5))
  
  genes_to_test <- fu$arabidopsis_homolog
  
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'MF')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'U_MF.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'UP_MF.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #BP 
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'BP')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'U_BP.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'UP_BP.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #CC
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'CC')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'U_CC.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'UP_CC.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #ALL
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'ALL')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'U_ALL.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'U_ALL.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #down genes
  #MF
  genes_to_test = fd$arabidopsis_homolog
  
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'MF')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'D_MF.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'D_MF.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #BP 
  
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'BP')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'D_BP.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'D_BP.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #CC
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'CC')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'D_CC.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'D_CC.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
  #ALL
  GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.At.tair.db", keyType = "TAIR", ont = 'ALL')
  
  j <- as.data.frame(GO_results)
  
  write.csv(j, file = paste(y,'D_ALL.csv',sep=''), row.names = FALSE)
  
  fit <- plot(barplot(GO_results, showCategory = 15))
  
  png(paste(y,'D_ALL.png',sep=''), res = 250, width = 1400, height = 1800)
  print(fit)
  dev.off()
  
}

GO('42','go_roots_42_','R')
GO('61','go_roots_61_','R')
GO('41', 'go_roots_41_','R')
GO('21', 'go_roots_21_','R')

GO('42', 'go_L_42_',"L")
GO('61','go_L_61_','L')
GO('41', 'go_L_41_','L')
GO('21', 'go_L_21_','L')






