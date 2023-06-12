#BiocManager::install("limma")

setwd("/home/keerthana2/AML/ILLUMINA/")
library(limma)
library(GEOquery)
library(dplyr)
library(tidyverse)

#untar the compressed files
unzip<- read.ilmn("/home/keerthana2/AML/ILLUMINA/GSE30029/GSE30029_non-normalized.txt.gz",probeid = "ID_REF", expr ="48" )

#background correction & normalization
exp<- backgroundCorrect(unzip, method = 'normexp')
exp<-normalizeBetweenArrays(exp)

#extracting expression values
exp_values<-as.data.frame(exp[[2]]) 

#assigning var to colname
var<-colnames(exp_values)

#renaming columns
colnames(exp_values)<- paste0('48',var)

#seriesmatrix download
gse <- getGEO("GSE30029", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE30029_series_matrix.txt.gz@featureData@data

# subset
###############################
feature..data <- feature.data[,c(1,11,14,28)]
feature..data$Synonyms<- ifelse(feature..data$Synonyms=="", NA, feature..data$Synonyms)
feature..data$Entrez_Gene_ID<- ifelse(feature..data$Entrez_Gene_ID=="", NA, feature..data$Entrez_Gene_ID)
feature..data$Symbol<- ifelse(feature..data$Symbol=="", NA, feature..data$Symbol)
feature..data<- feature..data[!(is.na(feature..data$Symbol) & is.na(feature..data$Entrez_Gene_ID)),]

###############################

normalized.expr_final <- exp_values %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature..data, by = 'ID')

#for 21 samples:::: 1,23:25,(2:22)
normalized.expr_final <- normalized.expr_final[,c(1,123:125,2:122)]

#no.of probes
no_of_probes<-unique(normalized.expr_final$ID)
probe<-length(no_of_probes)

#no.of genes
genes<- normalized.expr_final[,c(2,3)]
final_gene<- nrow(unique(genes))

# fetch pheno data
pheno.data<- gse$GSE30029_series_matrix.txt.gz@phenoData@data
pheno..data <- pheno.data[,c(2,6,10)]
pheno..data<- as.data.frame(t(pheno..data))

#save phenodata
write.table(pheno..data,"/home/keerthana2/AML/ILLUMINA/GSE30029/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(normalized.expr_final,"/home/keerthana2/AML/ILLUMINA/GSE30029/norm.txt",sep=" ",row.names=F)

pheno <- read.table("/home/keerthana2/AML/ILLUMINA/GSE30029/pheno.txt")
norm <- read.table("/home/keerthana2/AML/ILLUMINA/GSE30029/norm.txt", fill=T)

GSE30029<- bind_rows(pheno, norm)
write.table(GSE30029,"/home/keerthana2/AML/ILLUMINA/GSE30029/Final_Matrix_GSE30029.txt",sep=" ",row.names=F, col.names = F)


