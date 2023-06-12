setwd("/home/keerthana2/AML/AGILENT/")
#loading libraries
library(limma)
library(GEOquery)
library(dplyr)
library(tidyverse)
library(stringr)


#untar the compressed files
untar("GSE19853/GSE19853_RAW.tar", exdir="GSE19853/TXT")

#read the target.txt file
Target<-readTargets(file = "GSE19853/Targets.txt", sep = "\t")

#reading the files
file<-read.maimages(Target, source = "agilent" , green.only = T, other.columns = "gIsWellAboveBG")

#background correction and normaliztion
normalized_exp<- backgroundCorrect(file, method = "normexp")
normalized_exp<-normalizeBetweenArrays(normalized_exp)

#extracting the expression values
normalized_exp_final<- as.data.frame(normalized_exp[[1]])

#setting column names
names(normalized_exp_final) <- substring(names(normalized_exp_final), 43,51)

#extracting probenames from genes
genes<- as.data.frame(normalized_exp[[3]])
ProbeName<-as.data.frame(genes$ProbeName)
colnames(ProbeName)<- "ID"

#merging ProbeName and norm_exp_final
normalized.exp.final<- cbind(ProbeName, normalized_exp_final)

#seriesmatrix download
gse <- getGEO("GSE19853", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE19853_series_matrix.txt.gz@featureData@data

feature..data <- feature.data[,c(4,10)]

#read file homo
file_homo<- read.csv("~/AML/Homo genes.csv")
feature..data$GENE_SYMBOL <- ifelse(feature..data$GENE_SYMBOL=="",NA, feature..data$GENE_SYMBOL)
feature..data<- feature..data[!(is.na(feature..data$GENE_SYMBOL)),]

colnames(file_homo)[2]<- "GENE_SYMBOL"
feature...data<- merge(feature..data,file_homo, by="GENE_SYMBOL", all.x= TRUE)
final_feature_data<- feature...data[,c(1,2,3)]
colnames(final_feature_data)[2]<-"ID"

normalized.expr_final <-inner_join(normalized.exp.final, final_feature_data, by = 'ID')

#rearranging
normalized.expr_final <- normalized.expr_final[,c(1,42,43,2:41)]

#no.of probes
no_of_probes<-unique(normalized.expr_final$ID)
probe<-length(no_of_probes)

#no.of genes
genes<- normalized.expr_final[,c(2,3)]
final_gene<- nrow(unique(genes))

# fetch pheno data
pheno.data<- gse$`GSE19853_series_matrix.txt.gz`@phenoData@data
pheno..data<- pheno.data[,c(2,6,34,33,35)]
pheno..data<- as.data.frame(t(pheno..data))

#save phenodata
write.table(pheno..data,"/home/keerthana2/AML/AGILENT/GSE19853/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(normalized.expr_final,"/home/keerthana2/AML/AGILENT/GSE19853/norm.txt",sep=" ",row.names=F)

##read pheno and norm data
pheno <- read.table("/home/keerthana2/AML/AGILENT/GSE19853/pheno.txt")
norm <- read.table("/home/keerthana2/AML/AGILENT/GSE19853/norm.txt", fill=T)

GSE19853<- bind_rows(pheno, norm)

#save final matrix
write.table(GSE19853,"/home/keerthana2/AML/AGILENT/GSE19853/Final_Matrix_GSE19853.txt",sep=" ",row.names=F, col.names = F)
