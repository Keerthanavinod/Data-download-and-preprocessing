setwd("/home/keerthana2/AML/Methylation profiling data/data/")
library(GEOquery)
library(data.table)
library(dplyr)

#untar
#untar("/home/keerthana2/AML/Methylation profiling data/data/GSE19515/GSE19515_RAW.tar", exdir = "GSE19515/unzip/")


#gse download
gse <- getGEO("GSE19515", GSEMatrix = TRUE)

#feature data extraction
featuredata<-gse$GSE19515_series_matrix.txt.gz@featureData@data
feature.data<- featuredata[,c(16,5,6,11,7,9)]
feature.data<- tibble::rownames_to_column(feature.data, "ID")

#beta value extraction
beta_values<- as.data.frame(gse$GSE19515_series_matrix.txt.gz@assayData$exprs)
beta_values<- tibble::rownames_to_column(beta_values, "ID")

#merge beta values and feature data
merged<- merge(beta_values, feature.data, by ="ID")
final_data<- merged[,c(29:34,2:28)]
colnames(final_data)[1]<- "ID"
#no.of probes
no_of_probes<-unique(final_data$ID)
probe<-length(no_of_probes)

#no.of genes
genes<- final_data[,c(2,3)]
final_gene<- nrow(unique(genes))

# fetch pheno data
phenodata<- gse$GSE19515_series_matrix.txt.gz@phenoData@data
pheno.data<- phenodata[,c(2,8,11)]
pheno..data<- as.data.frame(t(pheno.data))

#save phenodata
write.table(pheno..data,"/home/keerthana2/AML/Methylation profiling data/data/GSE19515/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(final_data,"/home/keerthana2/AML/Methylation profiling data/data/GSE19515/beta_values.txt",sep=" ",row.names=F)

##read pheno and beta data
pheno <- read.table("/home/keerthana2/AML/Methylation profiling data/data/GSE19515/pheno.txt")
beta<- read.table("/home/keerthana2/AML/Methylation profiling data/data/GSE19515/beta_values.txt", fill=T)

GSE19515<- bind_rows(pheno, beta)
GSE19515<-GSE19515%>% mutate_all(na_if, "")

##save final matrix
write.table(GSE19515,"/home/keerthana2/AML/Methylation profiling data/data/GSE19515/Final_Matrix_GSE19515.txt",sep=" ",row.names=F, col.names = F)
