setwd("/home/keerthana2/AML/Methylation profiling data/data/")
library(GEOquery)
library(data.table)
library(dplyr)

#untar
#untar("/home/keerthana2/AML/Methylation profiling data/data/GSE111915/GSE111915_RAW.tar", exdir = "GSE111915/unzip/")


#gse download
gse <- getGEO("GSE111915", GSEMatrix = TRUE)

#feature data extraction
featuredata<-gse$GSE111915_series_matrix.txt.gz@featureData@data
feature.data<- featuredata[,c(3,14,8,4)]

#beta value extraction
beta_values<- as.data.frame(gse$GSE111915_series_matrix.txt.gz@assayData$exprs)
#beta_values<- tibble::rownames_to_column(beta_values, "ID")

#merge beta values and feature data
final_data<- cbind(feature.data, beta_values)
#merged<- merge(beta_values, feature.data, by="ID")
#final_data<- merged[,c(1,83:86,2:82)]

#no.of probes
no_of_probes<-unique(final_data$`PROBE ID`)
probe<-length(no_of_probes)

#no.of genes
genes<- final_data[,c(2,3)]
final_gene<- nrow(unique(genes))

# fetch pheno data
phenodata<- gse$GSE111915_series_matrix.txt.gz@phenoData@data
pheno.data<- phenodata[,c(2,8,11,13)]
pheno..data<- as.data.frame(t(pheno.data))

#save phenodata
write.table(pheno..data,"/home/keerthana2/AML/Methylation profiling data/data/GSE111915/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(final_data,"/home/keerthana2/AML/Methylation profiling data/data/GSE111915/beta_values.txt",sep=" ",row.names=F)


#read pheno and beta data
pheno <- read.table("/home/keerthana2/AML/Methylation profiling data/data/GSE111915/pheno.txt")
beta<- read.table("/home/keerthana2/AML/Methylation profiling data/data/GSE111915/beta_values.txt", fill=T)

GSE111915<- bind_rows(pheno, beta)
GSE111915<-GSE111915%>% mutate_all(na_if, "")

##save final matrix
write.table(GSE111915,"/home/keerthana2/AML/Methylation profiling data/data/GSE111915/Final_Matrix_GSE111915.txt",sep=" ",row.names=F, col.names = F)
