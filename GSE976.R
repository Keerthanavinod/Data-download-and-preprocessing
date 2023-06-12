if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pd.clariom.s.human")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pd.hugene.1.0.st.v1")

detach("package:oligo", unload = TRUE)


setwd("~/AML/data")
#library(oligo)
library(affy)
library(GEOquery)
library(dplyr)
library(tidyverse)

setwd("~/AML/data")
#untar the compressed files
untar("~/AML/data/GSE976/GSE976_RAW.tar", exdir="GSE976/CEL")


#listing all .cel files in the folder
#list.files("~/AML/")
#list.files("~/AML/data/CEL/GSE196107")

#read .cel files by oligo
#celfiles <- list.files("~/AML/data/GSE976/CEL", full = TRUE)
#rawData <- read.celfiles(celfiles)

#read .cel files by affy
rawData<-ReadAffy(celfile.path ="GSE976/CEL") 

## ----rma normalize-----------------------------------------------------------------
normalized.data <- rma(rawData)
normalized.expr <- as.data.frame(exprs(normalized.data))
names(normalized.expr) <- substring(names(normalized.expr),1,8)


##download GSE matrix file
gse <- getGEO("GSE976", GSEMatrix = TRUE)

# fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE976_series_matrix.txt.gz@featureData@data

# subset
###############################
feature..data <- feature.data[,c(1,7,11,12)]
###############################

normalized.expr_final <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature..data, by = 'ID')

#for 21 samples:::: 1,23:25,(2:22)
normalized.expr_final <- normalized.expr_final[,c(1,26:28,2:25)]

normalized.expr_final<- normalized.expr_final[!(normalized.expr_final$`Gene Symbol`=="" & normalized.expr_final$ENTREZ_GENE_ID==""),]

normalized.expr_final$ENTREZ_GENE_ID<- ifelse(normalized.expr_final$ENTREZ_GENE_ID=="", NA, normalized.expr_final$ENTREZ_GENE_ID)

#no.of probes
no_of_probes<-unique(normalized.expr_final$ID)
probe<-length(no_of_probes)

#no.of genes
genes<- normalized.expr_final[,c(3,4)]
genes<- strsplit(genes$`Gene Symbol`,split= '///' )
unlisting_gene<- data.frame(matrix(unlist(genes), byrow=TRUE),stringsAsFactors=FALSE)
final_gene<- nrow(unique(unlisting_gene))

# fetch pheno data
pheno.data<- gse$GSE976_series_matrix.txt.gz@phenoData@data
pheno..data <- pheno.data[,c(2,6,8)]
pheno..data<- as.data.frame(t(pheno..data))

#save phenodata
write.table(pheno..data,"~/AML/data/GSE976/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(normalized.expr_final,"~/AML/data/GSE976/norm.txt",sep=" ",row.names=F)

pheno <- read.table("~/AML/data/GSE976/pheno.txt")
norm <- read.table("~/AML/data/GSE976/norm.txt", fill=T)

GSE976<- bind_rows(pheno, norm)
write.table(GSE976,"~/AML/data/GSE976/Final_Matrix_GSE976.txt",sep=" ",row.names=F, col.names = F)


