setwd("/home/keerthana2/AML/High_throughput_sequencing/data/")
library(dplyr)
library(GEOquery)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(data.table)
unzip<- untar("GSE48843/GSE48843_RAW.tar", exdir = "GSE48843/unzip")

#read .txt.gz files
data <- lapply(list.files('GSE48843/unzip/', pattern = '.txt',  full.names = TRUE),
               function(data) data.table::fread(file = data))


n<- length(data)
df<- data[[1]]
df<- df[,c(4,5)]
colnames(df)[c(1,2)]<-c('Gene','1')
df<- df %>% distinct(Gene,.keep_all = T)

##loop to get all sample information
for (i in 2:n){
  df2<-data[[i]]
  df2<-df2[,c(4,5)]
  colnames(df2)[c(1,2)]<-c("Gene",i)
  df2<- df2 %>% distinct(Gene,.keep_all = T)
  df<-full_join(df,df2,by='Gene')
}

#Convert from gene.symbol to ensembl.gene
genes.name<- df$Gene
x<-select(org.Hs.eg.db,
          keys = genes.name,
          columns=c("ENTREZID","SYMBOL"),
          keytype="SYMBOL")
x<-x %>% distinct(SYMBOL,.keep_all = T)
colnames(x)[1]<-"GENE SYMBOL"
colnames(df)[1]<-"GENE SYMBOL"
df_new<-merge(df,x,by="GENE SYMBOL",all.x=T)

##rearranging the columns
df_new<- df_new[,c(1,34,2:33)]

##converting to tpm counts
final_data<- df_new[,-c(1,2)]
df_tpm<-final_data%>% mutate(across(everything(), ~(./sum(.))*10**6))
df_tpm<-df_tpm %>% mutate(across(where(is.numeric), ~ round(., digits = 8)))
df_new<- cbind(df_new[,c(1,2)],df_tpm)

#gse download
gse <- getGEO("GSE48843", GSEMatrix = TRUE)

# fetch pheno data
pheno.data<- gse$GSE48843_series_matrix.txt.gz@phenoData@data
pheno..data <- pheno.data[,c(2,12)]
pheno..data<- as.data.frame(t(pheno..data))
sample<- colnames(pheno..data)
colnames(df_new)[3:34]<-sample
colnames(df_new)[c(1,2)]<-c("Gene.Symbol","ENTREZ_ID")


#save phenodata and exprs data
write.table(pheno..data,"GSE48843/pheno.txt",sep=" ",row.names=T, col.names =F)
write.table(df_new,"GSE48843/exprs.txt",sep=" ",row.names=F)

##read phenodata and exprs data
pheno <- read.table("GSE48843/pheno.txt")
exprs <- read.table("GSE48843/exprs.txt", fill=T, header = T)
gene<-unique(exprs$Gene.Symbol)
gene<-length(gene)
exprs1 <- read.table("GSE48843/exprs.txt", fill=T)
final_data<- bind_rows(pheno, exprs1)
final_data[final_data==0]<-"NA"

#save final matrix
write.table(final_data,"GSE48843/Final_Matrix_GSE48843.txt",row.names=F, col.names = F)
