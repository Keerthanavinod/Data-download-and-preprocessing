library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidyr)
setwd("/home/keerthana/AML/data/")


# Load the expression data
expr_data <-fread("/home/keerthana/AML/data/edited final matrices/AML_GSE976.txt")
expr_data<- row_to_names(expr_data,4, remove_row = T, remove_rows_above = T)
expr_data<-expr_data %>%
  mutate(`Gene.Symbol` = strsplit(as.character(`Gene.Symbol`), " /// ")) %>%
  unnest(`Gene.Symbol`)
expr_data<-expr_data[,c(3,5:22)]
expr_data<-expr_data%>% mutate_at(c(2:19), as.numeric)
expr_data<- aggregate(expr_data[, 2:19], by = list(expr_data$Gene.Symbol), FUN = mean)
colnames(expr_data)[1]<-'Gene.Symbol'


#taking the rowmeans
expr_data$Exp_mean<- rowMeans(expr_data[,c(2:19)])

#separate into quantiles adn pick the top quantile
expr_data$Quantiles<- cut(expr_data$Exp_mean , breaks=quantile(expr_data$Exp_mean),
                          labels=1:4, include.lowest=TRUE)

#filter topquantile
top_quantile<- filter(expr_data, expr_data$Quantiles=="4")

#remove HKGs 
read<-fread("~/Downloads/House-keeping genes - Sheet1.csv")
top_quantile<-filter(top_quantile,!(top_quantile$Gene.Symbol %in% read$GENES))

top_quantile<- filter(top_quantile, substr(`Gene.Symbol`,1,2)!= "RP")
top_quantile<- filter(top_quantile, substr(`Gene.Symbol`,1,3)!= "SNO")
top_quantile<- filter(top_quantile, substr(`Gene.Symbol`,1,3)!= "LOC")
top_quantile<- filter(top_quantile, substr(`Gene.Symbol`,1,5)!= "GAPDH")
top_quantile<- filter(top_quantile, substr(`Gene.Symbol`,1,3)!= "ACT")


#merge with crispr dep data
dep<- fread("/home/keerthana2/AML/DEPMAP/crispr_gene_dependency_final.csv")
dep<- dep[-c(1:6),]
#dependency<- row_to_names(dep,1, remove_row = T, remove_rows_above = T)
colnames(dep)[1]<-"Gene.Symbol"
Gene.Symbol<-dep$Gene.Symbol
dependency<-dep[,c(2:27)] %>% mutate_if(is.character, as.numeric)
dependency<- round(dependency[,c(1:26)], digits=10)
dependency<-cbind(Gene.Symbol,dependency)

merged<- merge(top_quantile, dependency, by="Gene.Symbol")


#CALCULATE DEP MEAN
merged$Dep_mean<- rowMeans(merged[,c(22:47)])

#select required columns
final_data<- merged[,c(1,20,48,22:47)]

#count the number of column values >0.7 rowise
final_data$counts <- rowSums(final_data[,c(4:29)] > 0.7)

#filter values of counts>6
final_data<-filter(final_data, final_data$counts>=13)
final_data<- arrange(final_data, desc(counts))

#save the file
write.csv(final_data, "/home/keerthana/AML/GENE_ANALYSIS/RMA/GSE976.csv", row.names = F)

