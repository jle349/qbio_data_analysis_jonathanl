library(TCGAbiolinks)
library(SummarizedExperiment)

#Exercise 2.1
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
#GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)

counts = assays(sum_exp)$"HTSeq - Counts"

dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")

head(colData(sum_exp)) #rows in colData are patients, and columns are cancer data
head(assays(sum_exp)$"HTSeq - Counts") #rows in sum_exp are patients and columns are gene expression counts

colnames(colData(sum_exp)) #age at diagnosis

colData(sum_exp)$age_at_diagnosis[1:10]

colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis / 365

colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis >= 50, "Old", "Young")

#2.2
head(rowData(sum_exp))
dim(rowData(sum_exp))

"MSH2" %in% rowData(sum_exp)$external_gene_name #true
"MSH6" %in% rowData(sum_exp)$external_gene_name #true

#2.3
assays(sum_exp)$"HTSeq - Counts"[20:25,30:35]

#gene A = MSH2
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "MSH2")
sum(geneA_id_mask)
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask]

#gene B = MSH6
geneB_id_mask = (rowData(sum_exp)$external_gene_name == "MSH6")
sum(geneB_id_mask)
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask]

#ensembl gene ID is a row in assays(sum_exp)$"HTSeq - Counts"

min(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,])
max(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,])

summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB,])

plot(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,],
     assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB,],
     xlab = "MSH2", # remember to rename axes!
     ylab = "MSH6"
)

bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)

age_cat_no_NAs = colData(sum_exp)$age_category[!is.na(colData(sum_exp)$age_category)]
length(age_cat_no_NAs)

num_rows_colData = dim(colData(sum_exp))[1]
num_NAs = num_na

length_age_cat_no_NAs = length(age_cat_no_NAs)

num_rows_colData - num_NAs == length_age_cat_no_NAs #true

#copy of colData with NAs removed from age category
colData_no_age_NA = colData(sum_exp)[!is.na(colData(sum_exp)$age_category),]

length(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA,]) #521 patients
#copy of assays data with same NAs removed
counts_data_no_age_na = assays(sum_exp)$"HTSeq - Counts"[,!bool_age_na]

identical(rownames(colData_no_age_NA), colnames(counts_data_no_age_na)) #true

geneA_counts = counts_data_no_age_na[ensembl_geneA,]

length(age_cat_no_NAs) == length(geneA_counts) #true

boxplot(geneA_counts ~ colData_no_age_NA$age_category, xlab = "Age Category", ylab = "Counts for MSH2")

#3.1
#1
#The HTSeqs - Counts dataframe is accessed using assays(sum_exp)$"HTSeq - Counts"
#The rows represent each individual gene and the columns represent the counts for each patient sample.

#2
#Use rowData(sum_exp). This function gives you more information about each row, specifically the external gene name associated with each ensembl gene id

#3
#Use colData(sum_exp). This function gives you more information about each column and describes additional data about each patient beyond the gene counts