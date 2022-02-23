setwd("/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/analysis_data/")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(DESeq2)

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
sum_exp <- GDCprepare(query)

patients_data = colData(sum_exp)
NA_patients = is.na(patients_data$age_at_index)

counts = assays(sum_exp)$"HTSeq - Counts"


patients_data = patients_data[!NA_patients,] #removed NAs

counts = counts[,!NA_patients]

patients_data$age_category = ifelse(patients_data$age_at_index < 50, "young", "old")
patients_data$age_category = factor(patients_data$age_category, levels = c("young", "old"))

counts_row_sums = rowSums(counts)

low_counts_mask = ifelse(counts_row_sums >= 10, TRUE, FALSE)
sum(low_counts_mask)

counts_copy = counts
counts = counts[low_counts_mask,]

dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patients_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

row_order = order(results$padj)
results = results[row_order,]

#ENSG00000200087 is a gene that is more highly expressed in old patients, has it has a log2FoldChange of -3.11991.
#This gene's full name is SNORA73B and functions as a regulator of chromatin function.

log2FoldChange_threshold = 1
padj_threshold = 0.05

genes_pass_log_threshold = ((results$log2FoldChange > log2FoldChange_threshold) || (results$log2FoldChange < -1 * log2FoldChange_threshold))

genes_pass_padj_threshold = results$log2FoldChange < padj_threshold

results_subset = results[genes_pass_log_threshold & genes_pass_padj_threshold,]

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

plot(x = results_subset$log2FoldChange,
     y = -log10(results_subset$padj),
     xlab = "Log 2 Fold Change (Young over Old)", # be sure the specify that it's young over old!
     ylab = "Negative log 10 Adjusted p-value",
     pch = 20)

abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot

write.csv(x = results,
          file = "C:/Users/jonathanle/Documents/qbio/week6_DESeq2/results.csv",
          row.names = FALSE)
