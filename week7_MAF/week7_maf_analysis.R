library(TCGAbiolinks)
library(maftools)
library(ggplot2)
setwd("/Users/jonathanle/Documents/qbio/analysis_data/")

clinic = data.table::fread("/Users/jonathanle/Documents/qbio/analysis_data/coad_clinical_data.csv",
                           data.table = F)

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] = "Tumor_Sample_Barcode"

#colnames(clinic) has length 76
#colnames(clinic) == "bcr_patient_barcode" is also length 76. It contains bools on whether each column name is equal to "bcr_patient_barcode"
#there is one true

mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

maf_dataframe = data.table::fread("/Users/jonathanle/Documents/qbio/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F)

maf_object = read.maf(maf = maf_dataframe, clinicalData = clinic, isTCGA = T)


oncoplot(maf = maf_object,
         top = 10) 

ggsave("~/Documents/qbio/qbio_data_analysis_jonathanl/week7_MAF/oncoplot_10.png")

#The KRAS gene provides instructions for making a protein called K-Ras that is part of a signaling pathway known as the RAS/MAPK pathway. 
#The protein relays signals from outside the cell to the cell's nucleus. 
#These signals instruct the cell to grow and divide (proliferate) or to mature and take on specialized functions (differentiate).

clinic = maf_object@clinical.data
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "young", "old")

young_patient_ids = clinic[clinic$age_category == "young",]

young_maf = subsetMaf(maf_object, young_patient_ids$Tumor_Sample_Barcode)

old_maf = subsetMaf(maf_object, clinic[clinic$age_category == "old",]$Tumor_Sample_Barcode)

coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young Patients MAF", 
           m2Name = "Old Patients MAF")

ggsave("/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/week7_MAF/cooncoplot_young_old.png")


lollipopPlot(maf_object, gene = "KRAS")

ggsave("/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/week7_MAF/KRAS_lollipop_plot.png")

lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = "Young Patients MAF",
              m2_name = "Old Patients MAF",
              gene = "KRAS")

ggsave("/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/week7_MAF/KRAS_lollipop2_plot.png")

b = 7
c = 2
d = 35
e = 37
f = 42

geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

geneB_maf <- subsetMaf(maf = maf_object,
                       genes = "KRAS")

mut_bc_geneA = geneA_maf@clinical.data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@clinical.data$Tumor_Sample_Barcode

num_mut_geneA = length(mut_bc_geneA) #213
num_mut_geneB = length(mut_bc_geneB) #163

mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB) #78

num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB

num_neither_mutation = 397 - num_mut_geneA - num_mut_geneB + num_mut_geneAB

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneA_only,
                         num_mut_geneB_only,
                         num_neither_mutation), 
                         nrow=2)

# view the contingency table
contig_table

fe_results <- fisher.test(contig_table)
fe_results

#p-value = 0.06543. Since this is > 0.05, we cannot reject the Null hypothesis, so it is not likely that the 
#TP53 and KRAS have co-mutations (or that they are dependent on each other)


if (!require(foreach)) install.packages("foreach")
if (!require(doParallel)) install.packages("doParallel")
if (!require(reshape2)) install.packages("reshape2")

library(foreach)
library(doParallel)

# how many computations you can do at once
registerDoParallel(detectCores())

# turns the data frame into the following structure:
# rows: patient IDs
# columns: genes
# cells: number of mutations in that gene
maf_object@data$Index = 1
maf_cast = reshape2::dcast(maf_object@data, Tumor_Sample_Barcode ~ Hugo_Symbol)

# Compute the signifiances
ComputePairwiseMutation_Parallel = function(maf_cast, geneA){
  
  # Tumor_Sample_Barcode is the first column name
  all_genes <- unique(colnames(maf_cast)[-1])
  # TRUE if has mutation, FALSE otherwise
  # the important speedup: this is much faster than identifying names
  # and calling intersect
  has_mutA = (maf_cast[[geneA]] != 0)
  n = length(all_genes)
  
  # with parallelization
  foreach(i = 1:n, .combine = "c") %dopar% {
    gene = all_genes[i]
    has_mutB = (maf_cast[[gene]] != 0)
    
    # we can use the built-in table() function
    contig_table = table(has_mutA, has_mutB)
    fe_results <- fisher.test(contig_table)
    pval = fe_results$p.value
    
    names(pval) = gene
    return(pval)
  }
}

p_vals = ComputePairwiseMutation_Parallel(maf_cast, "KRAS")


results = data.frame(
  gene = names(p_vals),
  pvalues = p_vals
)

# 2. adjust p-values and create the padj column
# BH is the Benjamini-Hochberg method
results$padj = p.adjust(p = results$pvalues, method = "BH")

# 3. sort your results data frame
row_order = order(results$pvalues)
results = results[row_order,]

# 4. create your sig_results data frame
sig_results = results[results$padj < 0.05,]

data.table::fwrite(results,
                   file = "/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/week7_MAF/KRAS_pairwise_results.csv")
