setwd("/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/analysis_data/")
library(TCGAbiolinks)
library(maftools)
library(DESeq2)
library(ggplot2)
library(survival)
library(survminer)

### Mutation Analysis
clinic = data.table::fread("/Users/jonathanle/Documents/qbio/analysis_data/coad_clinical_data.csv",
                           data.table = F)

colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] = "Tumor_Sample_Barcode"
clinic$days_to_death = ifelse(is.na(clinic$days_to_death), clinic$days_to_last_follow_up, clinic$days_to_death)
clinic$death_event = ifelse(clinic$vital_status == "Dead", 1, 0)
clinic$age_category = ifelse(clinic$age_at_initial_pathologic_diagnosis < 50, "young", "old")


maf_dataframe = data.table::fread("/Users/jonathanle/Documents/qbio/analysis_data/GDCdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F)

maf_object = read.maf(maf = maf_dataframe, clinicalData = clinic, isTCGA = T)

oncoplot(maf = maf_object, top = 5) #get oncoplot of top 5 most mutated genes

young_patient_ids = clinic[clinic$age_category == "young",]

young_maf = subsetMaf(maf_object, young_patient_ids$Tumor_Sample_Barcode)
old_maf = subsetMaf(maf_object, clinic[clinic$age_category == "old",]$Tumor_Sample_Barcode)


coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young Patients MAF", 
           m2Name = "Old Patients MAF",
           outer_mar = 3,
           barcode_mar = 0.5,
           gene_mar = 0.3,
           titleFontSize = 0.6)


### DESEq Analysis

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

APC_mask = rowData(sum_exp)$external_gene_name == "APC"
TTN_mask = rowData(sum_exp)$external_gene_name == "TTN"
TP53_mask = rowData(sum_exp)$external_gene_name == "TP53"
KRAS_mask = rowData(sum_exp)$external_gene_name == "KRAS"

APC_ensembl_id = rowData(sum_exp)$ensembl_gene_id[APC_mask]
TTN_ensembl_id = rowData(sum_exp)$ensembl_gene_id[TTN_mask]
TP53_ensembl_id = rowData(sum_exp)$ensembl_gene_id[TP53_mask]
KRAS_ensembl_id = rowData(sum_exp)$ensembl_gene_id[KRAS_mask]


#get counts for top 3 most mutated genes
top_4_mutated_genes_counts = counts[c(APC_ensembl_id,TTN_ensembl_id, TP53_ensembl_id, KRAS_ensembl_id),]

dds = DESeqDataSetFromMatrix(countData = top_4_mutated_genes_counts,
                             colData = patients_data,
                             design = ~age_category)

dds_obj = DESeq(dds)
resultsNames(dds_obj)
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "young", "old"))

row_order = order(results$padj)
results = results[row_order,]

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

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

###Survival Plots
#filter each clinic dataframe to get only patients that have that mutation
KRAS_maf = subsetMaf(maf_object, genes = "KRAS")
KRAS_maf_patient_IDs = KRAS_maf@clinical.data$Tumor_Sample_Barcode
KRAS_clinic = clinic[clinic$Tumor_Sample_Barcode %in% c(KRAS_maf_patient_IDs),]

APC_maf = subsetMaf(maf_object, genes = "APC")
APC_maf_patient_IDs = APC_maf@clinical.data$Tumor_Sample_Barcode
APC_clinic = clinic[clinic$Tumor_Sample_Barcode %in% c(APC_maf_patient_IDs),]

TTN_maf = subsetMaf(maf_object, genes = "TTN")
TTN_maf_patient_IDs = TTN_maf@clinical.data$Tumor_Sample_Barcode
TTN_clinic = clinic[clinic$Tumor_Sample_Barcode %in% c(TTN_maf_patient_IDs),]

TP53_maf = subsetMaf(maf_object, genes = "TP53")
TP53_maf_patient_IDs = TP53_maf@clinical.data$Tumor_Sample_Barcode
TP53_clinic = clinic[clinic$Tumor_Sample_Barcode %in% c(TP53_maf_patient_IDs),]


#survival plots

surv_object <- Surv(time = KRAS_clinic$days_to_death, 
                    event = KRAS_clinic$death_event)

KRAS_fit <- surv_fit(surv_object ~ KRAS_clinic$age_category, data = KRAS_clinic)

survplot = ggsurvplot(KRAS_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p_kras = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_kras

surv_object <- Surv(time = APC_clinic$days_to_death, 
                    event = APC_clinic$death_event)

APC_fit <- surv_fit(surv_object ~ APC_clinic$age_category, data = APC_clinic)

survplot = ggsurvplot(APC_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p_APC = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_APC

surv_object <- Surv(time = TTN_clinic$days_to_death, 
                    event = TTN_clinic$death_event)

TTN_fit <- surv_fit(surv_object ~ TTN_clinic$age_category, data = TTN_clinic)

survplot = ggsurvplot(TTN_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p_TTN = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_TTN

surv_object <- Surv(time = TP53_clinic$days_to_death, 
                    event = TP53_clinic$death_event)

TP53_fit <- surv_fit(surv_object ~ TP53_clinic$age_category, data = TP53_clinic)

survplot = ggsurvplot(TP53_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p_TP53 = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p_TP53



