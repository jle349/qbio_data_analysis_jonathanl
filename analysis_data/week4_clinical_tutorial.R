clin_query <- GDCquery(project = "TCGA-COAD", data.category = "Clinical", file.type = "xml")
# Only use this line ONCE! Comment out after you have downloaded the data. 
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

clinic_copy = clinic

str(clinic)
head(clinic)
colnames(clinic)
clinic$lymph_node_examined_count


plot(clinic$age_at_initial_pathologic_diagnosis, clinic$weight, xlab = "Age at Diagnosis", ylab = "weight")

unique(clinic$race_list)
par(mar=c(10,5,1,1))
boxplot(clinic$age_at_initial_pathologic_diagnosis ~ clinic$race_list,
        las = 2, 
        cex.axis = 0.5)

clinic_copy$race_list = ifelse(clinic_copy$race_list == "", "No data", as.character(clinic_copy$race_list))
boxplot(clinic_copy$age_at_initial_pathologic_diagnosis ~ clinic_copy$race_list,
        las = 2, 
        cex.axis = 0.5)

min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)

summary(clinic$age_at_initial_pathologic_diagnosis)

clinic_copy$age_category = ifelse(clinic_copy$age_at_initial_pathologic_diagnosis >= 50, "old", "young")

old_clinic = clinic_copy[clinic_copy$age_category == "old",]
young_clinic = clinic_copy[clinic_copy$age_category == "young",]

num_old = nrow(old_clinic)
num_young = nrow(young_clinic)
nrow(clinic)
#they are equal

young_patient_ids = clinic_copy[clinic_copy$age_category == "young",]$patient_id

# This is a good sanity check (making sure you have the same number of patients)
length(young_patient_ids) == num_young #Sanity checks are always important!

old_patient_ids = clinic_copy[clinic_copy$age_category == "old",]$patient_id

length(old_patient_ids) == num_old


###Exercise 2.8
clinic[1,1] #first row, first column
clinic[1,]  #entire first row
clinic[2:5,] #rows 2 - 5
clinic[,3] #just 3rd column


young_clinic_one_line = clinic_copy[clinic_copy$age_at_initial_pathologic_diagnosis < 50,]

identical(dim(young_clinic), dim(young_clinic_one_line))

clinic_copy$days_to_death = ifelse(is.na(clinic$days_to_death), clinic_copy$days_to_last_follow_up, clinic_copy$days_to_death)

clinic_copy$death_event = ifelse(clinic_copy$vital_status == "Dead", 1, 0)

surv_object <- Surv(time = clinic_copy$days_to_death, 
                    event = clinic_copy$death_event)

race_fit <- surv_fit( surv_object ~ clinic_copy$race_list, data = clinic_copy)

survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")

p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p

ggsave("../kmplot_by_race.png", plot = p, width = 12, height = 9)


write.csv(clinic, "/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/analysis_data/coad_clinical_data.csv", row.names = F)

