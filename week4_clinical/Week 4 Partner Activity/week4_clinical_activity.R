setwd("/Users/jonathanle/Documents/qbio/qbio_data_analysis_jonathanl/analysis_data/")

clinic = read.csv("coad_clinical_data.csv")

sum(is.na(clinic$stage_event_pathologic_stage)) #no NAs

clinic_copy = clinic #make copy of clinic

#replace NAs in days to death column with days to last follow up
clinic_copy$days_to_death = ifelse(is.na(clinic$days_to_death), clinic_copy$days_to_last_follow_up, clinic_copy$days_to_death)

#remove NAs
clinic_copy = clinic_copy[!(is.na(clinic_copy$days_to_death)),]
#remove blank stages
clinic_copy = clinic_copy[clinic_copy$stage_event_pathologic_stage != "",]
#remove patients with no weight
clinic_copy = clinic_copy[!is.na(clinic_copy$weight),]
#remove patients with no height
clinic_copy = clinic_copy[!is.na(clinic_copy$height),]

#calculate BMI
clinic_copy$BMI = clinic_copy$weight / (clinic_copy$height/100)^2
#create categorical variable for BMI (for survival plot)
clinic_copy$weight_status = ifelse(clinic_copy$BMI >= 25, "overweight", 
                                   ifelse(clinic_copy$BMI < 25 & clinic_copy$BMI >= 18.5, "normal", "underweight"))

#create boxplot comparing BMI with pathological staging
boxplot(clinic_copy$BMI ~ clinic_copy$stage_event_pathologic_stage)

#survival analysis for pathogical staging
clinic_copy$death_event = ifelse(clinic_copy$vital_status == "Dead", 1, 0)

surv_object <- Surv(time = clinic_copy$days_to_death, 
                    event = clinic_copy$death_event)

pathological_staging_fit <- surv_fit(surv_object ~ clinic_copy$stage_event_pathologic_stage, data = clinic_copy)

survplot = ggsurvplot(pathological_staging_fit, 
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

ggsave("../kmplot_by_pathological_stage.png", plot = p, width = 12, height = 9)

#survival analysis for BMI
clinic_copy$death_event = ifelse(clinic_copy$vital_status == "Dead", 1, 0)

surv_object <- Surv(time = clinic_copy$days_to_death, 
                    event = clinic_copy$death_event)

BMI_fit <- surv_fit(surv_object ~ clinic_copy$weight_status, data = clinic_copy)

survplot = ggsurvplot(BMI_fit, 
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

ggsave("../kmplot_by_BMI.png", plot = p, width = 12, height = 9)



