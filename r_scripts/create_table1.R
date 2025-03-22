
library(tidyverse)
library(table1)
library(sjPlot)
library(survival)
library(survminer)

#data prep
df<-read.csv("./data/df_allrace_new.csv", row.names=1, check.names=FALSE) %>% 
  mutate(kappa_lambda_clean_cat= relevel(factor(kappa_lambda_clean_cat), ref = "normal"), 
         age_diag_cat= relevel(factor(age_diag_cat), ref= "32-59"), 
         response_abstracted_cr = case_when(response_abstracted_cr == "<CR" ~ 0,
                                            response_abstracted_cr == ">=CR" ~ 1,
                                            TRUE ~ NA_real_),
         response_abstracted_cr= factor(response_abstracted_cr, levels= c(0,1)), 
         iss_derived= as.character(iss_derived), 
         iss_bi= ifelse(iss_derived_m %in% c("2", "3"), "2/3", iss_derived_m))

#merge treatment information 
df$mrn<- as.character(df$mrn)
regimen<- readRDS("Y:/bw_folder/dissertation/5hmc_tx/analysis_drug/treatment_labeling/df_response.RData") %>% 
  select(mrn, dtq) %>% 
  mutate(dtq_3= ifelse((dtq== "quad"| is.na(dtq)), "quad/other", dtq))
df<- df %>% 
  mutate(mrn= as.character(mrn)) %>% 
  left_join(regimen, by="mrn")


#1.frequency of IMWG response ------------------------------------------------

table(df$response_abstracted, useNA="ifany")
round(prop.table(table(df$response_abstracted, useNA="ifany")), 3)

#by treatment received 
cross_tab<- table(df$response_abstracted, df$dtq, useNA = "ifany") 
percent<- prop.table(cross_tab, 2)*100

formatted_table <- matrix(paste(cross_tab, "(", round(percent, 1), "%)", sep=""), 
                          nrow = nrow(cross_tab), 
                          dimnames = dimnames(cross_tab))
formatted_table_df <- as.data.frame(formatted_table)

formatted_table_df


#2.population characteristics, table 1 ------------------------------------------------
#freq and proportions 
table1(~ age_diag_calculated+ age_diag_cat+ race_composite + sex_emr + educ_cat + iss_derived_m + as.character(ldh_bi) 
       + as.character(egfr_bi) + as.character(kappa_lambda_clean_cat)+ as.character(dtq_3)+ asct_abstracted
       | response_abstracted_cr, 
       data= df[!is.na(df$response_abstracted_cr),])

#get age-adjusted ORs 
df$iss_derived<- as.character(df$iss_derived)

#response 
table(df$response_abstracted_cr, useNA="ifany")
prop.table(table(df$response_abstracted_cr, useNA="ifany"))
#age
mod<- glm(response_abstracted_cr~ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

mod<- glm(response_abstracted_cr~ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#race
df$race_composite<- factor(df$race_composite)
df$race_composite <- relevel(df$race_composite, ref = "White")
mod<- glm(response_abstracted_cr~ race_composite+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#sex
mod<- glm(response_abstracted_cr~ sex_emr+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#education 
mod<- glm(response_abstracted_cr~ educ_cat+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#iss 
mod<- glm(response_abstracted_cr~ iss_derived+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#ldh 
mod<- glm(response_abstracted_cr~ ldh_bi+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#egfr
mod<- glm(response_abstracted_cr~ egfr_bi+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)
#sflc
mod<- glm(response_abstracted_cr~ kappa_lambda_clean_cat+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#dtq_3
mod<- glm(response_abstracted_cr~ dtq_3+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#asct
mod<- glm(response_abstracted_cr~ asct_abstracted+ age_diag_calculated, 
          data= df, family = binomial("logit"))
tab_model(mod, transform = "exp", digits = 1)

#3. km plots of response and os and pfs for validity ------------------------------------------------

#km plot by treatment response 
os_tx <- survfit(Surv(survival_months_composite, vital_status_composite) ~ response_abstracted_cr, data = df, conf.type = "log-log")
pfs_tx <- survfit(Surv(pfs_months, pfs_status) ~ response_abstracted_cr, data = df, conf.type = "log-log")

kms <- list()

# Generate the plot
kms[[1]] <- ggsurvplot(os_tx, data = df,
                       censor.shape = "|", censor.size = 10, pval = TRUE, pval.size = 7,
                       font.tickslab = c(25, "plain"),
                       subtitle = "OS", font.subtitle = c(25, "bold.italic", "darkblue"),
                       palette = c("red", "blue"),
                       xlim = c(0, 135), xlab = "Time since diagnosis (months)", font.x = c(25, "bold.italic", "darkblue"), break.x.by = 24,
                       font.y = c(25, "bold.italic", "darkblue"),
                       ggtheme = theme_light(base_size = 25), 
                       risk.table = TRUE,
                       risk.table.y.text = FALSE, 
                       risk.table.title.font = list(size = 10, face = "bold.italic", color = "red"),
                       legend="none")
kms[[2]] <- ggsurvplot(pfs_tx, data = df,
                       censor.shape = "|", censor.size = 10, pval = TRUE, pval.size = 7,
                       font.tickslab = c(25, "plain"),
                       subtitle = "PFS", font.subtitle = c(25, "bold.italic", "darkblue"),
                       palette = c("red", "blue"),
                       xlim = c(0, 135), xlab = "Time since diagnosis (months)", font.x = c(25, "bold.italic", "darkblue"), break.x.by = 24,
                       font.y = c(25, "bold.italic", "darkblue"),
                       ggtheme = theme_light(base_size = 25), 
                       risk.table = TRUE,
                       risk.table.y.text = FALSE, 
                       risk.table.title.font = list(size = 10, face = "bold.italic", color = "red"),
                       legend.title = "Response",
                       legend.labs = c("<CR", ">=CR"))

# Adjust the risk table theme
kms[[1]]$table <- kms[[1]]$table + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title= element_text(size = 20, face = "bold", color = "purple"))
kms[[2]]$table <- kms[[2]]$table + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank(), plot.title= element_text(size = 20, face = "bold", color = "purple"))

# Adjust the plot theme and add custom colors
kms[[2]]$plot <- kms[[2]]$plot + 
  theme(legend.position = c(0.95, 0.95),  # Right bottom position
        legend.justification = c("right", "top"),
        legend.key.size = unit(0.5, "cm"),  # Adjust legend key size
        legend.text = element_text(size = 10), # Adjust legend text size
        legend.title = element_text(size = 10), 
        legend.background = element_rect(fill = "beige"))

km_tx_os_pfs <- arrange_ggsurvplots(kms, print = FALSE, ncol = 2)

par(oma = c(1,1,1,1))
jpeg("./plots/km_tx_os_pfs.jpg", width = 1000, height = 600)
km_tx_os_pfs
dev.off()

#4. create a table 1 for population-specific factors in EAs and AAs 



df<- df %>% 
  filter(!is.na(race_composite))
table1(~ age_diag_calculated+ age_diag_cat + sex_emr + educ_cat + iss_derived_m + as.character(ldh_bi) 
       + as.character(egfr_bi) + as.character(kappa_lambda_clean_cat)+ as.character(dtq_3)+ asct_abstracted + response_abstracted_cr
       | race_composite, 
       data= df[!is.na(df$race_composite),])

chisq.test(table(df$response_abstracted_cr, as.character(df$ldh_bi)))$p.value
chisq.test(table(df$response_abstracted_cr, as.character(df$egfr_bi)))$p.value
chisq.test(table(df$response_abstracted_cr, as.character(df$kappa_lambda_clean_cat)))$p.value


