
library(tidyverse)
library(ggplot2)
library(pROC) 
library(ggplot2)
library(dvmisc)

#1. extract final genes and CIs and plot 
final_gene<- read.csv("./output/final_gene.csv", row.names=1) %>% 
  arrange(desc(abs(Coefficient))) %>% 
  mutate(across(c(Coefficient, Lower_CI, Upper_CI), ~ round(.,2))) %>% 
  filter(!Gene %in% c("sex_emr", "age_diag_calculated", "(Intercept)"))

#plot
coef_plot<- ggplot(final_gene, aes(x=reorder(Gene, abs(Coefficient)), y=Coefficient))+ 
  geom_point(color="red", size=3)+ 
  geom_errorbar(aes(ymin= Lower_CI, ymax= Upper_CI), width=0.4)+ 
  ylim(-6,6.5)+
  theme_bw()+
  coord_flip()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkblue", alpha= 0.7, size=1) + 
  labs(title= paste("25-gene panel"), y= "Coefficient", x= "Gene")+ 
  theme(
    plot.title = element_text(size = 20, face = "bold.italic", color = "red"),  # Increase title size
    axis.title.x = element_text(size = 20, face = "bold.italic", color = "darkblue"),  
    axis.title.y = element_text(size = 20, face = "bold.italic", color = "darkblue"),  
    axis.text.x = element_text(size = 18),  
    axis.text.y = element_text(size = 18)   
  )

par(oma = c(1,1,1,1))
jpeg("./plots/coef_plot.jpg", width = 1000, height = 600)
coef_plot
dev.off()

#2. generate roc curves with different models in testing set

#read wp score in the testing set 
score_test<- read.csv("./output/score_test.csv", row.names=1) 
score_train<- read.csv("./output/score_train.csv", row.names=1) 
score<- rbind(score_test, score_train)

#function to generate logistic regression models 
generate_models_probs <- function(df) {
  # logistic regression models
  logit.wp_score <- glm(response_abstracted_cr ~ wp_score, data = df, family = "binomial")
  logit.age.sex.iss <- glm(response_abstracted_cr ~ age_diag_calculated + iss_derived_m + sex_emr, data = df, family = "binomial")
  logit.age.sex.iss.wp_score <- glm(response_abstracted_cr ~ age_diag_calculated + iss_derived_m + sex_emr + wp_score, data = df, family = "binomial")
  
  # predicted probabilities (for ROC curve)
  prob.wp_score <- predict(logit.wp_score, type = "response")
  prob.age.sex.iss <- predict(logit.age.sex.iss, type = "response")
  prob.age.sex.iss.wp_score <- predict(logit.age.sex.iss.wp_score, type = "response")
  
  # return the list of predicted probabilities
  list(
    prob.wp_score = prob.wp_score,
    prob.age.sex.iss = prob.age.sex.iss,
    prob.age.sex.iss.wp_score = prob.age.sex.iss.wp_score
  )
}

#function to create roc curves 
generate_roc_plot_logit <- function(df, models_probs) {
  colors <- c("red", "blue", "orange")  # Colors for ROC curves
  
  # Store model names
  model_names <- c("Model 1: WP Score", "Model 2: Age + Sex + ISS", "Model 1 + Model 2")
  
  # Initialize legend text
  legend_text <- c()
  
  # Loop through models and plot ROC curves
  for (i in seq_along(models_probs)) {
    model_probs <- models_probs[[i]]
    
    # Generate ROC curve with confidence intervals
    roc_curve <- roc(df$response_abstracted_cr, model_probs, ci = TRUE)
    
    if (i == 1) {
      # For the first model, initialize the plot with plot.roc()
      plot.roc(
        roc_curve, 
        col = colors[i], 
        lwd = 3, 
        xlab = "1 - Specificity", 
        ylab = "Sensitivity",
        legacy.axes = TRUE, 
        asp = NA
      )
    } else {
      # For subsequent models, add to the existing plot
      plot.roc(
        roc_curve, 
        col = colors[i], 
        lwd = 3, 
        add = TRUE, 
        legacy.axes = TRUE
      )
    }
    
    # Construct legend text with AUC and CI
    auc <- round(roc_curve$auc, 1)
    ci_lower <- round(roc_curve$ci[1], 1)
    ci_upper <- round(roc_curve$ci[2], 1)
    legend_text <- c(legend_text, paste0(model_names[i], "; AUC: ", auc, " (", ci_lower, "-", ci_upper, ")"))
  }
  
  # Add legend
  legend("bottomright", legend = legend_text, col = colors, lwd = 3, cex = 1.3, bty = "n")
}


models_probs<- generate_models_probs(score_test)
par(oma = c(1,1,1,1))
jpeg("./plots/roc.jpg", width = 800, height = 700)
generate_roc_plot_logit(score_test, models_probs)
dev.off()

#3. generate roc curves in different drug subgroups 
#rbind score_train and score_test
score_train<- read.csv("./output/score_train.csv", row.names=1)
score_all<- rbind(score_train, score_test)

#create subgroups by drug class
df_response<- readRDS("Y:/bw_folder/dissertation/5hmc_tx/analysis_drug/treatment_labeling/df_response.RData")
#only examine in testing set for now
score_test<- score_test %>% 
  mutate(mrn= as.character(mrn)) %>% 
  left_join(df_response %>% 
              select(mrn, dtq), by="mrn")

score_doub<- score_test %>% 
  filter(dtq== "doub")
score_trip<- score_test %>% 
  filter(dtq== "trip")

#apply the functions to different subgroups 
#for doublet 
model_probs_doub<- generate_models_probs(score_doub)
par(oma = c(1,1,1,1))
jpeg("./plots/roc_doub.jpg", width = 800, height = 700)
generate_roc_plot_logit(score_doub, model_probs_doub)
dev.off()

#for triplet
model_probs_trip<- generate_models_probs(score_trip)
par(oma = c(1,1,1,1))
jpeg("./plots/roc_trip.jpg", width = 800, height = 700)
generate_roc_plot_logit(score_trip, model_probs_trip)
dev.off()

#4. generate roc curves in different racial subgroups 
score_ea<- score_test %>% 
  filter(race_composite== "White")
score_aa<- score_test %>% 
  filter(race_composite== "Black/African-American")

#apply the functions to different subgroups 
#for ea 
model_probs_ea<- generate_models_probs(score_ea)
par(oma = c(1,1,1,1))
jpeg("./plots/roc_ea.jpg", width = 800, height = 700)
generate_roc_plot_logit(score_ea, model_probs_ea)
dev.off()

#for aa

model_probs_aa<- generate_models_probs(score_aa)
par(oma = c(1,1,1,1))
jpeg("./plots/roc_aa.jpg", width = 800, height = 700)
generate_roc_plot_logit(score_aa, model_probs_aa)
dev.off()

#5. comparing AUCs

mod1 <- roc(score_test$response_abstracted_cr, models_probs$prob.wp_score, ci = TRUE)
mod2<-  roc(score_test$response_abstracted_cr, models_probs$prob.age.sex.iss, ci = TRUE)
mod3<-  roc(score_test$response_abstracted_cr, models_probs$prob.age.sex.iss.wp_score, ci = TRUE)

roc.test(mod1, mod2, method="d")
roc.test(mod2, mod3, method="d")
roc.test(mod1, mod3, method="d")

#comparing in triplet subgroup 
mod1 <- roc(score_trip$response_abstracted_cr, model_probs_trip$prob.wp_score, ci = TRUE)
mod2<-  roc(score_trip$response_abstracted_cr, model_probs_trip$prob.age.sex.iss, ci = TRUE)
mod3<-  roc(score_trip$response_abstracted_cr, model_probs_trip$prob.age.sex.iss.wp_score, ci = TRUE)

roc.test(mod1, mod2, method="d")
roc.test(mod2, mod3, method="d")
roc.test(mod1, mod3, method="d")
