library(tidyverse)
library(readxl)
library(latex2exp)
library(patchwork)
library(pROC)
library(xtable)
library(predtools)
library(arsenal)
library(gridExtra)
library(leaps)
library(bestglm)

# function for assessing in-sample performance with logistic regression
# all available variables are used, with dengue status as the response variable

in_sample_performance <- function(cur_data, threshold){
  yhat <- rep(NA, nrow(cur_data))
  sites <- unique(cur_data$Site)
  if(length(sites) > 1){ # leave-one-site-out CV
    for(j in 1:length(sites)){
      train <- cur_data %>%
        filter(Site != sites[j]) %>%
        select(-c(Site, Dataset))
      
      test <- cur_data %>%
        filter(Site == sites[j]) %>%
        select(-c(Site, Dataset))
      
      logistic_mod <- glm(Dengue ~ ., 
                          family = binomial,
                          data = train)
      
      yhat[cur_data$Site == sites[j]] <- 
        predict(logistic_mod, newdata = test, type = "response")
      
    }
  } else { # 10-fold CV
    folds <- sample(1:10, nrow(cur_data), replace=T)
    
    for(j in 1:10){
      train <- cur_data[folds != j,] %>%
        select(-c(Site, Dataset))
      test <- cur_data[folds == j,] %>%
        select(-c(Site, Dataset))
      
      logistic_mod <- glm(Dengue ~ ., 
                          family = binomial,
                          data = train)
      
      yhat[folds == j] <- 
        predict(logistic_mod, newdata = test, type = "response")
    }
  }
  
  sensitivity_ci <- binom.test(sum(yhat > threshold & 
                                     cur_data$Dengue == 1), 
                               sum(cur_data$Dengue == 1))$conf.int
  sensitivity <- sum((yhat > threshold) & 
                       cur_data$Dengue == 1)/sum(cur_data$Dengue)
  
  specificity_ci <- binom.test(sum(yhat <= threshold & 
                                     cur_data$Dengue == 0), 
                               sum(cur_data$Dengue == 0))$conf.int
  specificity <- sum((yhat <= threshold) & 
                       cur_data$Dengue == 0)/sum(1 - cur_data$Dengue)
  
  ppv_ci <- binom.test(sum(yhat > threshold & 
                             cur_data$Dengue == 1), 
                       sum(yhat > threshold))$conf.int
  ppv <- sum((yhat > threshold) & 
               cur_data$Dengue == 1)/sum(yhat > threshold)
  
  
  npv_ci <- binom.test(sum(yhat <= threshold & 
                             cur_data$Dengue == 0), 
                       sum(yhat <= threshold))$conf.int
  npv <- sum((yhat <= threshold) & 
               cur_data$Dengue == 0)/sum(yhat <= threshold)
  
  auc_ci <- round(c(ci.auc(roc(cur_data$Dengue, yhat))), 2)
  auc <- paste(auc_ci[2], " (", auc_ci[1], ", ", auc_ci[3], ")", sep="")
  
  sensitivity <- paste(round(sensitivity, 3), " (", 
                       round(sensitivity_ci[1], 3), ", ",
                       round(sensitivity_ci[2], 3), ")",
                       sep = "")
  specificity <- paste(round(specificity, 3), " (", 
                       round(specificity_ci[1], 3), ", ",
                       round(specificity_ci[2], 3), ")",
                       sep = "")
  ppv <- paste(round(ppv, 3), " (", 
               round(ppv_ci[1], 3), ", ",
               round(ppv_ci[2], 3), ")", 
               sep = "")
  npv <- paste(round(npv, 3), " (", 
               round(npv_ci[1], 3), ", ",
               round(npv_ci[2], 3), ")", 
               sep = "")
  
  return(c(sensitivity, specificity, ppv, npv, auc))
}

generalizability <- function(train, test, thresh){
  # fit a model on the training data
  logistic_mod <- glm(Dengue ~ ., 
                      family = binomial,
                      data = train)
  
  # predict on the test data
  preds <- predict(logistic_mod, newdata = test, type = "response")
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  sensitivity <- sum(preds > thresh & test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  specificity <- sum(preds <= thresh & test$Dengue == 0)/sum(test$Dengue == 0)
  
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    ppv <- paste(round(ppv, 3), " (", 
                 round(ppv_ci[1], 3), ", ",
                 round(ppv_ci[2], 3), ")", 
                 sep = "")
  }
  
  
  npv_ci <- binom.test(sum(preds <= thresh & 
                             test$Dengue == 0), 
                       sum(preds <= thresh))$conf.int
  
  npv <- sum(preds <= thresh & test$Dengue == 0)/sum(preds <= thresh)
  
  sensitivity <- paste(round(sensitivity, 3), " (", 
                       round(sensitivity_ci[1], 3), ", ",
                       round(sensitivity_ci[2], 3), ")",
                       sep = "")
  specificity <- paste(round(specificity, 3), " (", 
                       round(specificity_ci[1], 3), ", ",
                       round(specificity_ci[2], 3), ")",
                       sep = "")
  
  npv <- paste(round(npv, 3), " (", 
               round(npv_ci[1], 3), ", ",
               round(npv_ci[2], 3), ")", 
               sep = "")
  auc <- paste(auc_ci[2], " (", auc_ci[1], ", ", auc_ci[3], ")", sep="")
  
  return(c(sensitivity, specificity, ppv, npv, auc))
}



# Trying a different subset

data02 <- read_excel("../data/set-02/dengue-data-02.xlsx") %>%
  mutate(BP = (as.numeric(VSBPDIA) + as.numeric(VSBPSYS))/2) %>%
  dplyr::select(Site, Sex, `Age (Converted from Days to Years)`,
                LBLEUKRS, LBTHRRS, LBLIMFRS, VSTEMP, BP, LBHEMORS, `FINAL Category`) %>%
  rename(Age = `Age (Converted from Days to Years)`,
         WBC = LBLEUKRS, 
         PLT = LBTHRRS,
         LYMPH = LBLIMFRS,
         Temp = VSTEMP,
         Hb = LBHEMORS,
         Dengue = `FINAL Category`) %>%
  mutate(Dengue = ifelse(Dengue == "Dengue", 1, 0),
         Sex = 2 - Sex) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 2")

data03 <- read_excel("../data/set-03/dengue-data-03.xlsx") %>%
  mutate(Site = 1) %>%
  dplyr::select(Site, sex, age, wbc, plt, lymph, bt, mbp, hb, dengue) %>%
  rename(Age = age, WBC = wbc, PLT = plt, Dengue = dengue,
         LYMPH = lymph, Temp = bt, BP = mbp,
         Hb = hb, Sex = sex) %>%
  mutate(Hb = Hb/10) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 3")

data05 <- read_excel("../data/set-05/dengue-data-05.xls") %>%
  mutate(Site = 1,
         BP = (VitalSign_DBP + VitalSign_SBP)/2) %>%
  dplyr::select(Site, Gender, AgeEnrol, Result_WBC, Result_platelet, 
                Lymphocyte, VitalSign_Temp, BP,
                Result_Hb,
                ConfirmDengue_YesNo_FourCriteria) %>%
  rename(Age = AgeEnrol, WBC = Result_WBC, PLT = Result_platelet, 
         LYMPH = Lymphocyte, Temp = VitalSign_Temp,
         Hb = Result_Hb,
         Dengue = ConfirmDengue_YesNo_FourCriteria,
         Sex = Gender) %>%
  mutate(Dengue = ifelse(Dengue == "No dengue", 0, 1),
         Sex = ifelse(Sex == "Male", 1, 0),
         Dataset = "Dataset 5") %>%
  drop_na()



#### Data 2

data_2_results <- data.frame(train = c(), 
                             test = c(),
                             model = c(),
                             sense = c(),
                             spec = c(),
                             ppv = c(),
                             npv = c(),
                             auc = c())

# best subset selection

X <- data02 %>%
  select(-c(Site, Dengue, Dataset))
y <- data02$Dengue
bss_res <- bestglm(cbind(X, y), IC = "AIC")

## data2 and data3: Age, WBC, PLT, LYMPH, Temp, Hb; data 5: WBC, PLT, LYMPH, Temp

# Age, WBC, PLT, LYMPH, Temp, Hb
names(bss_res$BestModel$coefficients)

# in-sample performance
data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 2", "Original", 
                          data02 %>%
                            select(Site, Age, WBC, PLT, Dengue, Dataset) %>%
                            in_sample_performance(0.45)))

colnames(data_2_results) <- c("train", "test", "model", "sens", "spec", "ppv", "npv", "auc")


data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 2", "Full", 
                          data02 %>%
                            in_sample_performance(0.45)))

data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 2", "BSS", 
                          data02 %>%
                            select(Site, Dataset, Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            in_sample_performance(0.45)))

# a small model is also beneficial in terms of interpretation and usability for practitioners


# generalizability

data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 3", "Original", 
                          data02 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03, 0.45)))

data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 3", "Full", 
                          data02 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03, 0.45)))

data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 3", "BSS", 
                          data02 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03, 0.45)))


data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 5", "Original", 
                          data02 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data05, 0.45)))

data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 5", "Full", 
                          data02 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data05, 0.45)))

data_2_results <- rbind(data_2_results,
                        c("Dataset 2", "Dataset 5", "BSS", 
                          data02 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data05, 0.45)))





#### Data 3

data_3_results <- data.frame(train = c(), 
                             test = c(),
                             model = c(),
                             sense = c(),
                             spec = c(),
                             ppv = c(),
                             npv = c(),
                             auc = c())

# best subset selection

X <- data03 %>%
  select(-c(Site, Dengue, Dataset))
y <- data03$Dengue
bss_res <- bestglm(cbind(X, y), IC = "AIC")

# Age, WBC, PLT, LYMPH, Temp, Hb
names(bss_res$BestModel$coefficients)

# in-sample performance
data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 3", "Original", 
                          data03 %>%
                            select(Site, Age, WBC, PLT, Dengue, Dataset) %>%
                            in_sample_performance(0.21)))

colnames(data_3_results) <- c("train", "test", "model", "sens", "spec", "ppv", "npv", "auc")


data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 3", "Full", 
                          data03 %>%
                            in_sample_performance(0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 3", "BSS", 
                          data03 %>%
                            select(Site, Dataset, Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            in_sample_performance(0.21)))

# a small model is also beneficial in terms of interpretation and usability for practitioners


# generalizability

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 2", "Original", 
                          data03 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data02, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 2", "Full", 
                          data03 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data02, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 2", "BSS", 
                          data03 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data02, 0.21)))


data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 5", "Original", 
                          data03 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data05, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 5", "Full", 
                          data03 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data05, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 5", "BSS", 
                          data03 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data05, 0.21)))








### Data 5


data_5_results <- data.frame(train = c(), 
                             test = c(),
                             model = c(),
                             sense = c(),
                             spec = c(),
                             ppv = c(),
                             npv = c(),
                             auc = c())

# best subset selection

X <- data05 %>%
  select(-c(Site, Dengue, Dataset))
y <- data05$Dengue
bss_res <- bestglm(cbind(X, y), IC = "AIC")

# WBC, PLT,  LYMPH, Temp
names(bss_res$BestModel$coefficients)

# in-sample performance
data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 5", "Original", 
                          data05 %>%
                            select(Site, Age, WBC, PLT, Dengue, Dataset) %>%
                            in_sample_performance(0.43)))

colnames(data_5_results) <- c("train", "test", "model", "sens", "spec", "ppv", "npv", "auc")


data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 5", "Full", 
                          data05 %>%
                            in_sample_performance(0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 5", "BSS", 
                          data05 %>%
                            select(Site, Dataset, Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            in_sample_performance(0.43)))

# a small model is also beneficial in terms of interpretation and usability for practitioners


# generalizability

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 2", "Original", 
                          data05 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data02, 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 2", "Full", 
                          data05 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data02, 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 2", "BSS", 
                          data05 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data02, 0.43)))



data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 2 (Age > 16)", "Original", 
                          data05 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data02 %>% filter(Age > 16), 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 2 (Age > 16)", "Full", 
                          data05 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data02 %>% filter(Age > 16), 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 2 (Age > 16)", "BSS", 
                          data05 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data02 %>% filter(Age > 16), 0.43)))




data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 3", "Original", 
                          data05 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03, 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 3", "Full", 
                          data05 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03, 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 3", "BSS", 
                          data05 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03, 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 3 (Age > 16)", "Original", 
                          data05 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03 %>% filter(Age > 16), 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 3 (Age > 16)", "Full", 
                          data05 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03 %>% filter(Age > 16), 0.43)))

data_5_results <- rbind(data_5_results,
                        c("Dataset 5", "Dataset 3 (Age > 16)", "BSS", 
                          data05 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03 %>% filter(Age > 16), 0.43)))



### Supplementary Table 5

# convert to LaTeX for paper
rbind(data_2_results, data_3_results, data_5_results) %>%
  xtable() %>%
  print(include.rownames = F)
