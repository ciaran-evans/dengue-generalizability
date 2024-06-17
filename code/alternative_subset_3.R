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


data01 <- read_excel("../data/set-01/dengue-data-01.xls") %>%
  dplyr::select(SiteNo, Sex, Age, WBC, PLT, LYMCount, HCT, 
                AST, ALT, Lab_Confirmed_Dengue) %>%
  rename(Site = SiteNo,
         Dengue = Lab_Confirmed_Dengue,
         LYMPH = LYMCount) %>%
  mutate(Dengue = 2 - Dengue,
         Sex = 2 - Sex) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 1")

data03 <- read_excel("../data/set-03/dengue-data-03.xlsx") %>%
  mutate(Site = 1) %>%
  dplyr::select(Site, sex, age, wbc, plt, lymph, hct, ast, alt, dengue) %>%
  rename(Age = age, WBC = wbc, PLT = plt, Dengue = dengue,
         LYMPH = lymph, HCT = hct, AST = ast, ALT = alt,
         Sex = sex) %>%
  mutate(HCT = 100 * HCT) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 3")


data04 <- read_excel("../data/set-04/dengue-data-04.xlsx") %>%
  dplyr::select(sex, age2, wbc_m3, wbc_m1, platelets_m3, platelets_m1, 
                lymph_m3, lymph_m1, max_hct_m3, max_hct_m1,
                ast_m3, ast_m1, alt_m3, alt_m1, dengue) %>%
  pivot_longer(c(wbc_m3, wbc_m1, platelets_m3, platelets_m1,
                 lymph_m3, lymph_m1, max_hct_m3, max_hct_m1,
                 ast_m3, ast_m1, alt_m3, alt_m1),
               names_to = c(".value", "Day"),
               names_sep = "_m") %>%
  rename(Dengue = dengue, WBC = wbc, PLT = platelets,
         LYMPH = lymph, HCT = max_hct, AST = ast, ALT = alt,
         Sex = sex) %>%
  separate(age2, c("age1", "age2")) %>%
  mutate(Age = (as.numeric(age1) + as.numeric(age2))/2,
         WBC = WBC/1000, PLT = PLT/1000,
         Sex = 2 - Sex,
         Dataset = paste("Dataset 4, Day -", Day, sep=""),
         Site = 1) %>%
  dplyr::select(Site, Sex, Age, WBC, PLT, 
                LYMPH, HCT, AST, ALT, Dengue, Dataset)

data04_1 <- data04 %>%
  filter(Dataset == "Dataset 4, Day -1")

data04_3 <- data04 %>%
  filter(Dataset == "Dataset 4, Day -3")



#### Data 1

data_1_results <- data.frame(train = c(), 
                             test = c(),
                             model = c(),
                             sense = c(),
                             spec = c(),
                             ppv = c(),
                             npv = c(),
                             auc = c())

# best subset selection

X <- data01 %>%
  select(-c(Site, Dengue, Dataset))
y <- data01$Dengue
bss_res <- bestglm(cbind(X, y), IC = "AIC")

## data 1: Age, WBC, PLT, LYMPH, HCT, AST
## data 3: Age, PLT, LYMPH, HCT
## data 4 day 1: Age, WBC, PLT, LYMPH, ALT
## data 4 day 3: Age, WBC, LYMPH, AST

names(bss_res$BestModel$coefficients)

# in-sample performance
data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 1", "Original", 
                          data01 %>%
                            select(Site, Age, WBC, PLT, Dengue, Dataset) %>%
                            in_sample_performance(0.33)))

colnames(data_1_results) <- c("train", "test", "model", "sens", "spec", "ppv", "npv", "auc")


data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 1", "Full", 
                          data01 %>%
                            in_sample_performance(0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 1", "BSS", 
                          data01 %>%
                            select(Site, Dataset, Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            in_sample_performance(0.33)))

# a small model is also beneficial in terms of interpretation and usability for practitioners


# generalizability

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 3", "Original", 
                          data01 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 3", "Full", 
                          data01 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 3", "BSS", 
                          data01 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 3 (Age < 16)", "Original", 
                          data01 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03 %>% filter(Age < 16), 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 3 (Age < 16)", "Full", 
                          data01 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03 %>% filter(Age < 16), 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 3 (Age < 16)", "BSS", 
                          data01 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03 %>% filter(Age < 16), 0.33)))



data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 4, Day -3", "Original", 
                          data01 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data04_3, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 4, Day -3", "Full", 
                          data01 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data04_3, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 4, Day -3", "BSS", 
                          data01 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data04_3, 0.33)))



data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 4, Day -1", "Original", 
                          data01 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data04_1, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 4, Day -1", "Full", 
                          data01 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data04_1, 0.33)))

data_1_results <- rbind(data_1_results,
                        c("Dataset 1", "Dataset 4, Day -1", "BSS", 
                          data01 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data04_1, 0.33)))






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

# Age, PLT, LYMPH, HCT
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
                        c("Dataset 3", "Dataset 1", "Original", 
                          data03 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data01, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 1", "Full", 
                          data03 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data01, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 1", "BSS", 
                          data03 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data01, 0.21)))


data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 4, Day -3", "Original", 
                          data03 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data04_3, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 4, Day -3", "Full", 
                          data03 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data04_3, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 4, Day -3", "BSS", 
                          data03 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data04_3, 0.21)))


data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 4, Day -1", "Original", 
                          data03 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data04_1, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 4, Day -1", "Full", 
                          data03 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data04_1, 0.21)))

data_3_results <- rbind(data_3_results,
                        c("Dataset 3", "Dataset 4, Day -1", "BSS", 
                          data03 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data04_1, 0.21)))








### Data 4, Day -3


data_4_3_results <- data.frame(train = c(), 
                             test = c(),
                             model = c(),
                             sense = c(),
                             spec = c(),
                             ppv = c(),
                             npv = c(),
                             auc = c())

# best subset selection

X <- data04_3 %>%
  select(-c(Site, Dengue, Dataset))
y <- data04_3$Dengue
bss_res <- bestglm(cbind(X, y), IC = "AIC")

#  Age, WBC, LYMPH, AST
names(bss_res$BestModel$coefficients)

# in-sample performance
data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 4, Day -3", "Original", 
                          data04_3 %>%
                            select(Site, Age, WBC, PLT, Dengue, Dataset) %>%
                            in_sample_performance(0.59)))

colnames(data_4_3_results) <- c("train", "test", "model", "sens", "spec", "ppv", "npv", "auc")


data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 4, Day -3", "Full", 
                          data04_3 %>%
                            in_sample_performance(0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 4, Day -3", "BSS", 
                          data04_3 %>%
                            select(Site, Dataset, Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            in_sample_performance(0.59)))

# a small model is also beneficial in terms of interpretation and usability for practitioners


# generalizability

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 1", "Original", 
                          data04_3 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data01, 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 1", "Full", 
                          data04_3 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data01, 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 1", "BSS", 
                          data04_3 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data01, 0.59)))


data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 3", "Original", 
                          data04_3 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03, 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 3", "Full", 
                          data04_3 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03, 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 3", "BSS", 
                          data04_3 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03, 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 3 (Age < 16)", "Original", 
                          data04_3 %>%
                            select(Age, WBC, PLT, Dengue) %>%
                            generalizability(data03 %>% filter(Age < 16), 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 3 (Age < 16)", "Full", 
                          data04_3 %>%
                            select(-c(Site, Dataset)) %>%
                            generalizability(data03 %>% filter(Age < 16), 0.59)))

data_4_3_results <- rbind(data_4_3_results,
                        c("Dataset 4, Day -3", "Dataset 3 (Age < 16)", "BSS", 
                          data04_3 %>%
                            select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                            generalizability(data03 %>% filter(Age < 16), 0.59)))





### Data 4, Day -1


data_4_1_results <- data.frame(train = c(), 
                               test = c(),
                               model = c(),
                               sense = c(),
                               spec = c(),
                               ppv = c(),
                               npv = c(),
                               auc = c())

# best subset selection

X <- data04_1 %>%
  select(-c(Site, Dengue, Dataset))
y <- data04_1$Dengue
bss_res <- bestglm(cbind(X, y), IC = "AIC")

#  Age, WBC, LYMPH, PLT, LYMPH, ALT
names(bss_res$BestModel$coefficients)

# in-sample performance
data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 4, Day -1", "Original", 
                            data04_1 %>%
                              select(Site, Age, WBC, PLT, Dengue, Dataset) %>%
                              in_sample_performance(0.64)))

colnames(data_4_1_results) <- c("train", "test", "model", "sens", "spec", "ppv", "npv", "auc")


data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 4, Day -1", "Full", 
                            data04_1 %>%
                              in_sample_performance(0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 4, Day -1", "BSS", 
                            data04_1 %>%
                              select(Site, Dataset, Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                              in_sample_performance(0.64)))

# a small model is also beneficial in terms of interpretation and usability for practitioners


# generalizability

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 1", "Original", 
                            data04_1 %>%
                              select(Age, WBC, PLT, Dengue) %>%
                              generalizability(data01, 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 1", "Full", 
                            data04_1 %>%
                              select(-c(Site, Dataset)) %>%
                              generalizability(data01, 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 1", "BSS", 
                            data04_1 %>%
                              select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                              generalizability(data01, 0.64)))


data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 3", "Original", 
                            data04_1 %>%
                              select(Age, WBC, PLT, Dengue) %>%
                              generalizability(data03, 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 3", "Full", 
                            data04_1 %>%
                              select(-c(Site, Dataset)) %>%
                              generalizability(data03, 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 3", "BSS", 
                            data04_1 %>%
                              select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                              generalizability(data03, 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 3 (Age < 16)", "Original", 
                            data04_1 %>%
                              select(Age, WBC, PLT, Dengue) %>%
                              generalizability(data03 %>% filter(Age < 16), 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 3 (Age < 16)", "Full", 
                            data04_1 %>%
                              select(-c(Site, Dataset)) %>%
                              generalizability(data03 %>% filter(Age < 16), 0.64)))

data_4_1_results <- rbind(data_4_1_results,
                          c("Dataset 4, Day -1", "Dataset 3 (Age < 16)", "BSS", 
                            data04_1 %>%
                              select(Dengue, names(bss_res$BestModel$coefficients)[-1]) %>%
                              generalizability(data03 %>% filter(Age < 16), 0.64)))



### Supplementary Table 6


# convert to LaTeX for paper
rbind(data_1_results, data_3_results, data_4_3_results, data_4_1_results) %>%
  xtable() %>%
  print(include.rownames = F)
