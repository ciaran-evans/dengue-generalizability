library(tidyverse)
library(readxl)
library(latex2exp)
library(patchwork)
library(pROC)
library(xtable)
library(predtools)
library(arsenal)
library(gridExtra)
library(e1071)
library(rpart)

### Import data

data01 <- read_excel("../data/set-01/dengue-data-01.xls") %>%
  dplyr::select(SiteNo, Sex, Age, WBC, PLT, LYMCount, NEUCount, Temp, HCT, 
                Abdo, Weight, Vomiting, Lab_Confirmed_Dengue) %>%
  rename(Site = SiteNo,
         Dengue = Lab_Confirmed_Dengue,
         LYMPH = LYMCount,
         NEUT = NEUCount) %>%
  mutate(Dengue = 2 - Dengue,
         Vomiting = 2- Vomiting,
         Abdo = 2 - Abdo,
         Sex = 2 - Sex) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 1")


data03 <- read_excel("../data/set-03/dengue-data-03.xlsx") %>%
  mutate(Site = 1) %>%
  dplyr::select(Site, sex, age, wbc, plt, lymph, neut, bt, hct, "abdominal pain", 
                bw, "Vomiting/nausea", dengue) %>%
  rename(Age = age, WBC = wbc, PLT = plt, Dengue = dengue,
         LYMPH = lymph, NEUT = neut, Temp = bt, HCT = hct,
         Abdo = `abdominal pain`,
         Weight = bw, Vomiting = `Vomiting/nausea`,
         Sex = sex) %>%
  mutate(HCT = 100 * HCT) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 3")

data05 <- read_excel("../data/set-05/dengue-data-05.xls") %>%
  mutate(Site = 1) %>%
  dplyr::select(Site, Gender, AgeEnrol, Result_WBC, Result_platelet, 
                Lymphocyte, Neutrophils, VitalSign_Temp,
                Result_Hematocrit, GIT_abdominal_pain,
                Exam_weight, GIT_vomiting,
                ConfirmDengue_YesNo_FourCriteria) %>%
  rename(Age = AgeEnrol, WBC = Result_WBC, PLT = Result_platelet, 
         LYMPH = Lymphocyte, NEUT = Neutrophils, Temp = VitalSign_Temp,
         HCT = Result_Hematocrit, Abdo = GIT_abdominal_pain,
         Weight = Exam_weight, Vomiting = GIT_vomiting,
         Dengue = ConfirmDengue_YesNo_FourCriteria,
         Sex = Gender) %>%
  mutate(Dengue = ifelse(Dengue == "No dengue", 0, 1),
         Abdo = ifelse(Abdo == "No", 0, 1),
         Vomiting = ifelse(Vomiting == "No", 0, 1),
         Sex = ifelse(Sex == "Male", 1, 0),
         Dataset = "Dataset 5") %>%
  drop_na()

# combine the datasets into one
full_data <- data01 %>%
  rbind(data03) %>%
  rbind(data05) %>%
  mutate(Dengue = as.factor(Dengue))



##### Part 1: without adjusting age ranges


# Logistic regression

dataset_names <- unique(full_data$Dataset)
thresholds <- c(0.33, 0.21, 0.43)

# look at each training/test combination
# Keep the pairs where training and test are the same, as a baseline
# for comparison
generalizability_results <- expand.grid(dataset_names,
                                        dataset_names,
                                        c("Training threshold", 
                                          "Test threshold")) %>%
  rename(train_data = Var1, test_data = Var2, threshold = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA, AUC = NA,
         sensitivity_lower = NA, sensitivity_upper = NA,
         specificity_lower = NA, specificity_upper = NA,
         ppv_lower = NA, ppv_upper = NA,
         npv_lower = NA, npv_upper = NA,
         auc_lower = NA, auc_upper = NA)

for(r in 1:nrow(generalizability_results)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == generalizability_results$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == generalizability_results$test_data[r])
  
  # fit a model on the training data
  logistic_mod <- glm(Dengue ~ ., 
                      family = binomial,
                      data = train %>% select(-c(Site, Dataset)))
  
  # predict on the test data
  preds <- predict(logistic_mod, newdata = test, type = "response")
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  generalizability_results$AUC[r] <- auc_ci[2]
  generalizability_results$auc_lower[r] <- auc_ci[1]
  generalizability_results$auc_upper[r] <- auc_ci[3]
  
  # use either the threshold from the training set or from the test set
  if(generalizability_results$threshold[r] == "Training threshold"){
    thresh = thresholds[which(dataset_names == 
                                generalizability_results$train_data[r])]
  } else {
    thresh = thresholds[which(dataset_names == 
                                generalizability_results$test_data[r])]
  }
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  generalizability_results$sensitivity_lower[r] <- sensitivity_ci[1]
  generalizability_results$sensitivity_upper[r] <- sensitivity_ci[2]
  generalizability_results$Sensitivity[r] <- sum(preds > thresh & 
                                                   test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  generalizability_results$specificity_lower[r] <- specificity_ci[1]
  generalizability_results$specificity_upper[r] <- specificity_ci[2]
  generalizability_results$Specificity[r] <- sum(preds <= thresh & 
                                                   test$Dengue == 0)/sum(test$Dengue == 0)
  
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
    generalizability_results$ppv_lower[r] <- NA
    generalizability_results$ppv_upper[r] <- NA
    generalizability_results$PPV[r] <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    
    generalizability_results$ppv_lower[r] <- ppv_ci[1]
    generalizability_results$ppv_upper[r] <- ppv_ci[2]
    generalizability_results$PPV[r] <- ppv
  }                                    
  
  
  
  if(sum(preds <= thresh) == 0){
    generalizability_results$npv_lower[r] <- NA
    generalizability_results$npv_upper[r] <- NA
    generalizability_results$NPV[r] <- NA
  } else {
    npv_ci <- binom.test(sum(preds <= thresh & 
                               test$Dengue == 0), 
                         sum(preds <= thresh))$conf.int
    
    generalizability_results$npv_lower[r] <- npv_ci[1]
    generalizability_results$npv_upper[r] <- npv_ci[2]
    generalizability_results$NPV[r] <- sum(preds <= thresh & 
                                                 test$Dengue == 0)/sum(preds <= thresh)
  }
  
}

# don't want to train and test on different days for Dataset 4
generalizability_results <- generalizability_results %>% 
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))



# SVM

# look at each training/test combination
# Keep the pairs where training and test are the same, as a baseline
# for comparison
generalizability_results_svm <- expand.grid(dataset_names,
                                            dataset_names,
                                            c("Training threshold", 
                                              "Test threshold")) %>%
  rename(train_data = Var1, test_data = Var2, threshold = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA, AUC = NA,
         sensitivity_lower = NA, sensitivity_upper = NA,
         specificity_lower = NA, specificity_upper = NA,
         ppv_lower = NA, ppv_upper = NA,
         npv_lower = NA, npv_upper = NA,
         auc_lower = NA, auc_upper = NA)

for(r in 1:nrow(generalizability_results_svm)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == generalizability_results_svm$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == generalizability_results_svm$test_data[r])
  
  # fit a model on the training data
  
  svm_mod <- svm(Dengue ~ ., train %>% select(-c(Site, Dataset)), kernel = "radial",
                 probability = T)
  
  # predict on the test data
  preds <- attr(predict(svm_mod, newdata = test, probability = T), 
                "probabilities")[,1]
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  generalizability_results_svm$AUC[r] <- auc_ci[2]
  generalizability_results_svm$auc_lower[r] <- auc_ci[1]
  generalizability_results_svm$auc_upper[r] <- auc_ci[3]
  
  # use either the threshold from the training set or from the test set
  if(generalizability_results_svm$threshold[r] == "Training threshold"){
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_svm$train_data[r])]
  } else {
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_svm$test_data[r])]
  }
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  generalizability_results_svm$sensitivity_lower[r] <- sensitivity_ci[1]
  generalizability_results_svm$sensitivity_upper[r] <- sensitivity_ci[2]
  generalizability_results_svm$Sensitivity[r] <- sum(preds > thresh & 
                                                       test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  generalizability_results_svm$specificity_lower[r] <- specificity_ci[1]
  generalizability_results_svm$specificity_upper[r] <- specificity_ci[2]
  generalizability_results_svm$Specificity[r] <- sum(preds <= thresh & 
                                                       test$Dengue == 0)/sum(test$Dengue == 0)
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
    generalizability_results_svm$ppv_lower[r] <- NA
    generalizability_results_svm$ppv_upper[r] <- NA
    generalizability_results_svm$PPV[r] <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    
    generalizability_results_svm$ppv_lower[r] <- ppv_ci[1]
    generalizability_results_svm$ppv_upper[r] <- ppv_ci[2]
    generalizability_results_svm$PPV[r] <- ppv
  }
  # ppv_ci <- binom.test(sum(preds > thresh & 
  #                            test$Dengue == 1), 
  #                      sum(preds > thresh))$conf.int
  
  # generalizability_results_svm$ppv_lower[r] <- ppv_ci[1]
  # generalizability_results_svm$ppv_upper[r] <- ppv_ci[2]
  # generalizability_results_svm$PPV[r] <- sum(preds > thresh & 
  #                                          test$Dengue == 1)/sum(preds > thresh)
  
  
  
  if(sum(preds <= thresh) == 0){
    generalizability_results_svm$npv_lower[r] <- NA
    generalizability_results_svm$npv_upper[r] <- NA
    generalizability_results_svm$NPV[r] <- NA
  } else {
    npv_ci <- binom.test(sum(preds <= thresh & 
                               test$Dengue == 0), 
                         sum(preds <= thresh))$conf.int
    
    generalizability_results_svm$npv_lower[r] <- npv_ci[1]
    generalizability_results_svm$npv_upper[r] <- npv_ci[2]
    generalizability_results_svm$NPV[r] <- sum(preds <= thresh & 
                                                 test$Dengue == 0)/sum(preds <= thresh)
  }
  
  # npv_ci <- binom.test(sum(preds <= thresh & 
  #                            test$Dengue == 0), 
  #                      sum(preds <= thresh))$conf.int
  # 
  # generalizability_results_svm$npv_lower[r] <- npv_ci[1]
  # generalizability_results_svm$npv_upper[r] <- npv_ci[2]
  # generalizability_results_svm$NPV[r] <- sum(preds <= thresh & 
  #                                          test$Dengue == 0)/sum(preds <= thresh)
  
}

# don't want to train and test on different days for Dataset 4
generalizability_results_svm <- generalizability_results_svm %>% 
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))




# Decision tree

generalizability_results_dt <- expand.grid(dataset_names,
                                           dataset_names,
                                           c("Training threshold", 
                                             "Test threshold")) %>%
  rename(train_data = Var1, test_data = Var2, threshold = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA, AUC = NA,
         sensitivity_lower = NA, sensitivity_upper = NA,
         specificity_lower = NA, specificity_upper = NA,
         ppv_lower = NA, ppv_upper = NA,
         npv_lower = NA, npv_upper = NA,
         auc_lower = NA, auc_upper = NA)

for(r in 1:nrow(generalizability_results_dt)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == generalizability_results_dt$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == generalizability_results_dt$test_data[r])
  
  # fit a model on the training data
  
  dt_mod <- rpart(Dengue ~ ., data = train %>% select(-c(Site, Dataset)), method = "class")
  
  # predict on the test data
  preds <- attr(predict(svm_mod, newdata = test, probability = T), 
                "probabilities")[,1]
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  generalizability_results_dt$AUC[r] <- auc_ci[2]
  generalizability_results_dt$auc_lower[r] <- auc_ci[1]
  generalizability_results_dt$auc_upper[r] <- auc_ci[3]
  
  # use either the threshold from the training set or from the test set
  if(generalizability_results_dt$threshold[r] == "Training threshold"){
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_dt$train_data[r])]
  } else {
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_dt$test_data[r])]
  }
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  generalizability_results_dt$sensitivity_lower[r] <- sensitivity_ci[1]
  generalizability_results_dt$sensitivity_upper[r] <- sensitivity_ci[2]
  generalizability_results_dt$Sensitivity[r] <- sum(preds > thresh & 
                                                      test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  generalizability_results_dt$specificity_lower[r] <- specificity_ci[1]
  generalizability_results_dt$specificity_upper[r] <- specificity_ci[2]
  generalizability_results_dt$Specificity[r] <- sum(preds <= thresh & 
                                                      test$Dengue == 0)/sum(test$Dengue == 0)
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
    generalizability_results_dt$ppv_lower[r] <- NA
    generalizability_results_dt$ppv_upper[r] <- NA
    generalizability_results_dt$PPV[r] <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    
    generalizability_results_dt$ppv_lower[r] <- ppv_ci[1]
    generalizability_results_dt$ppv_upper[r] <- ppv_ci[2]
    generalizability_results_dt$PPV[r] <- ppv
  }
  # ppv_ci <- binom.test(sum(preds > thresh & 
  #                            test$Dengue == 1), 
  #                      sum(preds > thresh))$conf.int
  
  # generalizability_results_dt$ppv_lower[r] <- ppv_ci[1]
  # generalizability_results_dt$ppv_upper[r] <- ppv_ci[2]
  # generalizability_results_dt$PPV[r] <- sum(preds > thresh & 
  #                                          test$Dengue == 1)/sum(preds > thresh)
  
  
  
  if(sum(preds <= thresh) == 0){
    generalizability_results_dt$npv_lower[r] <- NA
    generalizability_results_dt$npv_upper[r] <- NA
    generalizability_results_dt$NPV[r] <- NA
  } else {
    npv_ci <- binom.test(sum(preds <= thresh & 
                               test$Dengue == 0), 
                         sum(preds <= thresh))$conf.int
    
    generalizability_results_dt$npv_lower[r] <- npv_ci[1]
    generalizability_results_dt$npv_upper[r] <- npv_ci[2]
    generalizability_results_dt$NPV[r] <- sum(preds <= thresh & 
                                                test$Dengue == 0)/sum(preds <= thresh)
  }
  
  # npv_ci <- binom.test(sum(preds <= thresh & 
  #                            test$Dengue == 0), 
  #                      sum(preds <= thresh))$conf.int
  # 
  # generalizability_results_dt$npv_lower[r] <- npv_ci[1]
  # generalizability_results_dt$npv_upper[r] <- npv_ci[2]
  # generalizability_results_dt$NPV[r] <- sum(preds <= thresh & 
  #                                          test$Dengue == 0)/sum(preds <= thresh)
  
}

# don't want to train and test on different days for Dataset 4
generalizability_results_dt <- generalizability_results_dt %>% 
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))



p1 <- generalizability_results %>%
  mutate(method = "Logistic regression") %>%
  rbind(generalizability_results_svm %>%
          mutate(method = "SVM")) %>%
  rbind(generalizability_results_dt %>%
          mutate(method = "Decision tree")) %>%
  filter(threshold == "Training threshold") %>%
  ggplot(aes(x = train_data, y = AUC, color = method, shape = method)) + 
  geom_errorbar(aes(ymin=auc_lower, ymax=auc_upper), width=0.2) +
  geom_point() +
  facet_wrap(~test_data) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7)) +
  labs(x = "Training dataset", y = "AUC for each test dataset", color = "Method",
       shape = "Method",
       title = "No age range adjustment")





### Part II: With age range adjustment

# Logistic regression

generalizability_results <- expand.grid(dataset_names,
                                        dataset_names,
                                        c("Training threshold", 
                                          "Test threshold")) %>%
  rename(train_data = Var1, test_data = Var2, threshold = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA, AUC = NA,
         sensitivity_lower = NA, sensitivity_upper = NA,
         specificity_lower = NA, specificity_upper = NA,
         ppv_lower = NA, ppv_upper = NA,
         npv_lower = NA, npv_upper = NA,
         auc_lower = NA, auc_upper = NA)

for(r in 1:nrow(generalizability_results)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == generalizability_results$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == generalizability_results$test_data[r])
  
  if(generalizability_results$train_data[r] %in% c("Dataset 1",
                                                   "Dataset 4, Day -1",
                                                   "Dataset 4, Day -3") & 
     generalizability_results$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age < 16)
  } else if (generalizability_results$train_data[r] == "Dataset 5" & 
             generalizability_results$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age > 16)
  }
  
  # fit a model on the training data
  logistic_mod <- glm(Dengue ~ ., 
                      family = binomial,
                      data = train %>% select(-c(Site, Dataset)))
  
  # predict on the test data
  preds <- predict(logistic_mod, newdata = test, type = "response")
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  generalizability_results$AUC[r] <- auc_ci[2]
  generalizability_results$auc_lower[r] <- auc_ci[1]
  generalizability_results$auc_upper[r] <- auc_ci[3]
  
  # use either the threshold from the training set or from the test set
  if(generalizability_results$threshold[r] == "Training threshold"){
    thresh = thresholds[which(dataset_names == 
                                generalizability_results$train_data[r])]
  } else {
    thresh = thresholds[which(dataset_names == 
                                generalizability_results$test_data[r])]
  }
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  generalizability_results$sensitivity_lower[r] <- sensitivity_ci[1]
  generalizability_results$sensitivity_upper[r] <- sensitivity_ci[2]
  generalizability_results$Sensitivity[r] <- sum(preds > thresh & 
                                                   test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  generalizability_results$specificity_lower[r] <- specificity_ci[1]
  generalizability_results$specificity_upper[r] <- specificity_ci[2]
  generalizability_results$Specificity[r] <- sum(preds <= thresh & 
                                                   test$Dengue == 0)/sum(test$Dengue == 0)
  
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
    generalizability_results$ppv_lower[r] <- NA
    generalizability_results$ppv_upper[r] <- NA
    generalizability_results$PPV[r] <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    
    generalizability_results$ppv_lower[r] <- ppv_ci[1]
    generalizability_results$ppv_upper[r] <- ppv_ci[2]
    generalizability_results$PPV[r] <- ppv
  }                               
  
  
  
  if(sum(preds <= thresh) == 0){
    generalizability_results$npv_lower[r] <- NA
    generalizability_results$npv_upper[r] <- NA
    generalizability_results$NPV[r] <- NA
  } else {
    npv_ci <- binom.test(sum(preds <= thresh & 
                               test$Dengue == 0), 
                         sum(preds <= thresh))$conf.int
    
    generalizability_results$npv_lower[r] <- npv_ci[1]
    generalizability_results$npv_upper[r] <- npv_ci[2]
    generalizability_results$NPV[r] <- sum(preds <= thresh & 
                                                 test$Dengue == 0)/sum(preds <= thresh)
  }
  
}

# don't want to train and test on different days for Dataset 4
generalizability_results <- generalizability_results %>% 
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))






# SVM

# look at each training/test combination
# Keep the pairs where training and test are the same, as a baseline
# for comparison
generalizability_results_svm <- expand.grid(dataset_names,
                                            dataset_names,
                                            c("Training threshold", 
                                              "Test threshold")) %>%
  rename(train_data = Var1, test_data = Var2, threshold = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA, AUC = NA,
         sensitivity_lower = NA, sensitivity_upper = NA,
         specificity_lower = NA, specificity_upper = NA,
         ppv_lower = NA, ppv_upper = NA,
         npv_lower = NA, npv_upper = NA,
         auc_lower = NA, auc_upper = NA)

for(r in 1:nrow(generalizability_results_svm)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == generalizability_results_svm$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == generalizability_results_svm$test_data[r])
  
  if(generalizability_results$train_data[r] %in% c("Dataset 1",
                                                   "Dataset 4, Day -1",
                                                   "Dataset 4, Day -3") & 
     generalizability_results$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age < 16)
  } else if (generalizability_results$train_data[r] == "Dataset 5" & 
             generalizability_results$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age > 16)
  }
  
  # fit a model on the training data
  
  svm_mod <- svm(Dengue ~ ., data = train %>% select(-c(Site, Dataset)), kernel = "radial",
                 probability = T)
  
  # predict on the test data
  preds <- attr(predict(svm_mod, newdata = test, probability = T), 
                "probabilities")[,1]
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  generalizability_results_svm$AUC[r] <- auc_ci[2]
  generalizability_results_svm$auc_lower[r] <- auc_ci[1]
  generalizability_results_svm$auc_upper[r] <- auc_ci[3]
  
  # use either the threshold from the training set or from the test set
  if(generalizability_results_svm$threshold[r] == "Training threshold"){
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_svm$train_data[r])]
  } else {
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_svm$test_data[r])]
  }
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  generalizability_results_svm$sensitivity_lower[r] <- sensitivity_ci[1]
  generalizability_results_svm$sensitivity_upper[r] <- sensitivity_ci[2]
  generalizability_results_svm$Sensitivity[r] <- sum(preds > thresh & 
                                                       test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  generalizability_results_svm$specificity_lower[r] <- specificity_ci[1]
  generalizability_results_svm$specificity_upper[r] <- specificity_ci[2]
  generalizability_results_svm$Specificity[r] <- sum(preds <= thresh & 
                                                       test$Dengue == 0)/sum(test$Dengue == 0)
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
    generalizability_results_svm$ppv_lower[r] <- NA
    generalizability_results_svm$ppv_upper[r] <- NA
    generalizability_results_svm$PPV[r] <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    
    generalizability_results_svm$ppv_lower[r] <- ppv_ci[1]
    generalizability_results_svm$ppv_upper[r] <- ppv_ci[2]
    generalizability_results_svm$PPV[r] <- ppv
  }
  # ppv_ci <- binom.test(sum(preds > thresh & 
  #                            test$Dengue == 1), 
  #                      sum(preds > thresh))$conf.int
  
  # generalizability_results_svm$ppv_lower[r] <- ppv_ci[1]
  # generalizability_results_svm$ppv_upper[r] <- ppv_ci[2]
  # generalizability_results_svm$PPV[r] <- sum(preds > thresh & 
  #                                          test$Dengue == 1)/sum(preds > thresh)
  
  
  
  if(sum(preds <= thresh) == 0){
    generalizability_results_svm$npv_lower[r] <- NA
    generalizability_results_svm$npv_upper[r] <- NA
    generalizability_results_svm$NPV[r] <- NA
  } else {
    npv_ci <- binom.test(sum(preds <= thresh & 
                               test$Dengue == 0), 
                         sum(preds <= thresh))$conf.int
    
    generalizability_results_svm$npv_lower[r] <- npv_ci[1]
    generalizability_results_svm$npv_upper[r] <- npv_ci[2]
    generalizability_results_svm$NPV[r] <- sum(preds <= thresh & 
                                                 test$Dengue == 0)/sum(preds <= thresh)
  }
  
  # npv_ci <- binom.test(sum(preds <= thresh & 
  #                            test$Dengue == 0), 
  #                      sum(preds <= thresh))$conf.int
  # 
  # generalizability_results_svm$npv_lower[r] <- npv_ci[1]
  # generalizability_results_svm$npv_upper[r] <- npv_ci[2]
  # generalizability_results_svm$NPV[r] <- sum(preds <= thresh & 
  #                                          test$Dengue == 0)/sum(preds <= thresh)
  
}

# don't want to train and test on different days for Dataset 4
generalizability_results_svm <- generalizability_results_svm %>% 
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))




# Decision tree

generalizability_results_dt <- expand.grid(dataset_names,
                                           dataset_names,
                                           c("Training threshold", 
                                             "Test threshold")) %>%
  rename(train_data = Var1, test_data = Var2, threshold = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, PPV = NA, NPV = NA, AUC = NA,
         sensitivity_lower = NA, sensitivity_upper = NA,
         specificity_lower = NA, specificity_upper = NA,
         ppv_lower = NA, ppv_upper = NA,
         npv_lower = NA, npv_upper = NA,
         auc_lower = NA, auc_upper = NA)

for(r in 1:nrow(generalizability_results_dt)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == generalizability_results_dt$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == generalizability_results_dt$test_data[r])
  
  if(generalizability_results$train_data[r] %in% c("Dataset 1",
                                                   "Dataset 4, Day -1",
                                                   "Dataset 4, Day -3") & 
     generalizability_results$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age < 16)
  } else if (generalizability_results$train_data[r] == "Dataset 5" & 
             generalizability_results$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age > 16)
  }
  
  # fit a model on the training data
  
  dt_mod <- rpart(Dengue ~ ., data = train %>% select(-c(Site, Dataset)), method = "class")
  
  # predict on the test data
  preds <- attr(predict(svm_mod, newdata = test, probability = T), 
                "probabilities")[,1]
  
  ## calculate performance metrics for predictions on the test data
  
  auc_ci <- round(c(ci.auc(roc(test$Dengue, preds))), 3)
  
  generalizability_results_dt$AUC[r] <- auc_ci[2]
  generalizability_results_dt$auc_lower[r] <- auc_ci[1]
  generalizability_results_dt$auc_upper[r] <- auc_ci[3]
  
  # use either the threshold from the training set or from the test set
  if(generalizability_results_dt$threshold[r] == "Training threshold"){
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_dt$train_data[r])]
  } else {
    thresh = thresholds[which(dataset_names == 
                                generalizability_results_dt$test_data[r])]
  }
  
  sensitivity_ci <- binom.test(sum(preds > thresh & 
                                     test$Dengue == 1), 
                               sum(test$Dengue == 1))$conf.int
  
  generalizability_results_dt$sensitivity_lower[r] <- sensitivity_ci[1]
  generalizability_results_dt$sensitivity_upper[r] <- sensitivity_ci[2]
  generalizability_results_dt$Sensitivity[r] <- sum(preds > thresh & 
                                                      test$Dengue == 1)/sum(test$Dengue == 1)
  
  
  specificity_ci <- binom.test(sum(preds <= thresh & 
                                     test$Dengue == 0), 
                               sum(test$Dengue == 0))$conf.int
  
  generalizability_results_dt$specificity_lower[r] <- specificity_ci[1]
  generalizability_results_dt$specificity_upper[r] <- specificity_ci[2]
  generalizability_results_dt$Specificity[r] <- sum(preds <= thresh & 
                                                      test$Dengue == 0)/sum(test$Dengue == 0)
  
  if(sum(preds > thresh) == 0){
    ppv <- NA
    generalizability_results_dt$ppv_lower[r] <- NA
    generalizability_results_dt$ppv_upper[r] <- NA
    generalizability_results_dt$PPV[r] <- NA
  } else {
    ppv_ci <- binom.test(sum(preds > thresh & 
                               test$Dengue == 1), 
                         sum(preds > thresh))$conf.int
    
    ppv <- sum(preds > thresh & test$Dengue == 1)/sum(preds > thresh)
    
    generalizability_results_dt$ppv_lower[r] <- ppv_ci[1]
    generalizability_results_dt$ppv_upper[r] <- ppv_ci[2]
    generalizability_results_dt$PPV[r] <- ppv
  }
  # ppv_ci <- binom.test(sum(preds > thresh & 
  #                            test$Dengue == 1), 
  #                      sum(preds > thresh))$conf.int
  
  # generalizability_results_dt$ppv_lower[r] <- ppv_ci[1]
  # generalizability_results_dt$ppv_upper[r] <- ppv_ci[2]
  # generalizability_results_dt$PPV[r] <- sum(preds > thresh & 
  #                                          test$Dengue == 1)/sum(preds > thresh)
  
  
  
  if(sum(preds <= thresh) == 0){
    generalizability_results_dt$npv_lower[r] <- NA
    generalizability_results_dt$npv_upper[r] <- NA
    generalizability_results_dt$NPV[r] <- NA
  } else {
    npv_ci <- binom.test(sum(preds <= thresh & 
                               test$Dengue == 0), 
                         sum(preds <= thresh))$conf.int
    
    generalizability_results_dt$npv_lower[r] <- npv_ci[1]
    generalizability_results_dt$npv_upper[r] <- npv_ci[2]
    generalizability_results_dt$NPV[r] <- sum(preds <= thresh & 
                                                test$Dengue == 0)/sum(preds <= thresh)
  }
  
  # npv_ci <- binom.test(sum(preds <= thresh & 
  #                            test$Dengue == 0), 
  #                      sum(preds <= thresh))$conf.int
  # 
  # generalizability_results_dt$npv_lower[r] <- npv_ci[1]
  # generalizability_results_dt$npv_upper[r] <- npv_ci[2]
  # generalizability_results_dt$NPV[r] <- sum(preds <= thresh & 
  #                                          test$Dengue == 0)/sum(preds <= thresh)
  
}

# don't want to train and test on different days for Dataset 4
generalizability_results_dt <- generalizability_results_dt %>% 
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))



p2 <- generalizability_results %>%
  mutate(method = "Logistic regression") %>%
  rbind(generalizability_results_svm %>%
          mutate(method = "SVM")) %>%
  rbind(generalizability_results_dt %>%
          mutate(method = "Decision tree")) %>%
  filter(threshold == "Training threshold") %>%
  ggplot(aes(x = train_data, y = AUC, color = method, shape = method)) + 
  geom_errorbar(aes(ymin=auc_lower, ymax=auc_upper), width=0.2) +
  geom_point() +
  facet_wrap(~test_data) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7)) +
  labs(x = "Training dataset", y = "AUC for each test dataset", color = "Method",
       shape = "Method",
       title = "With age range adjustment")


### Supplementary Figure 4: comparison of the three different models 
### (logistic regression, SVM, decision tree) on each pair of training and test
### datasets, with and without age range adjustments

pdf(file = "../images/supplementary_figure_4.pdf",
    width = 12, height = 8)

grid.arrange(p1, p2, ncol=1)

dev.off()
