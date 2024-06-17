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

# Dataset 1 (Tuan et al., 2015)
data01 <- read_excel("../data/set-01/dengue-data-01.xls") %>%
  dplyr::select(SiteNo, Age, WBC, PLT, Lab_Confirmed_Dengue) %>%
  rename(Site = SiteNo,
         Dengue = Lab_Confirmed_Dengue) %>%
  mutate(Dengue = 2 - Dengue) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 1")

# Dataset 2 (Gasem et al., 2020)
data02 <- read_excel("../data/set-02/dengue-data-02.xlsx") %>%
  dplyr::select(Site, `Age (Converted from Days to Years)`, LBLEUKRS, LBTHRRS, `FINAL Category`) %>%
  rename(Age = `Age (Converted from Days to Years)`,
         WBC = LBLEUKRS, PLT = LBTHRRS,
         Dengue = `FINAL Category`) %>%
  mutate(Dengue = ifelse(Dengue == "Dengue", 1, 0)) %>%
  drop_na() %>%
  mutate(Dataset = "Dataset 2")

# Dataset 3 (Saito et al., 2022)
data03 <- read_excel("../data/set-03/dengue-data-03.xlsx") %>%
  dplyr::select(age, wbc, plt, dengue) %>%
  drop_na() %>%
  rename(Age = age, WBC = wbc, PLT = plt, Dengue = dengue) %>%
  mutate(Site = 1, Dataset = "Dataset 3")

# Dataset 4 (Park et al., 2018)
data04 <- read_excel("../data/set-04/dengue-data-04.xlsx") %>%
  dplyr::select(age2, wbc_m3, wbc_m1, platelets_m3, platelets_m1, dengue) %>%
  pivot_longer(c(wbc_m3, wbc_m1, platelets_m3, platelets_m1),
               names_to = c(".value", "Day"),
               names_sep = "_m") %>%
  rename(Dengue = dengue, WBC = wbc, PLT = platelets) %>%
  separate(age2, c("age1", "age2")) %>%
  mutate(Age = (as.numeric(age1) + as.numeric(age2))/2,
         WBC = WBC/1000, PLT = PLT/1000,
         Dataset = paste("Dataset 4, Day -", Day, sep=""),
         Site = 1) %>%
  dplyr::select(Site, Age, WBC, PLT, Dengue, Dataset)

# Dataset 5 (Ngim et al., 2021)
data05 <- read_excel("../data/set-05/dengue-data-05.xls") %>%
  dplyr::select(AgeEnrol, Result_WBC, Result_platelet, ConfirmDengue_YesNo_FourCriteria) %>%
  rename(Age = AgeEnrol, WBC = Result_WBC, PLT = Result_platelet, 
         Dengue = ConfirmDengue_YesNo_FourCriteria) %>%
  mutate(Dengue = ifelse(Dengue == "No dengue", 0, 1),
         Site = 1,
         Dataset = "Dataset 5") 

# combine the five datasets into one
full_data <- data01 %>%
  rbind(data02) %>%
  rbind(data03) %>%
  rbind(data04) %>%
  rbind(data05) %>%
  mutate(Dengue = as.factor(Dengue))



##### Part 1: without adjusting age ranges


# Logistic regression

dataset_names <- unique(full_data$Dataset)
#thresholds <- c(0.33, 0.45, 0.21, 0.62, 0.71, 0.43)
thresholds <- c(0.33, 0.45, 0.21, 0.59, 0.64, 0.43)

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
  logistic_mod <- glm(Dengue ~ Age + WBC + PLT, 
                      family = binomial,
                      data = train)
  
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
  
  
  ppv_ci <- binom.test(sum(preds > thresh & 
                             test$Dengue == 1), 
                       sum(preds > thresh))$conf.int
  
  generalizability_results$ppv_lower[r] <- ppv_ci[1]
  generalizability_results$ppv_upper[r] <- ppv_ci[2]
  generalizability_results$PPV[r] <- sum(preds > thresh & 
                                           test$Dengue == 1)/sum(preds > thresh)
  
  npv_ci <- binom.test(sum(preds <= thresh & 
                             test$Dengue == 0), 
                       sum(preds <= thresh))$conf.int
  
  generalizability_results$npv_lower[r] <- npv_ci[1]
  generalizability_results$npv_upper[r] <- npv_ci[2]
  generalizability_results$NPV[r] <- sum(preds <= thresh & 
                                           test$Dengue == 0)/sum(preds <= thresh)
  
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
  
  svm_mod <- svm(Dengue ~ Age + WBC + PLT, data = train, kernel = "radial",
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
  
  dt_mod <- rpart(Dengue ~ Age + WBC + PLT, data = train, method = "class")
  
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
  logistic_mod <- glm(Dengue ~ Age + WBC + PLT, 
                      family = binomial,
                      data = train)
  
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
  
  
  ppv_ci <- binom.test(sum(preds > thresh & 
                             test$Dengue == 1), 
                       sum(preds > thresh))$conf.int
  
  generalizability_results$ppv_lower[r] <- ppv_ci[1]
  generalizability_results$ppv_upper[r] <- ppv_ci[2]
  generalizability_results$PPV[r] <- sum(preds > thresh & 
                                           test$Dengue == 1)/sum(preds > thresh)
  
  npv_ci <- binom.test(sum(preds <= thresh & 
                             test$Dengue == 0), 
                       sum(preds <= thresh))$conf.int
  
  generalizability_results$npv_lower[r] <- npv_ci[1]
  generalizability_results$npv_upper[r] <- npv_ci[2]
  generalizability_results$NPV[r] <- sum(preds <= thresh & 
                                           test$Dengue == 0)/sum(preds <= thresh)
  
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
  
  svm_mod <- svm(Dengue ~ Age + WBC + PLT, data = train, kernel = "radial",
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
  
  dt_mod <- rpart(Dengue ~ Age + WBC + PLT, data = train, method = "class")
  
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


### Supplementary Figure 3: comparison of the three different models 
### (logistic regression, SVM, decision tree) on each pair of training and test
### datasets, with and without age range adjustments

pdf(file = "../images/supplementary_figure_3.pdf",
    width = 12, height = 12)

grid.arrange(p1, p2, ncol=1)

dev.off()
