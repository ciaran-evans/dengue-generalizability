library(tidyverse)
library(readxl)
library(latex2exp)
library(patchwork)
library(pROC)
library(xtable)
library(predtools)
library(arsenal)
library(gridExtra)

### functions for label shift estimation

# fixed model method adapted from Saerens et al.
label_shift_estimate_search <- function(a, b, test_f, train_prop){
  possible_vals <- seq(a, b, by = 0.001)
  diffs <- c()
  for(v in possible_vals){
    w1 <- v/train_prop
    w2 <- (1-v)/(1 - train_prop)
    corrected_preds <- w1*test_f/(w1*test_f + w2*(1 - test_f))
    diffs <- c(diffs, abs(mean(corrected_preds) - v))
  }
  return(possible_vals[which.min(diffs)])
}

# confusion matrix method adapted from Lipton et al.
label_shift_estimate_lipton <- function(train_y, train_f, test_f, thresh){
  train_props <- c(1 - mean(train_y), mean(train_y))
  train_conf <- table(train_f > thresh, train_y)/length(train_f)
  test_pred_props <- c(1-mean(test_f > thresh), 
                       mean(test_f > thresh))
  w = solve(train_conf) %*% test_pred_props
  muy = train_props * w
  return(muy[2])
}


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
  rbind(data05)


### Supplementary table 1: Summary statistics for each dataset

# specifying summary statistics to calculate for numeric variables
# (mean, sd, median, q1, q3, range)
controls <- tableby.control(
  test=FALSE, 
  total=FALSE,
  numeric.stats = c("meansd", "medianq1q3", "range"),
  stats.labels = list(
    meansd = "Mean (SD)",
    medianq1q3 = "Median",
    range = "Range"
  ),
  digits = 3
)

# customize the labels of the table
labels <- list(
  Dengue = "Dengue",
  Age = "Age",
  WBC = "White Blood Cell Count (WBC)",
  PLT = "Platelet Count (PLT)"
)


# calculating the summary statistics
sum_table <- full_data %>%
  dplyr::select(Dataset, Dengue, Age, WBC, PLT) %>%
  mutate(Dengue = ifelse(Dengue == 0, "Negative", "Positive")) %>%
  tableby(Dataset ~ .,
                   data = .,
                   control = controls)

# convert to LaTeX format
summary(sum_table,
        labelTranslations = labels,
        title = "Summary Statistic of All Datasets",
        digits = 1,
        text = NULL) %>%
  as.data.frame() %>%
  xtable() %>%
  print(include.rownames=F)



### Figure 1: Feature distributions for each dataset
### (Distributions of Age, WBC, and PLT broken down by
### dengue status and dataset)

p1 <- full_data %>%
  filter(Dengue == 1) %>%
  ggplot(aes(x = Age, color = Dataset)) +
  geom_density() +
  theme_bw() +
  labs(title = "Dengue patients",
       x = "Age (years)") +
  theme(plot.title = element_text(size = 11))

p2 <- full_data %>%
  filter(Dengue == 1) %>%
  ggplot(aes(x = WBC, color = Dataset)) +
  geom_density() +
  theme_bw() +
  labs(x = TeX("WBC ($\\times 10^3 / mm^3$)"))

p3 <- full_data %>%
  filter(Dengue == 1) %>%
  ggplot(aes(x = PLT, color = Dataset)) +
  geom_density() +
  theme_bw() +
  labs(x = TeX("PLT ($\\times 10^3 / mm^3$)"))


p4 <- full_data %>%
  filter(Dengue == 0) %>%
  ggplot(aes(x = Age, color = Dataset)) +
  geom_density() +
  theme_bw() +
  labs(title = "Non-dengue patients",
       x = "Age (years)") +
  theme(plot.title = element_text(size = 11))

p5 <- full_data %>%
  filter(Dengue == 0) %>%
  ggplot(aes(x = WBC, color = Dataset)) +
  geom_density() +
  theme_bw() +
  labs(x = TeX("WBC ($\\times 10^3 / mm^3$)"))

p6 <- full_data %>%
  filter(Dengue == 0) %>%
  ggplot(aes(x = PLT, color = Dataset)) +
  geom_density() +
  theme_bw() +
  labs(x = TeX("PLT ($\\times 10^3 / mm^3$)"))


pdf(file = "../images/figure_1.pdf",
    width = 12, height = 8)

(p1 + p2 + p3)/(p4 + p5 + p6) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dev.off()



### Table 1: in-sample performance for each dataset
### Summarize performance when model is trained and evaluated on the
### same dataset. Cross-validation is used to avoid training and testing
### on the same observations

# set a seed for reproducibility
set.seed(1)

dataset_names <- unique(full_data$Dataset)
#thresholds <- c(0.33, 0.45, 0.21, 0.62, 0.71, 0.43)
thresholds <- c(0.33, 0.45, 0.21, 0.59, 0.64, 0.43)

# vectors to store in-sample performance metrics for
# each dataset
sensitivity <- rep(NA, length(dataset_names))
specificity <- rep(NA, length(dataset_names))
ppv <- rep(NA, length(dataset_names))
npv <- rep(NA, length(dataset_names))
auc <- rep(NA, length(dataset_names))

# calculate in-sample performance metrics for each dataset
for(i in 1:length(dataset_names)){
  
  cur_data <- full_data %>%
    filter(Dataset == dataset_names[i])
  
  yhat <- rep(NA, nrow(cur_data))
  
  if(dataset_names[i] %in% c("Dataset 1", "Dataset 2")){
    # leave-one-site-out CV (multiple sites in these Datasets)
    
    sites <- unique(cur_data$Site)
    
    for(j in 1:length(sites)){
      train <- cur_data %>%
        filter(Site != sites[j])
      
      test <- cur_data %>%
        filter(Site == sites[j])
      
      logistic_mod <- glm(Dengue ~ Age + WBC + PLT, 
                          family = binomial,
                          data = train)
      
      yhat[cur_data$Site == sites[j]] <- 
        predict(logistic_mod, newdata = test, type = "response")
      
    }
  } else {
    # 10-fold CV
    
    folds <- sample(1:10, nrow(cur_data), replace=T)
    
    for(j in 1:10){
      train <- cur_data[folds != j,]
      test <- cur_data[folds == j,]
      
      logistic_mod <- glm(Dengue ~ Age + WBC + PLT, 
                          family = binomial,
                          data = train)
      
      yhat[folds == j] <- 
        predict(logistic_mod, newdata = test, type = "response")
    }
  }
  
  sensitivity[i] <- sum((yhat > thresholds[i]) & 
                          cur_data$Dengue == 1)/sum(cur_data$Dengue)
  specificity[i] <- sum((yhat <= thresholds[i]) & 
                          cur_data$Dengue == 0)/sum(1 - cur_data$Dengue)
  ppv[i] <- sum((yhat > thresholds[i]) & 
                  cur_data$Dengue == 1)/sum(yhat > thresholds[i])
  npv[i] <- sum((yhat <= thresholds[i]) & 
                  cur_data$Dengue == 0)/sum(yhat <= thresholds[i])
  
  auc_ci <- round(c(ci.auc(roc(cur_data$Dengue, yhat))), 2)
  auc[i] <- paste(auc_ci[2], " (", auc_ci[1], ", ", auc_ci[3], ")", sep="")
}

# store the in-sample results
in_sample_results <- data.frame(Dataset = dataset_names,
           Threshold = thresholds,
           Sensitivity = round(sensitivity, 3),
           Specificity = round(specificity, 3),
           PPV = round(ppv, 3),
           NPV = round(npv, 3),
           AUC = auc) 

# convert to LaTeX for paper
in_sample_results %>%
  xtable() %>%
  print(include.rownames = F)



### Results: Generalizability
### (Figure 2 and Supplementary table 2)
### Goal: assess model performance on new data for
### each training/test pair


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

### Figure 2: model performance for each pair of training/test datasets

p1 <- generalizability_results %>%
  select(train_data, test_data, AUC, auc_lower, auc_upper) %>%
  distinct() %>%
  ggplot(aes(x = train_data, y = AUC, group=1)) +
  geom_errorbar(aes(ymin=auc_lower, ymax=auc_upper), width=0.2) +
  geom_point() +
  facet_wrap(~test_data) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7)) +
  labs(x = "Training dataset", y = "AUC for each test dataset") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


p2 <- generalizability_results %>%
  select(train_data, test_data, Sensitivity, sensitivity_lower, 
         sensitivity_upper, threshold) %>%
  distinct() %>%
  ggplot(aes(x = train_data, y = Sensitivity, color = threshold)) +
  geom_point(position=position_dodge(width=0.75),
             aes(shape = threshold)) +
  geom_errorbar(aes(ymin=sensitivity_lower, ymax=sensitivity_upper), width=0.2,
                position=position_dodge(width=0.75)) +
  facet_wrap(~test_data) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7)) +
  labs(x = "Training dataset", y = "Sensitivity for each test dataset", color = "",
       shape = "") +
  theme(legend.position="bottom") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


p3 <- generalizability_results %>%
  select(train_data, test_data, Specificity, specificity_lower, 
         specificity_upper, threshold) %>%
  distinct() %>%
  ggplot(aes(x = train_data, y = Specificity, color = threshold, 
             shape=threshold)) +
  geom_point(position=position_dodge(width=0.75),
             aes(shape = threshold)) +
  geom_errorbar(aes(ymin=specificity_lower, ymax=specificity_upper), width=0.2,
                position=position_dodge(width=0.75)) +
  facet_wrap(~test_data) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 7)) +
  labs(x = "Training dataset", y = "Specificity for each test dataset", color = "",
       shape = "") +
  theme(legend.position="bottom")


pdf(file = "../images/figure_2.pdf",
    width = 12, height = 12)

p1/p2/p3 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

dev.off()


### Supplementary table 2: full performance metrics for each combination of
### training/test datasets

generalizability_results %>%
  distinct() %>%
  arrange(train_data, test_data) %>%
  filter(train_data != test_data | threshold == "Training threshold") %>%
  rename("Training data" = train_data,
         "Test data" = test_data,
         "Threshold" = threshold) %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  mutate(Sensitivity = paste(Sensitivity, " (",
                             sensitivity_lower, ", ",
                             sensitivity_upper, ")", sep=""),
         Specificity = paste(Specificity, " (",
                             specificity_lower, ", ",
                             specificity_upper, ")", sep=""),
         PPV = paste(PPV, " (",
                     ppv_lower, ", ",
                     ppv_upper, ")", sep=""),
         NPV = paste(NPV, " (",
                     npv_lower, ", ",
                     npv_upper, ")", sep=""),
         AUC = paste(AUC, " (",
                     auc_lower, ", ",
                     auc_upper, ")", sep="")) %>%
  select("Training data", "Test data", Threshold,
         Sensitivity, Specificity, PPV, NPV, AUC) %>%
  xtable() %>%
  print(include.rownames=F)


### Results: Restricting age ranges


training_sets <- c("Dataset 1", "Dataset 4, Day -3",
                   "Dataset 4, Day -1", "Dataset 5")
test_sets <- c("Dataset 2", "Dataset 3")
filtered <- c("Unrestricted", "Restricted")
#training_thresholds <- c(0.33, 0.62, 0.71, 0.43)
training_thresholds <- c(0.33, 0.59, 0.64, 0.43)

# combinations of training/test sets to further explore
# when restricting age range of the test set
age_filtered_results <- expand.grid(training_sets, test_sets,
                                    filtered) %>%
  rename(train_data = Var1, test_data = Var2, filtered = Var3) %>%
  mutate(Sensitivity = NA, Specificity = NA, AUC = NA)

for(r in 1:nrow(age_filtered_results)){
  # select the training and test data
  train <- full_data %>%
    filter(Dataset == age_filtered_results$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == age_filtered_results$test_data[r])
  
  # restrict age range on the test data
  # (> 16 if training data is Dataset 5, 
  # < 16 if training data is Dataset 1 or Dataset 4)
  if(age_filtered_results$filtered[r] == "Restricted" & 
     age_filtered_results$train_data[r] != "Dataset 5"){
    test <- test %>%
      filter(Age < 16)
  } else if(age_filtered_results$filtered[r] == "Restricted"){
    test <- test %>%
      filter(Age > 16)
  }
  
  # fit the logistic regression model on training data
  logistic_mod <- glm(Dengue ~ Age + WBC + PLT, 
                      family = binomial,
                      data = train)
  
  # predict on (possibly age-restricted) test data
  preds <- predict(logistic_mod, newdata = test, type = "response")
  
  # calculate performance metrics
  age_filtered_results$AUC[r] <- round(c(ci.auc(roc(test$Dengue, preds))), 3)[2]
  
  age_filtered_results$Sensitivity[r] <- sum(preds > training_thresholds[which(training_sets == 
                                                age_filtered_results$train_data[r])] & 
                       test$Dengue == 1)/sum(test$Dengue == 1)
  
  age_filtered_results$Specificity[r] <- sum(preds <= training_thresholds[which(training_sets == 
                                                age_filtered_results$train_data[r])] & 
                       test$Dengue == 0)/sum(test$Dengue == 0)
  
}


### Figure 3: Model performance with restricted age ranges

pdf(file = "../images/figure_3.pdf",
    width = 12, height = 6)

age_filtered_results %>%
  mutate(train_data = paste("Train:", train_data),
         test_data = paste("Test:", test_data)) %>%
  pivot_longer(cols = c(Sensitivity, Specificity, AUC),
               names_to = "metric",
               values_to = "performance") %>%
  ggplot(aes(x = filtered, y = performance, color = metric, group = metric)) +
  geom_point(aes(shape = metric), size=2) +
  geom_line() +
  facet_grid(test_data ~ train_data) +
  theme_bw() +
  labs(x = "Age range for test data",
       y = "Performance on test data",
       color = "", shape = "") +
  theme(legend.position = "bottom")

dev.off()


### Calibration plots
### (Figure 4, Supplementary figure 1, Supplementary figure 2)

### Figure 4: Example calibration plots for several training/test pairs

# train on Data 1, test on Data 3 (age restricted)

train <- full_data %>%
  filter(Dataset == "Dataset 1")

test <- full_data %>%
  filter(Dataset == "Dataset 3",
         Age < 16)

train_frac <- mean(train$Dengue)
test_frac <- mean(test$Dengue)

logistic_mod <- glm(Dengue ~ Age + WBC + PLT, data = train, family = binomial)

preds <- predict(logistic_mod, newdata = test, type = "response")

new_preds <- (test_frac/train_frac)*preds/((test_frac/train_frac)*preds + 
                                             (1-test_frac)/(1-train_frac)*(1-preds))

p1 <- (test %>%
  mutate(preds = preds) %>%
  as.data.frame() %>%
  calibration_plot(obs = "Dengue", pred = "preds",
                   y_lim = c(0, 1), x_lim = c(0, 1),
                   xlab = "Uncorrected predicted probability of dengue (training data)",
                   ylab = "True dengue rate (test data)",
                   title = "Train: Dataset 1, Test: Dataset 3 (Age < 16)"))$calibration_plot +
  theme(text = element_text(size = 9))

p2 <- (test %>%
  mutate(preds = new_preds) %>%
  as.data.frame() %>%
  calibration_plot(obs = "Dengue", pred = "preds",
                   y_lim = c(0, 1), x_lim = c(0, 1),
                   xlab = "Corrected predicted probability of dengue (training data)",
                   ylab = "True dengue rate (test data)"))$calibration_plot +
  theme(text = element_text(size = 9))


# train on Data 3, test on Data 5

train <- full_data %>%
  filter(Dataset == "Dataset 3")

test <- full_data %>%
  filter(Dataset == "Dataset 5")

train_frac <- mean(train$Dengue)
test_frac <- mean(test$Dengue)

logistic_mod <- glm(Dengue ~ Age + WBC + PLT, data = train, family = binomial)

preds <- predict(logistic_mod, newdata = test, type = "response")

new_preds <- (test_frac/train_frac)*preds/((test_frac/train_frac)*preds + 
                                             (1-test_frac)/(1-train_frac)*(1-preds))

p3 <- (test %>%
         mutate(preds = preds) %>%
         as.data.frame() %>%
         calibration_plot(obs = "Dengue", pred = "preds",
                          y_lim = c(0, 1), x_lim = c(0, 1),
                          xlab = "Uncorrected predicted probability of dengue (training data)",
                          ylab = "True dengue rate (test data)",
                          title = "Train: Dataset 3, Test: Dataset 5"))$calibration_plot +
  theme(text = element_text(size = 9))

p4 <- (test %>%
         mutate(preds = new_preds) %>%
         as.data.frame() %>%
         calibration_plot(obs = "Dengue", pred = "preds",
                          y_lim = c(0, 1), x_lim = c(0, 1),
                          xlab = "Corrected predicted probability of dengue (training data)",
                          ylab = "True dengue rate (test data)"))$calibration_plot +
  theme(text = element_text(size = 9))



# train on Data 4 Day -3, test on Data 1 

train <- full_data %>%
  filter(Dataset == "Dataset 4, Day -3")

test <- full_data %>%
  filter(Dataset == "Dataset 1")

train_frac <- mean(train$Dengue)
test_frac <- mean(test$Dengue)

logistic_mod <- glm(Dengue ~ Age + WBC + PLT, data = train, family = binomial)

preds <- predict(logistic_mod, newdata = test, type = "response")

new_preds <- (test_frac/train_frac)*preds/((test_frac/train_frac)*preds + 
                                             (1-test_frac)/(1-train_frac)*(1-preds))

p5 <- (test %>%
         mutate(preds = preds) %>%
         as.data.frame() %>%
         calibration_plot(obs = "Dengue", pred = "preds",
                          y_lim = c(0, 1), x_lim = c(0, 1),
                          xlab = "Uncorrected predicted probability of dengue (training data)",
                          ylab = "True dengue rate (test data)",
                          title = "Train: Dataset 4 Day -3, Test: Dataset 1"))$calibration_plot +
  theme(text = element_text(size = 9))

p6 <- (test %>%
         mutate(preds = new_preds) %>%
         as.data.frame() %>%
         calibration_plot(obs = "Dengue", pred = "preds",
                          y_lim = c(0, 1), x_lim = c(0, 1),
                          xlab = "Corrected predicted probability of dengue (training data)",
                          ylab = "True dengue rate (test data)"))$calibration_plot +
  theme(text = element_text(size = 9))

pdf(file = "../images/figure_4.pdf",
    width = 12, height = 7)

(p1 + p3 + p5)/(p2 + p4 + p6)

dev.off()




### Results: Estimating dengue prevalence in test data
### can we use the label shift assumption to estimate the 
### fraction of dengue patients?

dataset_names <- unique(full_data$Dataset)

# combinations of training/test pairs for estimating prevalence in new data
label_shift_adjustment <- expand.grid(dataset_names, dataset_names) %>%
  rename(train_data = Var1, test_data = Var2) %>%
  mutate(deviance_uncorrected = NA,
         deviance_corrected = NA,
         train_frac = NA,
         test_frac_true = NA,
         test_frac_uncorrected_est = NA,
         test_frac_search_est = NA,
         test_frac_lipton_est = NA)
label_shift_adjustment$test_data <- as.character(label_shift_adjustment$test_data)
label_shift_adjustment <- label_shift_adjustment %>%
  filter(!(train_data == "Dataset 4, Day -3" & test_data == "Dataset 4, Day -1"),
         !(train_data == "Dataset 4, Day -1" & test_data == "Dataset 4, Day -3"))

# also make calibration plots for each pair
calibration_plots_uncorrected <- list()
calibration_plots_corrected <- list()

for(r in 1:nrow(label_shift_adjustment)){
  # select training and test data
  train <- full_data %>%
    filter(Dataset == label_shift_adjustment$train_data[r])
  
  test <- full_data %>%
    filter(Dataset == label_shift_adjustment$test_data[r])
  
  # restrict ages for certain test sets
  if(label_shift_adjustment$train_data[r] %in% c("Dataset 1", "Dataset 4, Day -1",
                                                 "Dataset 4, Day -3")
     & label_shift_adjustment$test_data[r] %in% c("Dataset 2", "Dataset 3")){
    test <- test %>%
      filter(Age < 16)
    
    label_shift_adjustment$test_data[r] <- paste(label_shift_adjustment$test_data[r], 
                                                 "(Age < 16)")
  }
  
  # fit logistic regression model
  logistic_mod <- glm(Dengue ~ Age + WBC + PLT, 
                      family = binomial,
                      data = train)
  
  train_frac <- mean(train$Dengue)
  test_frac <- mean(test$Dengue)
  
  # uncorrected predictions on test data
  preds <- predict(logistic_mod, newdata = test, type = "response")
  
  # label shift-adjusted predictions on test data
  # (See Equation 1 in the manuscript)
  new_preds <- (test_frac/train_frac)*preds/((test_frac/train_frac)*preds + 
                                               (1-test_frac)/(1-train_frac)*(1-preds))
  
  # calculate average deviance before and after the label shift adjustment
  # (to summarize whether the label shift adjustment improves predictive performance)
  label_shift_adjustment$deviance_uncorrected[r] <- 
    -2*mean(test$Dengue*log(preds) + (1 - test$Dengue)*log(1 - preds))
  
  label_shift_adjustment$deviance_corrected[r] <- 
    -2*mean(test$Dengue*log(new_preds) + (1 - test$Dengue)*log(1 - new_preds))
  
  # calculate different estimates of the dengue prevalence in test data
  label_shift_adjustment$train_frac[r] <- train_frac
  label_shift_adjustment$test_frac_true[r] <- test_frac
  label_shift_adjustment$test_frac_uncorrected_est[r] <- mean(preds)
  label_shift_adjustment$test_frac_search_est[r] <- 
    label_shift_estimate_search(0.01, 0.99, preds, train_frac)
  label_shift_adjustment$test_frac_lipton_est[r] <- 
    label_shift_estimate_lipton(train$Dengue, 
                                predict(logistic_mod, type="response"),
                                preds, 0.5)
  
  # calibration plot for the uncorrected predictions
  calibration_plots_uncorrected[[r]] <- (test %>%
           mutate(preds = preds) %>%
           as.data.frame() %>%
           calibration_plot(obs = "Dengue", pred = "preds",
                            y_lim = c(0, 1), x_lim = c(0, 1),
                            xlab = "Uncorrected predicted probability of dengue (training data)",
                            ylab = "True dengue rate (test data)",
                            title = paste("Train: ",
                                          label_shift_adjustment$train_data[r],
                                          ", Test: ",
                                          label_shift_adjustment$test_data[r],
                                          sep = "")))$calibration_plot +
    theme(text = element_text(size = 9))
  
  # calibration plot for the corrected predictions
  calibration_plots_corrected[[r]] <- (test %>%
                                           mutate(preds = new_preds) %>%
                                           as.data.frame() %>%
                                           calibration_plot(obs = "Dengue", pred = "preds",
                                                            y_lim = c(0, 1), x_lim = c(0, 1),
                                                            xlab = "Corrected predicted probability of dengue (training data)",
                                                            ylab = "True dengue rate (test data)",
                                                            title = paste("Train: ",
                                                                          label_shift_adjustment$train_data[r],
                                                                          ", Test: ",
                                                                          label_shift_adjustment$test_data[r],
                                                                          sep = "")))$calibration_plot +
    theme(text = element_text(size = 9))
}


### Figure 5: Estimating prevalence of dengue in test data

pdf(file = "../images/figure_5.pdf",
    width = 12, height = 8)

label_shift_adjustment %>%
  pivot_longer(cols = c(test_frac_uncorrected_est, test_frac_search_est, 
                        test_frac_lipton_est),
               names_to = "method",
               values_to = "estimate") %>%
  mutate(method = fct_recode(method,  "Uncorrected estimate" = "test_frac_uncorrected_est",
                         "Fixed point estimate" = "test_frac_search_est",
                         "Discretization estimate" = "test_frac_lipton_est")) %>%
  ggplot(aes(x = train_data, y = estimate, color = method, shape = method)) +
  geom_point(position=position_dodge(width=0.5),
             size = 2) +
  geom_hline(aes(yintercept = test_frac_true)) +
  facet_wrap(~test_data) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, size = 7)) +
  labs(color = "", shape = "",
       x = "Training dataset", y = "Estimate of dengue prevalence in test data")

dev.off()


### Supplementary figure 1: uncorrected calibration plots
pdf(file = "../images/supplementary_figure_1.pdf",
    width = 26, height = 26)
do.call("grid.arrange", c(calibration_plots_uncorrected, ncol=6))
dev.off()

### Supplementary figure 2: corrected calibration plots
pdf(file = "../images/supplementary_figure_2.pdf",
    width = 26, height = 26)
do.call("grid.arrange", c(calibration_plots_corrected, ncol=6))
dev.off()


output <- label_shift_adjustment %>%
  arrange(train_data, test_data) %>%
  filter(train_data != test_data)

# 20/28 training/test pairs have lower deviance when we
# adjust the probabilities with the label shift correction
sum(output$deviance_corrected < output$deviance_uncorrected)

# summarize performance for each estimation method

output %>%
  summarize(fc_uncorrected = median(log(test_frac_uncorrected_est/test_frac_true)),
            fc_search = median(log(test_frac_search_est/test_frac_true)),
            fc_lipton = median(log(test_frac_lipton_est/test_frac_true), na.rm=T))

### Supplementary table 3: label shift adjustment and estimation for each
### training/test pair
output <- output %>%
  mutate_if(is.numeric, round, digits = 3) %>%
  rename("Training data" = train_data,
         "Test data" = test_data,
         "Deviance (uncorrected predictions)" = deviance_uncorrected,
         "Deviance (corrected predictions)" = deviance_corrected,
         "Training prevalence" = train_frac,
         "Test prevalence" = test_frac_true,
         "Est. test prevalence, average" = test_frac_uncorrected_est,
         "Est. test prevalence, discretization" = test_frac_lipton_est,
         "Est. test prevalence, fixed point" = test_frac_search_est) %>%
  xtable()

align(output) <- c( rep('l', 3),  rep('>{\\centering}p{1.5in}', 7))

output %>%
  print(include.rownames = F)



