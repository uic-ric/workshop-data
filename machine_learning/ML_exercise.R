# =============================================================================
# Machine Learning
# Research Informatics Core
# Fall 2026
#
# Generated from ML_exercise.Rmd by purl_handouts.R -- do not edit by hand.
# =============================================================================

# =============================================================================
# MORNING
# =============================================================================

# -----------------------------------------------------------------------------
# Create a dataset for machine learning
# -----------------------------------------------------------------------------

# --- Select cancer dataset for download (INSTRUCTOR DEMONSTRATION) ----------

library(TCGAretriever)

# Obtain a list of cancer studies from cBio
all_studies <- get_cancer_studies()

# Get list of breast cancer studies
brca_studies <- all_studies[all_studies$cancerTypeId == "brca",]
nrow(brca_studies)

# Find studies with > 1000 patients
brca_studies[brca_studies$allSampleCount > 1000,c("name", "allSampleCount", "studyId")]

# Define the cancer study id: brca_tcga
my_csid <- "brca_tcga"

# Obtain genetic profiles, type of assays and data available
brca_pro <- get_genetic_profiles(csid = my_csid)
brca_pro$molecularProfileId

# Obtain cases
brca_cas <- get_case_lists(csid = my_csid)
brca_cas$sampleListId

# Determine cases with mRNA data available
q_csid <- "brca_tcga"
q_cases <- "brca_tcga_rna_seq_v2_mrna"
rna_prf <- "brca_tcga_rna_seq_v2_mrna"

# Download clinical data
brca_cli <- get_clinical_data(csid = q_csid, case_list_id = q_cases)
head(colnames(brca_cli))

# Count cancer types
table(brca_cli$CANCER_TYPE_DETAILED)
# Download mRNA expression as TPM
all_brca_RNA <- fetch_all_tcgadata(case_list_id = q_cases, 
                                   gprofile_id = rna_prf, 
                                   mutations = FALSE)

# --- Normalization and filtering --------------------------------------------

# Read mRNA expression data
all_brca_RNA <- readRDS(url("https://wd.cri.uic.edu/machine_learning/brca_RNA_TPM.rds"))
all_brca_RNA[1:5,1:5]
# use the gene symbol as the rownames, remove the other feature metadata
rownames(all_brca_RNA) <- all_brca_RNA$hugoGeneSymbol
all_brca_RNA <- all_brca_RNA[,-1:-3]
all_brca_RNA[1:5,1:5]
all_brca_RNA.log2 <- log2(all_brca_RNA + 1)
# Notes on this command:
#               1. Create boolean matrix based on values bigger than 0
#       2. Sum the rows: gets total number of samples with >0 expression
#                                       3. Compare to 30% of the number
#                                          of samples (total columns)
keep <- rowSums(all_brca_RNA.log2>0) >= 0.3 * ncol(all_brca_RNA.log2)
table(keep)
all_brca_RNA.log2.filtered <- all_brca_RNA.log2[keep,]

# --- Make a balanced dataset ------------------------------------------------

library(tidyverse)
# Read clinical data
brca_cli <- readRDS(url("https://wd.cri.uic.edu/machine_learning/brca_clinical.rds"))
# Remove samples that are missing CANCER_TYPE_DETAILED metadata
brca_cli <- subset(brca_cli, !is.na(CANCER_TYPE_DETAILED))
brca_cli.ductal <- brca_cli[
  brca_cli$CANCER_TYPE_DETAILED=="Breast Invasive Ductal Carcinoma",]
brca_cli.lobular <- brca_cli[
  brca_cli$CANCER_TYPE_DETAILED=="Breast Invasive Lobular Carcinoma",]
num_ductal <- nrow(brca_cli.ductal)
num_ductal
num_lobular <- nrow(brca_cli.lobular)
num_lobular
set.seed(123)
ductal.samples <- sample(brca_cli.ductal$sampleId, num_lobular)
ductal.exp <- all_brca_RNA.log2.filtered[,ductal.samples]
lobular.samples <- brca_cli.lobular$sampleId
lobular.exp <- all_brca_RNA.log2.filtered[,lobular.samples]

# Combine ductal and lobular expression datasets
expr_combined <- cbind(ductal.exp, lobular.exp)
num_ductal <- num_lobular
labels <- c(rep("Ductal", num_ductal), rep("Lobular", num_lobular))
table(labels)

# Build data.frame for modeling: samples x genes
df <- as.data.frame(t(expr_combined))  # rows = samples, cols = genes
df <- tibble::rownames_to_column(df, var = "sample_id")
df$class <- factor(labels, levels = c("Ductal","Lobular"))

# Fix gene names to a format acceptable by R
colnames(df) <- make.names(colnames(df))

# --- Split into training and testing cohorts --------------------------------

df <- readRDS(url("https://wd.cri.uic.edu/machine_learning/ML_df.rds"))
library(tidyverse)
library(caret)
set.seed(123)
train_idx <- createDataPartition(df$class, p = 0.7, list = FALSE)
train_df  <- df[train_idx, ]
test_df   <- df[-train_idx, ]

# For convenience, create x (data) and y (class labels)
x_train <- train_df %>% select(-sample_id, -class)
y_train <- train_df$class
x_test  <- test_df  %>% select(-sample_id, -class)
y_test  <- test_df$class

# -----------------------------------------------------------------------------
# Feature selection
# -----------------------------------------------------------------------------

load(url("https://wd.cri.uic.edu/machine_learning/MLobjects.RData"))

# --- Variance filtering of features -----------------------------------------

# Compute variance per gene
gene_vars <- sort(apply(x_train, 2, var), decreasing=TRUE)
head(gene_vars)
# Keep top 500 variable genes (adjust as needed)
top_var_genes <- names(gene_vars)[1:500]
df_var_filtered <- x_train[, top_var_genes]
df_var_filtered[1:5,1:5]

# --- T-test filtering of features -------------------------------------------

# T-test to rank remaining features
ttest_ranks <- apply(x_train, 2, function(g) {
  t.test(g ~ y_train)$statistic
})
ttest_ranks <- sort(abs(ttest_ranks), decreasing=TRUE)
head(ttest_ranks)
# Pick the top 100 genes
top_ttest_genes <- names(ttest_ranks)[1:100]
df_ttest_filtered <- x_train[, top_ttest_genes]
df_ttest_filtered[1:5,1:5]

# -----------------------------------------------------------------------------
# Logistic regression model with RFE
# -----------------------------------------------------------------------------

library(caret)
library(ggplot2)

# RFE setup
rfe_ctrl <- rfeControl(functions = lrFuncs, 
                   method = "cv", # cross validation
                   number = 10) # number of folds

# Number of features to be retained in the updated model
sizes <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

# Run RFE
set.seed(123)
rfe_results <- rfe(
  x = df_ttest_filtered,
  y = y_train,
  sizes = sizes,
  rfeControl = rfe_ctrl
)

rfe_results

# Get best features
best_rfe_feats <- predictors(rfe_results)
best_rfe_feats

# Plot RFE results
rfe_results_df <- rfe_results$results

ggplot(rfe_results_df, aes(x = Variables, y = Accuracy)) +
  geom_line() +
  geom_point() +
  labs(y = "Accuracy", x = "Number of genes", title = "RFE logistic regression") +
  theme_classic()

# Make training and testing datasets using best features
Xtr <- x_train[, best_rfe_feats]
Xte <- x_test[, best_rfe_feats]

# -----------------------------------------------------------------------------
# Test cohort and cross-validation with LR, confusion matrix, ROC/AUC
# -----------------------------------------------------------------------------

load(url("https://wd.cri.uic.edu/machine_learning/MLobjects.RData"))
library(caret)
library(ggplot2)
library(pROC)

# Define cross-validation
lr_ctrl <- trainControl(
  method = "cv",
  number = 10,              # 10-fold CV
  classProbs = TRUE,
  summaryFunction = twoClassSummary, # Report ROC, Sensitivity, and Specificity
  savePredictions = "final"  # keep CV predictions
)

set.seed(123)
# Train logistic regression
lr_model <- train(
  x = Xtr,
  y = y_train,
  method = "glm",
  family = "binomial",
  trControl = lr_ctrl,
  metric = "ROC"
)
lr_model # CV results summary
# Cross-validation predictions
lr_train_preds <- lr_model$pred$pred
lr_train_obs <- lr_model$pred$obs

# Make prediction table of training data
lr_train_pred_table <- table(Predicted = lr_train_preds, Actual = lr_train_obs)
lr_train_pred_table

# Make confusion matrix for model on training data cross-validation
confusionMatrix(lr_train_preds, lr_train_obs, positive = "Lobular")

# Make test set predictions
lr_test_preds <- predict(lr_model, newdata = Xte)
head(lr_test_preds)

# Make prediction table of testing data
lr_test_pred_table <- table(Predicted = lr_test_preds, Actual = y_test)
lr_test_pred_table

# Make confusion matrix for model on testing data
confusionMatrix(lr_test_preds, y_test, positive = "Lobular")

# Calculate the ROC curve
prob_pred <- predict(lr_model, newdata = Xte, type = "prob")$Lobular
head(prob_pred)
lr_roc <- roc(response = y_test, predictor = prob_pred)
lr_auc <- auc(lr_roc)
lr_auc

# Plotting the ROC curve with ggroc
ggroc(lr_roc) +
  geom_abline(intercept=1, slope=1, linetype="dashed") +
  labs(title="ROC Logistic Regression: Ductal vs Lobular",
    subtitle=paste0("AUC = ", round(lr_auc, 3)))

# =============================================================================
# AFTERNOON
# =============================================================================

# -----------------------------------------------------------------------------
# Train an SVM model using the same predictors
# -----------------------------------------------------------------------------

library(caret)
library(ggplot2)
library(pROC)
load(url("https://wd.cri.uic.edu/machine_learning/MLobjects.RData"))
library(kernlab)
svm_ctrl <- trainControl(
  method = "cv",   # cross-validation
  number = 10,     # 10-fold CV
  classProbs = TRUE,
  summaryFunction = twoClassSummary, # Report ROC, Sensitivity, and Specificity
  savePredictions = "final"  # keep CV predictions
)

set.seed(123)
svm_model <- train(
  y_train ~ .,
  data = data.frame(y_train, Xtr),
  method = "svmLinear",
  trControl = svm_ctrl,
  tuneGrid = expand.grid(C = c(0.1, 1, 10)),
  metric = "ROC"
)
svm_model
# Test predictions
svm_probs <- predict(svm_model, newdata = Xte, type = "prob")
head(svm_probs)
svm_preds <- predict(svm_model, newdata = Xte)
head(svm_preds)
table(svm_preds)
accuracy <- mean(svm_preds == y_test)
accuracy
# ROC/AUC
svm_roc <- roc(response = y_test, predictor = svm_probs$Lobular)
svm_auc <- auc(svm_roc)
svm_auc
# Plotting the ROC curve with ggroc
ggroc(svm_roc) +
  geom_abline(intercept=1, slope=1, linetype="dashed")

# --- Comparison of two models -----------------------------------------------

# recall the performance of each model from cross validation
lr_model$results$ROC
max(svm_model$results$ROC)
# use the probabilities for each model
logistic_lob <- prob_pred
svm_lob <- svm_probs[,"Lobular"]
head(logistic_lob)
head(svm_lob)
# set paired=T for the signed rank test
wilcox.test(logistic_lob, svm_lob, paired=T)
# compare in a scatterplot
plot_df <- data.frame(Logistic=logistic_lob, SVM=svm_lob)
# basic R scatterplot
plot(plot_df)

# -----------------------------------------------------------------------------
# Train a neural net model with regularization
# -----------------------------------------------------------------------------

library(caret)
library(ggplot2)
library(pROC)
load(url("https://wd.cri.uic.edu/machine_learning/MLobjects.RData"))
nn_ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(123)
nn_model <- train(
  y_train ~ .,
  data = data.frame(y_train, Xtr),
  method = "nnet",
  trControl = nn_ctrl,
  metric = "ROC",
  tuneGrid = expand.grid(
    size = c(3, 5, 7),   # hidden layer sizes to try
    decay = c(0.001, 0.01, 0.1)  # regularization values
  ),
  trace = FALSE,
  maxit = 500
)
nn_model
# Predict on test set
nn_probs <- predict(nn_model, newdata = Xte, type = "prob")
head(nn_probs)
nn_preds <- predict(nn_model, newdata = Xte)
head(nn_preds)
table(nn_preds)
nn_accuracy <- mean(nn_preds == y_test)
nn_accuracy
# ROC/AUC
nn_roc <- roc(response = y_test, predictor = nn_probs$Lobular)
nn_auc <- auc(nn_roc)
nn_auc
# Plotting the ROC curve with ggroc
ggroc(nn_roc) +
  geom_abline(intercept=1, slope=1, linetype="dashed")

# -----------------------------------------------------------------------------
# Implement a random forest and a XGBoost model
# -----------------------------------------------------------------------------

library(caret)
library(ggplot2)
library(pROC)
load(url("https://wd.cri.uic.edu/machine_learning/MLobjects.RData"))

# --- Random forest ----------------------------------------------------------

rf_ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
  )

set.seed(123)
rf_model <- train(
  y_train ~ ., 
  data = data.frame(y_train, Xtr),
  method = "rf",
  trControl = rf_ctrl,
  metric = "ROC",
  tuneGrid = expand.grid(mtry = c(2,5))  # number of variables randomly sampled at each split
)
rf_model
# Predict on test set
rf_probs <- predict(rf_model, newdata = Xte, type = "prob")
head(rf_probs)
rf_preds <- predict(rf_model, newdata = Xte)
head(rf_preds)
table(rf_preds)
rf_accuracy <- mean(rf_preds == y_test)
rf_accuracy
# ROC/AUC
rf_roc <- roc(response = y_test, predictor = rf_probs$Lobular)
rf_auc <- auc(rf_roc)
rf_auc
# Plotting the ROC curve with ggroc
ggroc(rf_roc) +
  geom_abline(intercept=1, slope=1, linetype="dashed")

# --- XGBoost model ----------------------------------------------------------

library(caret)
library(pROC)
xgb_ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  verboseIter = FALSE
)

set.seed(123)
xgb_model <- train(
  y_train ~ ., 
  data = data.frame(y_train, Xtr),
  method = "xgbTree",
  trControl = xgb_ctrl,
  metric = "ROC",
  tuneGrid = expand.grid(
    nrounds = c(50, 100),    # The total number of trees built sequentially
    max_depth = c(2, 3, 4), # maximum depth of a tree
    eta = c(0.05, 0.1, 0.3), # learning rate, shrinkage
    gamma = 0, # pruning parameter: 0 -> tree keeps splitting as long as it reduces error
    colsample_bytree = 1, # Fraction of features randomly chosen to build each tree
    min_child_weight = 1, # Controls the minimum amount of data needed in a leaf node
    subsample = 1 # Fraction of training data sampled for growing each tree
  ),
  verbosity = 0     # hides xgboost internal warnings
)
xgb_model
xgb_probs <- predict(xgb_model, newdata = Xte, type = "prob")
head(xgb_probs)
xgb_preds <- predict(xgb_model, newdata = Xte)
head(xgb_preds)
table(xgb_preds)
xgb_accuracy <- mean(xgb_preds == y_test)
xgb_accuracy
# ROC/AUC
xgb_roc <- roc(response = y_test, predictor = xgb_probs$Lobular)
xgb_auc <- auc(xgb_roc)
xgb_auc
# Plotting the ROC curve with ggroc
ggroc(xgb_roc) +
  geom_abline(intercept=1, slope=1, linetype="dashed")

# -----------------------------------------------------------------------------
# Ridge regression model for Age
# -----------------------------------------------------------------------------

# --- Download DNA methylation data (INSTRUCTOR DEMONSTRATION) ---------------

library(TCGAretriever)

# The molecular profiles for the study include a methylation profile
brca_pro <- get_genetic_profiles(csid = "brca_tcga")
brca_pro[brca_pro$molecularAlterationType == "METHYLATION", "molecularProfileId"]

# Case list of samples with HM450 methylation data
brca_cas <- get_case_lists(csid = "brca_tcga")
brca_cas[grepl("HM450", brca_cas$name), c("sampleListId", "description")]

# Download methylation beta values for all genes (this takes a long time)
all_brca_meth <- fetch_all_tcgadata(case_list_id = "brca_tcga_methylation_hm450",
                                    gprofile_id = "brca_tcga_methylation_hm450",
                                    mutations = FALSE)

# --- Read in the methylation data -------------------------------------------

all_brca_meth <- readRDS(url(
  "https://wd.cri.uic.edu/machine_learning/brca_methylation_hm450.rds"))
brca_cli <- readRDS(url("https://wd.cri.uic.edu/machine_learning/brca_clinical.rds"))
all_brca_meth[1:5,1:5]
dim(all_brca_meth)
rownames(all_brca_meth) <- all_brca_meth$hugoGeneSymbol
all_brca_meth <- as.matrix(all_brca_meth[,-1:-3])
all_brca_meth[1:5,1:5]

# --- Preprocess the methylation data ----------------------------------------

keep <- rowSums(is.na(all_brca_meth)) == 0
table(keep)
all_brca_meth <- all_brca_meth[keep,]
range(all_brca_meth)
all_brca_meth.M <- log2(all_brca_meth / (1 - all_brca_meth))
range(all_brca_meth.M)
library(glmnet)
library(caret)
library(ggplot2)
meth.data <- t(all_brca_meth.M)
colnames(meth.data) <- make.names(colnames(meth.data))
# look up the age for each methylation sample
meth_age <- brca_cli$AGE[match(rownames(meth.data), brca_cli$sampleId)]
keep <- !is.na(meth_age)
table(keep)
meth_age <- meth_age[keep]
meth_age.data <- meth.data[keep, ]
dim(meth_age.data)
meth.gene_vars <- apply(meth_age.data, 2, var)
meth.top_var_genes <- names(sort(meth.gene_vars, decreasing = TRUE))[1:5000]
meth.df.var_filt <- meth_age.data[, meth.top_var_genes]
# make training and test sets
set.seed(123)
age_idx  <- createDataPartition(meth_age, p = 0.7, list = F)
age_x_train  <- meth.df.var_filt[age_idx, ]
age_y_train <- meth_age[age_idx]
age_x_test   <- meth.df.var_filt[-age_idx, ]
age_y_test <- meth_age[-age_idx]
set.seed(123)
cvfit <- cv.glmnet(age_x_train, age_y_train,
                   family = "gaussian", # regression with MSE
                   alpha = 0,   # Ridge,
                   nfolds = 5)  # 5-fold CV
best_lambda <- cvfit$lambda.min
best_lambda

lin_model <- glmnet(age_x_train, age_y_train, alpha = 0, lambda = best_lambda)
# Coefficients
head(lin_model$beta)
cv_mse <- cvfit$cvm[cvfit$lambda == best_lambda]
r2_cv <- 1 - cv_mse / var(age_y_train)
r2_cv
# Predictions
preds_lin <- predict(lin_model, newx = age_x_test)

# Performance (R², RMSE)
r2_lin <- cor(preds_lin, age_y_test)^2
r2_lin
rmse_lin <- sqrt(mean((preds_lin - age_y_test)^2))
rmse_lin

# Make a plot of regression
df_pred_ridge <- data.frame(
  Observed = age_y_test,
  Predicted = as.numeric(preds_lin)
)

# Scatter plot with 1:1 line
metrics_ridge <- paste0("R2 = ", round(r2_lin, 3), 
                       " RMSE = ", round(rmse_lin, 2))
ggplot(df_pred_ridge, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title = "Ridge Regression: Predicted vs Observed Age",
    subtitle = metrics_ridge,
    x = "Observed Age",
    y = "Predicted Age")

# --- Going further: probe-level methylation and normal tissue ---------------

# Obtaining and filtering the probe-level data (INSTRUCTOR DEMONSTRATION)

library(data.table)
library(caret)

# Download the probe-level HM450 matrix for TCGA BRCA from UCSC Xena (~780 MB)
xena_url <- "https://tcga.xenahubs.net/download/TCGA.BRCA.sampleMap/HumanMethylation450.gz"
download.file(xena_url, "brca_hm450_probes.tsv.gz")

# Clinical data (the same file used throughout the workshop)
brca_cli <- readRDS(url("https://wd.cri.uic.edu/machine_learning/brca_clinical.rds"))

# Read only the header to get the sample IDs
hm450_samples <- colnames(fread("brca_hm450_probes.tsv.gz", nrows = 0))[-1]

# Sample type from the last two characters, patient from the first 12
tissue <- ifelse(grepl("-01$", hm450_samples), "Tumor",
          ifelse(grepl("-11$", hm450_samples), "Normal", NA))
patient <- substr(hm450_samples, 1, 12)
age <- brca_cli$AGE[match(patient, substr(brca_cli$sampleId, 1, 12))]

# Keep tumor and normal samples from patients with a recorded age
keep_samples <- !is.na(tissue) & !is.na(age)
hm450_samples <- hm450_samples[keep_samples]
sample_table <- data.frame(sampleId = hm450_samples,
                           patientId = patient[keep_samples],
                           tissue = tissue[keep_samples],
                           AGE = age[keep_samples])
table(sample_table$tissue)

# Read the matrix for those samples only: one row per probe, one column per sample
hm450 <- fread("brca_hm450_probes.tsv.gz", select = c("sample", hm450_samples))
beta <- as.matrix(hm450[, -1])
rownames(beta) <- hm450$sample
rm(hm450)
dim(beta)
keep_probes <- rowSums(is.na(beta)) == 0 & grepl("^cg", rownames(beta))
table(keep_probes)
beta <- beta[keep_probes, ]
# 1. Top 10,000 variable probes, all samples
probe_vars <- apply(beta, 1, var)
top_var_probes <- names(sort(probe_vars, decreasing = TRUE))[1:10000]

# 2. Top 2,000 age-correlated probes in the tumor training split only
is_tumor <- sample_table$tissue == "Tumor"
tumor_age <- sample_table$AGE[is_tumor]
set.seed(123)
tumor_idx <- createDataPartition(tumor_age, p = 0.7, list = FALSE)
tumor_train_beta <- beta[, is_tumor][, tumor_idx]
probe_age_cor <- cor(t(tumor_train_beta), tumor_age[tumor_idx])[, 1]
top_cor_probes <- names(sort(abs(probe_age_cor), decreasing = TRUE))[1:2000]

# 3. Horvath 2013 clock probes and coefficients (Genome Biology 14:R115, Additional file 3)
horvath_url <- paste0("https://static-content.springer.com/esm/",
  "art%3A10.1186%2Fgb-2013-14-10-r115/MediaObjects/13059_2013_3156_MOESM3_ESM.csv")
horvath <- read.csv(horvath_url, skip = 2)
horvath_intercept <- horvath$CoefficientTraining[horvath$CpGmarker == "(Intercept)"]
horvath <- horvath[grepl("^cg", horvath$CpGmarker), c("CpGmarker", "CoefficientTraining")]
names(horvath) <- c("probe", "coefficient")
# 47 of the 353 clock probes have missing values in this data set and were removed above
horvath <- horvath[horvath$probe %in% rownames(beta), ]
nrow(horvath)

# Combine the three sets
selected_probes <- unique(c(top_var_probes, top_cor_probes, horvath$probe))
length(selected_probes)
probe_table <- data.frame(probe = selected_probes,
                          top_variable = selected_probes %in% top_var_probes,
                          tumor_age_correlated = selected_probes %in% top_cor_probes,
                          horvath_clock = selected_probes %in% horvath$probe)
meth_probes <- list(
  beta = beta[selected_probes, ],
  samples = sample_table,
  probes = probe_table,
  horvath = list(intercept = horvath_intercept, coefficients = horvath)
)
saveRDS(meth_probes, "brca_methylation_hm450_probes.rds")

# Read in the filtered probe-level data

meth_probes <- readRDS(url(
  "https://wd.cri.uic.edu/machine_learning/brca_methylation_hm450_probes.rds"))
names(meth_probes)
dim(meth_probes$beta)
meth_probes$beta[1:5,1:4]
head(meth_probes$samples)
table(meth_probes$samples$tissue)
head(meth_probes$probes)
probe.data <- t(meth_probes$beta)
probe_meta <- meth_probes$samples
dim(probe.data)

# Ridge regression in tumor samples

is_tumor <- probe_meta$tissue == "Tumor"
tumor.data <- probe.data[is_tumor, ]
tumor_age <- probe_meta$AGE[is_tumor]

set.seed(123)
tumor_idx <- createDataPartition(tumor_age, p = 0.7, list = F)
tumor_x_train <- tumor.data[tumor_idx, ]
tumor_y_train <- tumor_age[tumor_idx]
tumor_x_test <- tumor.data[-tumor_idx, ]
tumor_y_test <- tumor_age[-tumor_idx]

set.seed(123)
tumor_cvfit <- cv.glmnet(tumor_x_train, tumor_y_train, alpha = 0, nfolds = 5)
tumor_model <- glmnet(tumor_x_train, tumor_y_train, alpha = 0,
                      lambda = tumor_cvfit$lambda.min)
preds_tumor <- predict(tumor_model, newx = tumor_x_test)
r2_tumor <- cor(preds_tumor, tumor_y_test)^2
r2_tumor
rmse_tumor <- sqrt(mean((preds_tumor - tumor_y_test)^2))
rmse_tumor

# Ridge regression in adjacent normal tissue

is_normal <- probe_meta$tissue == "Normal"
normal.data <- probe.data[is_normal, ]
normal_age <- probe_meta$AGE[is_normal]
length(normal_age)

set.seed(123)
normal_idx <- createDataPartition(normal_age, p = 0.7, list = F)
normal_x_train <- normal.data[normal_idx, ]
normal_y_train <- normal_age[normal_idx]
normal_x_test <- normal.data[-normal_idx, ]
normal_y_test <- normal_age[-normal_idx]

set.seed(123)
normal_cvfit <- cv.glmnet(normal_x_train, normal_y_train, alpha = 0, nfolds = 5)
normal_model <- glmnet(normal_x_train, normal_y_train, alpha = 0,
                       lambda = normal_cvfit$lambda.min)
preds_normal <- predict(normal_model, newx = normal_x_test)
r2_normal <- cor(preds_normal, normal_y_test)^2
r2_normal
rmse_normal <- sqrt(mean((preds_normal - normal_y_test)^2))
rmse_normal

metrics_tumor <- paste0("Tumor: R2 = ", round(r2_tumor, 3),
                        ", RMSE = ", round(rmse_tumor, 2))
metrics_normal <- paste0("Normal: R2 = ", round(r2_normal, 3),
                         ", RMSE = ", round(rmse_normal, 2))
df_pred_probes <- rbind(
  data.frame(Observed = tumor_y_test, Predicted = as.numeric(preds_tumor),
             Tissue = metrics_tumor),
  data.frame(Observed = normal_y_test, Predicted = as.numeric(preds_normal),
             Tissue = metrics_normal)
)
# Order the panels: tumor first, then normal
df_pred_probes$Tissue <- factor(df_pred_probes$Tissue,
                                levels = c(metrics_tumor, metrics_normal))
ggplot(df_pred_probes, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  facet_wrap(~ Tissue) +
  labs(title = "Ridge Regression on Probe-level Methylation: Tumor vs Normal",
    x = "Observed Age",
    y = "Predicted Age")

# Apply a published epigenetic clock

clock <- meth_probes$horvath$coefficients
head(clock)
clock_intercept <- meth_probes$horvath$intercept

# Linear combination of the clock probes for every sample
clock_score <- clock_intercept + probe.data[, clock$probe] %*% clock$coefficient

# Horvath's inverse age transformation
horvath_age <- ifelse(clock_score <= 0, 21 * exp(clock_score) - 1, 21 * clock_score + 20)

df_clock <- data.frame(Observed = probe_meta$AGE,
                       DNAm_Age = as.numeric(horvath_age),
                       Tissue = probe_meta$tissue)
# Correlation with age in each tissue
by(df_clock, df_clock$Tissue, function(d) round(cor(d$Observed, d$DNAm_Age), 3))
ggplot(df_clock, aes(x = Observed, y = DNAm_Age)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  facet_wrap(~ Tissue) +
  labs(title = "Horvath Epigenetic Clock: Tumor vs Normal Tissue",
    x = "Observed Age",
    y = "DNA Methylation Age")

# -----------------------------------------------------------------------------
# Predict survival probability
# -----------------------------------------------------------------------------

brca_cli <- readRDS(url("https://wd.cri.uic.edu/machine_learning/brca_clinical.rds"))
brca.df.var_filt <- readRDS(url(
  "https://wd.cri.uic.edu/machine_learning/brca.df.var_filt.rds"))

# Get samples with complete cases of DFS_MONTHS and DFS_STATUS from the clinical data
keep <- complete.cases(brca_cli$DFS_MONTHS, 
                       brca_cli$DFS_STATUS) & brca_cli$DFS_MONTHS >= 0
DFS_MONTHS_clean <- brca_cli$DFS_MONTHS[keep]
summary(DFS_MONTHS_clean)
DFS_STATUS_clean <- brca_cli$DFS_STATUS[keep]
table(DFS_STATUS_clean)
brca.clean.df <- brca.df.var_filt[keep, ]
library(survival)
library(survcomp)
library(survminer)

# Adjust 0 times
DFS_MONTHS_clean[DFS_MONTHS_clean == 0 & 
                   DFS_STATUS_clean == "1:Recurred/Progressed"] <- 0.001

# Drop 0-month censored patients
keep2 <- !(DFS_MONTHS_clean == 0 & DFS_STATUS_clean == "0:DiseaseFree")

DFS_MONTHS_clean <- DFS_MONTHS_clean[keep2]
DFS_STATUS_clean <- DFS_STATUS_clean[keep2]
brca.clean.df <- brca.clean.df[keep2, ]

# Build survival object
status <- ifelse(DFS_STATUS_clean == "1:Recurred/Progressed", 1, 0)
head(status)
y <- Surv(time = DFS_MONTHS_clean, event = status)
head(y)
x <- as.matrix(brca.clean.df)
library(glmnet)

set.seed(123)
# subset data for a faster processing time
subsample <- sample(1:length(y),400)
x <- x[subsample, ]
y <- y[subsample]
DFS_MONTHS_clean <- DFS_MONTHS_clean[subsample]
status <- status[subsample]
# fit the model
coxfit <- cv.glmnet(
  x, 
  y, 
  family = "cox",
  alpha = 1,              # LASSO
  nfolds = 5              # 5-fold CV
)
coxfit
# Best lambda
best_lambda <- coxfit$lambda.min

# Refit model at best lambda
coxfinal <- glmnet(x, 
                   y, 
                   family = "cox", 
                   alpha = 1, 
                   lambda = best_lambda)

# Get nonzero coefficients
cox_selected_genes <- rownames(coef(coxfinal))[coef(coxfinal)[,1] != 0]
length(cox_selected_genes)
cox_selected_genes

# Predicted risk score
risk_scores <- predict(coxfinal, newx = x, type = "link")

cindex <- concordance.index(
  x = risk_scores, 
  surv.time = DFS_MONTHS_clean,
  surv.event = status
)

cindex$c.index

# Visualize Kaplan-Meier curves
risk_group <- ifelse(risk_scores > median(risk_scores), "High", "Low")
df.km <- data.frame(
  DFS_MONTHS = DFS_MONTHS_clean,
  status = status,
  risk_group = risk_group
)
fit_km <- survfit(Surv(DFS_MONTHS, status) ~ risk_group, data = df.km)
ggsurvplot(fit_km, data = df.km, pval = TRUE, risk.table = TRUE, risk.table.height = 0.25)

# =============================================================================
# TAKE-HOME EXERCISES
# =============================================================================

# -----------------------------------------------------------------------------
# Unbalanced model training (take home)
# -----------------------------------------------------------------------------

# --- Create unbalanced dataset ----------------------------------------------

library(tidyverse)

# Use complete set of ductal and lobular samples
ductal.samples.all <- brca_cli.ductal$sampleId
ductal.exp.all <- all_brca_RNA.log2.filtered[,ductal.samples.all]
lobular.samples.all <- brca_cli.lobular$sampleId
lobular.exp.all <- all_brca_RNA.log2.filtered[,lobular.samples.all]

# Combine ductal and lobular expression datasets and convert to matrix
expr_combined.all <- cbind(ductal.exp.all, lobular.exp.all)
n_ductal  <- ncol(ductal.exp.all)
n_lobular <- ncol(lobular.exp.all)
labels.all <- c(rep("Ductal", n_ductal), rep("Lobular", n_lobular))
table(labels.all)

# Build data.frame for modeling: samples x genes
df.all <- as.data.frame(t(expr_combined.all))  # rows = samples, cols = genes
df.all <- tibble::rownames_to_column(df.all, var = "sample_id")
df.all$class <- factor(labels.all, levels = c("Ductal","Lobular"))
colnames(df.all) <- make.names(colnames(df.all))

# --- Make a stratified training and testing set -----------------------------

library(caret)
library(dplyr)

set.seed(123)  # for reproducibility

# Make stratified split (70% train, 30% test)
# Note: When you pass a factor or outcome vector, createDataPartition() 
# automatically does stratified sampling

train_index <- createDataPartition(df.all$class, p = 0.7, list = FALSE)

train_data <- df.all[train_index, ]
test_data  <- df.all[-train_index, ]

# For convenience, create x and y
x_train.all <- train_data %>% select(-sample_id, -class)
y_train.all <- train_data$class
x_test.all  <- test_data %>% select(-sample_id, -class)
y_test.all  <- test_data$class

# Check class distribution
prop.table(table(df.all$class))       # overall
prop.table(table(train_data$class))   # training
prop.table(table(test_data$class))    # testing

# --- Feature elimination and RFE with logistic regression -------------------

library(caret)
library(ggplot2)

# T-test to rank remaining features
ttest_ranks.all <- apply(x_train.all, 2, function(g) {
  t.test(g ~ y_train.all)$statistic
})
ttest_ranks.all <- sort(abs(ttest_ranks.all), decreasing = TRUE)

# Pick the top 100 genes
top_ttest_genes.all <- names(ttest_ranks.all)[1:100]
df_ttest_filtered.all <- x_train.all[, top_ttest_genes.all]

# Number of features to be retained in the updated model
sizes <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

# Perform RFE
set.seed(123)
rfe_ctrl <- rfeControl(functions = lrFuncs, 
                   method = "cv", # cross validation
                   number = 10) # number of folds
rfe_results_all <- rfe(
  x = df_ttest_filtered.all,
  y = y_train.all,
  sizes = sizes,   # candidate numbers of features
  rfeControl = rfe_ctrl
)

best_rfe_feats_all<- predictors(rfe_results_all)
best_rfe_feats_all

# --- Calculate F1 statistics on test set ------------------------------------

library(caret)
library(MLmetrics)

# Subset training and testing data to the selected features
Xtr.all <- x_train.all[, best_rfe_feats_all]
Xte.all  <- x_test.all[, best_rfe_feats_all]

# Make training control function
lr_ctrl_F1 <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE,
  savePredictions = "final", 
  summaryFunction = function(data, lev = NULL, model = NULL) {
    f1 <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[2])
    c(F1 = f1)
    }
  )

# Train logistic regression with 10-fold CV
set.seed(123)
lr_model2 <- train(
  x = Xtr.all,
  y = y_train.all,
  method = "glm",
  family = "binomial",
  metric = "F1",
  trControl = lr_ctrl_F1
)

# Predict on the test set
y_pred <- predict(lr_model2, newdata = Xte.all)

# Test-set F1 score
f1_test <- F1_Score(y_pred, y_test.all, positive = "Lobular")
f1_test

# Cross-validation predictions
lr2_train_preds <- lr_model2$pred$pred
lr2_train_obs <- lr_model2$pred$obs

lr_model2_cm_cv <- confusionMatrix(lr2_train_preds, lr2_train_obs, positive = "Lobular")
lr_model2_cm_cv
lr_model2_cm_test <- confusionMatrix(y_pred, y_test.all, positive = "Lobular")
lr_model2_cm_test

# -----------------------------------------------------------------------------
# Random Forest regression (take home)
# -----------------------------------------------------------------------------

library(glmnet)
library(caret)
library(ggplot2)
# read in data
all_brca_meth <- readRDS(url(
  "https://wd.cri.uic.edu/machine_learning/brca_methylation_hm450.rds"))
brca_cli <- readRDS(url("https://wd.cri.uic.edu/machine_learning/brca_clinical.rds"))
rownames(all_brca_meth) <- all_brca_meth$hugoGeneSymbol
all_brca_meth <- as.matrix(all_brca_meth[,-1:-3])
# remove genes with missing values and convert to M-values
all_brca_meth <- all_brca_meth[rowSums(is.na(all_brca_meth)) == 0,]
all_brca_meth.M <- log2(all_brca_meth / (1 - all_brca_meth))
meth.data <- t(all_brca_meth.M)
colnames(meth.data) <- make.names(colnames(meth.data))
# match age by sample ID and remove samples with missing age
meth_age <- brca_cli$AGE[match(rownames(meth.data), brca_cli$sampleId)]
keep <- !is.na(meth_age)
meth_age <- meth_age[keep]
meth_age.data <- meth.data[keep, ]
# top variable genes
meth.gene_vars <- apply(meth_age.data, 2, var)
meth.top_var_genes <- names(sort(meth.gene_vars, decreasing = TRUE))[1:5000]
meth.df.var_filt <- meth_age.data[, meth.top_var_genes]
# make training and test sets
set.seed(123)
age_idx  <- createDataPartition(meth_age, p = 0.7, list = F)
age_x_train  <- meth.df.var_filt[age_idx, ]
age_y_train <- meth_age[age_idx]
age_x_test   <- meth.df.var_filt[-age_idx, ]
age_y_test <- meth_age[-age_idx]
# model training
rf_reg_ctrl <- trainControl(
  method = "cv",
  number = 5,
  savePredictions = "final"
)
set.seed(123)
rf_reg_model <- train(
  age_y_train ~ .,
  data = data.frame(age_y_train, age_x_train),
  method = "rf",
  tuneGrid = expand.grid(mtry=c(100, 500)), # with 5000 genes, try larger mtry values
  trControl = rf_reg_ctrl
)
rf_reg_model
# model evaluation
preds_rf <- predict(rf_reg_model, newdata = age_x_test)
# Performance (R², RMSE)
r2_rf <- cor(preds_rf, age_y_test)^2
r2_rf
rmse_rf <- sqrt(mean((preds_rf - age_y_test)^2))
rmse_rf
# Make a plot of regression
df_pred_rf <- data.frame(
  Observed = age_y_test,
  Predicted = as.numeric(preds_rf)
)
metrics_rf <- paste0("R2 = ", round(r2_rf, 3),
                       " RMSE = ", round(rmse_rf, 2))
ggplot(df_pred_rf, aes(x = Observed, y = Predicted)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  labs(title = "Random Forest Regression: Predicted vs Observed Age",
    subtitle = metrics_rf,
    x = "Observed Age",
    y = "Predicted Age")

