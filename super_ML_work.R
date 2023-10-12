#' Complete workflow for the supervised machine learning approach for generating a panel of analytes that discriminate two sample conditions.
#' super_ML_work provides the complete analysis of the omics dataset imported as a csv file
#' @param file the csv file that contains the omics data. Rownames should be in format, condition_1. First column should include the binary sample conditions.
#' @param metadata_df another csv file that contains the patient metadata, including patient factors and disease measures that are to be investigated
#' @param patient_factor the patient variable against which the omic features will be correlated.#'
#' @param train_test_split the modelling involves a holdout cross validation approach, involving the splitting of the data into training and testing splits.
#'     This designates the proportion used to train the model.
#' @param number_features User can provide a number of features to include in the model, generally based on the large_fs from the feature selection.
#' @param algorithm User can designate the supervised machine learning algorithm used in the final model. Generally selected from the best performing algorithm shown in the ROC curves and table.
#' @return A list containing the following components:
#'    \describe{
#'      \item{large_fs}{A bar plot showing the top 20 features from the feature selection process. Used to select the number of features.}
#'      \item{ft_actual}{A bar plot showing the selected features for the final model. Model was trained on 80% of the dataset, and evaluated in 20%.}
#'      \item{trainn_2}{The reduced dataframe containing the condition column and selected features}
#'      \item{met_table}{A table showing the performance metrics from seven commonly used supervised machine learning algorithms for the binary classification}
#'      \item{met_roc}{A plot of the ROC curves for each supervised machine learning algorithms. Complements the table showing the performance metrics}
#'      \item{met_modell}{The final tuned model with 10-fold cross validation repeated 10 times.}
#'      \item{final_roc_met}{The ROC curve for the final model evaluated within the test subset.}
#'      \item{partials}{Model-agnostic feature interpretation plots showing the inner workings of the model. Partial dependence plot and accumulated local effects plot shown for each feature.}
#'      \item{actuals}{Boxplots of the selected features' abundances/expression across the sample conditions. T-tests/Wilcoxon tests performed depending on the distribution of the features.}
#'    }
#' @examples
#' ## Read data from a CSV file
#' file <- read.csv(system.file("extdata", "example_data.csv", package = "markerHuntR"))
#' super_ML_work <- super_ML_work(file, metadata_df, train_test_split, number_features, algorithm)
#'
#' ## Alternatively, provide a data frame directly
#' file <- data.frame(...)
#' super_ML_work <- super_ML_work(file, metadata_df, train_test_split, number_features, algorithm)
#'
#' @export
#' @import tidyverse
#' @import ggplot2
#' @import ggpubr
#' @import caret
#' @import ranger
#' @import neuralnet
#' @import mboost
#' @import limma
#' @import parallel
#' @import doParallel
#' @import e1071
#' @import purrr
#' @import factoextra
#' @import xgboost
#' @import naivebayes
#' @import kknn
#' @import naivebayes
#' @import patchwork
#' @import discrim
#' @import klaR
#' @import ROCR
#' @import pROC
#' @import grid
#' @import gridExtra
#' @import MLeval
#' @import randomForest
#' @import data.table
#' @import ggrepel
#' @import DALEXtra
#' @import mice
#' @import stacks
#' @import VIM

super_ML_work <- function(file = file,
                          metadata_df= metadata_df,
                          patient_factor=patient_factor,
                          train_test_split=train_test_split,
                          number_features = number_features,
                          algorithm = algorithm) {
  if (is.character(file)) {
    ## Read data from a CSV file if a file path is provided
    peak_df <- read.csv(file, row.names = 1, check.names = FALSE)
  } else if (is.data.frame(file)) {
    ## Use the provided data frame as is
    peak_df <- file
  } else {
    stop("Invalid 'file' argument. Provide either a file path or a data frame.")
  }
  ## Check if patient_factor is a valid column in metadata_df
  if (!patient_factor %in% colnames(metadata_df)) {
    stop("The 'patient_factor' variable is not found in the metadata data frame.")
  }

  doParallel::registerDoParallel(cores=4)

  ### Data ----
  # insert data
  # file <- read.csv('/Users/cameronbest/postdoc/stroke/MD_stroke_df/20230927_stroke_peaks.csv',
  #                  check.names = F, row.names = 1)
  peak_df <- file

  # metadata_df <- read.csv('/Users/cameronbest/postdoc/stroke/MD_stroke_df/20230927_stroke_metadata.csv',
  #                      check.names = F, row.names = 1)
  metadata <- metadata_df
  metadata$Condition <- as.factor(metadata$Condition)

  # Define the patient_factor based on user input
  # Replace 'Disease_Activity' with the actual user-defined factor variable
  patient_factor <- "Disease_Activity"

  peak_df[1:5,1:5]
  metadata[1:5,]

  ## Wrangle
  combined <- left_join(metadata, peak_df, by=c('SampleID'))
  combined$patient_factor <- factor(combined[, patient_factor])
  combined[1:5,1:5]

  combined <- combined %>%
    dplyr::select(1, !!rlang::sym(patient_factor), names(peak_df)[-c(1:2)])

  names(combined)[2] <- 'Patient_Factor'
  combined$Patient_Factor <- as.factor(combined$Patient_Factor)

  condition_1 <- levels(combined$Patient_Factor)[1]
  condition_2 <- levels(combined$Patient_Factor)[2]
  condition_1;condition_2

  combined <- combined[,-1]

  ### Supervised Machine Learning
  set.seed(42)
  index <- createDataPartition(combined$Patient_Factor, p = train_test_split, list = FALSE)
  train_data <- combined[index, ]
  test_data <- combined[-index, ]

  # function for complete model development

  shared_metabolites <- function(train_data) {
    for (i in 1:10) {
      profile_function <- function(Profile_number) {
        set.seed(42)
        options(warn = -1)
        subsets <- c(2:10) # 10
        set.seed(42)
        ctrl <- rfeControl(
          functions = rfFuncs,
          method = "repeatedcv",
          number = 10,
          repeats = 5,
          verbose = FALSE
        )

        profile <- rfe(
          x = train_data[, -1], y = train_data$Patient_Factor,
          sizes = subsets,
          rfeControl = ctrl
        )
        profile$optVariables
        profile <- profile$variables
        profile <- profile %>%
          arrange(-Overall) %>%
          distinct(var, .keep_all = TRUE) %>%
          filter(Overall > 1)

        profile_df <- profile %>%
          group_by(var) %>%
          nest() %>%
          mutate(Profile = 1)
      }
      profile <- profile_function(i)
      assign(paste0("profile_test", i), profile)
    }
    df <- do.call(rbind, mget(ls(pattern = "profile_test")))

    df_dupl <- df %>%
      group_by(var) %>%
      filter(n() >= 5) %>%
      unnest()
    df_dupl <- df_dupl[, c(1, 4)]
    df_dupl <- df_dupl %>%
      aggregate(. ~ var, mean, na.rm = TRUE) %>%
      mutate(Round_overall = round(Overall, 1)) %>%
      arrange(-Overall)
  }

  ft_sel_profile <- shared_metabolites(train_data)

  shared_view <- ft_sel_profile %>%
    head(20)

  large_fs <- shared_view %>%
    ggplot(aes(x = Overall, y = reorder(var, Overall))) +
    geom_col(aes(fill = Overall)) +
    geom_text(aes(label = round(Overall, 2)),
      vjust = 0.8,
      hjust = 1.2,
      color = "white",
      size = 6
    ) +
    # scale_color_colorblind()+
    scale_fill_continuous(low = "light green", high = "navy") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.text.y = element_text(size = 16),
      title = element_text(size = 22)
    ) +
    labs(
      x = "Relative Importance",
      y = "Feature"
    )
  large_fs

  shared_actual <- ft_sel_profile %>%
    head(number_features)

  ft_actual <- shared_actual %>%
    ggplot(aes(x = Overall, y = reorder(var, Overall))) +
    geom_col(aes(fill = Overall)) +
    geom_text(aes(label = round(Overall, 2)),
      vjust = 0.8,
      hjust = 1.2,
      color = "white",
      size = 6
    ) +
    scale_fill_continuous(low = "light green", high = "navy") +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.text.y = element_text(size = 16),
      title = element_text(size = 22)
    ) +
    labs(
      x = "Relative Importance",
      y = "Putative Metabolite"
    )
  ft_actual

  trainn_2 <- train_data %>%
    dplyr::select(Patient_Factor, shared_actual$var) %>%
    as.data.frame()

  ## Modelling
  # Algorithm Selection
  # Function for algorithm tuning and selection
  multi_mods_test <- function() {
    control <- trainControl(
      method = "repeatedcv",
      number = 10,
      repeats = 10,
      summaryFunction = twoClassSummary,
      savePredictions = TRUE,
      classProbs = TRUE,
      verboseIter = TRUE,
      search = "random"
    )
    # train the SVM model
    set.seed(42)
    modelSvm <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "svmRadial",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )

    # train the LogM model
    set.seed(42)
    modelglm <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "glm",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )
    # train the LogM model
    set.seed(42)
    modelglm_boost <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "glmboost",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )
    # train the RF model
    set.seed(42)
    model_rf <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "ranger",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )
    # train the XGBoost model
    set.seed(42)
    model_xgb <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "xgbTree",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )
    # train the KNN model
    set.seed(42)
    model_knn <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "kknn",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )
    # train the naive Bayes model
    set.seed(42)
    model_naivebayes <- caret::train(Patient_Factor ~ .,
      data = trainn_2,
      method = "naive_bayes",
      metric = "ROC",
      tuneLength = 5,
      trControl = control
    )


    # collect resamples
    results <- resamples(list(
      SVM = modelSvm, RF = model_rf, LRM = modelglm, GLMB = modelglm_boost,
      XGBoost = model_xgb, KNN = model_knn,
      Naive.Bayes = model_naivebayes
    ))

    plott <- bwplot(results)
    comp_roc <- evalm(
      list(
        modelSvm, model_rf, modelglm, modelglm_boost,
        model_xgb, model_knn,
        model_naivebayes
      ),
      gnames = c("KNN", "NB", "RF", "LRM", "GLMB", "SVM", "XGB")
    )

    ml_eval_output <- as.data.frame(comp_roc$stdres)
    ml_eval_output$Measure <- rownames(ml_eval_output)
    ml_eval_output <- back_2_front(ml_eval_output)
    ml_eval_output <- flextable_only(ml_eval_output)
    listt <- list(comp_roc, ml_eval_output)
    return(listt)
  }

  alg_sel <- multi_mods_test()
  met_table <- alg_sel[[2]]
  met_roc <- alg_sel[[1]]$roc
  met_prg <- alg_sel[[1]]$prg

  # Model generation

  set.seed(42)
  fit_control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 10,
    summaryFunction = twoClassSummary,
    savePredictions = TRUE,
    classProbs = TRUE,
    verboseIter = TRUE,
    search = "random"
  )

  set.seed(42)
  met_modell <- caret::train(Patient_Factor ~ .,
    data = trainn_2,
    method = 'ranger',
    metric = "ROC",
    tuneLength = 15,
    trControl = fit_control
  )

  # testing subset
  # function
  model_performance <- function(model, test_df) {
    predictions <- predict(model, test_df)
    confusionMatrix(predictions, test_df$Patient_Factor)
    con_matr <- confusionMatrix(predictions, test_df$Patient_Factor)
    con_stats <- con_matr$overall

    pr <- prediction(as.numeric(predictions), as.numeric(test_df$Patient_Factor))
    prf <- performance(pr, measure = "tpr", x.measure = "fpr")
    auc <- performance(pr, measure = "auc")
    auc_val <- auc@y.values[[1]]
    result.predicted.prob <- predict(model, test_df, type = "prob") # Prediction

    Condition <- as.factor(condition_1)
    result.roc <- roc(test_df$Patient_Factor, result.predicted.prob[[1]]) # Draw ROC curve.
    list_pred <- list(model, con_stats, result.roc, con_matr)
    return(list_pred)
  }

  mod_5_guys <- model_performance(met_modell, test_data)
  roc_met_rfe <- mod_5_guys[[3]]
  roc_met_rfe$auc

  roc_met_df <- data.frame(roc_met_rfe$sensitivities, roc_met_rfe$specificities)
  names(roc_met_df) <- c("Sensitivity", "Specificity")
  roc_met_df$One_Minus_Spec <- 1 - roc_met_df$Specificity
  roc_met_df$Model <- "Metabolite-Only"

  roc_met_final <- roc_met_df %>%
    ggplot(aes(x = One_Minus_Spec, y = Sensitivity)) +
    geom_path(colour = "red", size = 1) +
    geom_abline(linetype = "solid", colour = "grey") +
    coord_equal() +
    theme_minimal() +
    labs(
      x = "False positive rate",
      y = "True positive"
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      legend.position = "none",
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20)
    )

  AUC_ROC_met <- ggplot() +
    geom_line(aes(x = c(-1:4), y = 0), colour = "red", size = 1.5) +
    xlim(-3, 50) +
    theme_void() +
    annotate("text",
      size = 6,
      x = 30,
      y = 0,
      label = paste0("AUC-ROC = ", round(roc_met_rfe$auc, 2)), colour = "black"
    )

  final_roc_met <- (roc_met_final + AUC_ROC_met + plot_layout(widths = c(2.5, 1)))
  final_roc_met

  # Feature Interpretation
  train_data_3 <- trainn_2
  train_data_3$Patient_Factor <- ifelse(train_data_3$Patient_Factor == condition_1, 0, 1)
  xplainer_rf <- DALEX::explain(
    model = met_modell,
    data = train_data_3[, -1],
    y = train_data_3$Patient_Factor,
    type = "classification"
  )

  ## attempt accumulated local effects plots to explain interaction of features in the plot
  pd_rf <- model_profile(
    explainer = xplainer_rf,
    type = "partial",
    variables = names(train_data_3)[2:ncol(train_data_3)]
  )
  ld_rf <- model_profile(
    explainer = xplainer_rf,
    type = "conditional",
    variables = names(train_data_3)[2:ncol(train_data_3)]
  )
  ale_rf <- model_profile(
    explainer = xplainer_rf,
    type = "accumulated",
    variables = names(train_data_3)[2:ncol(train_data_3)]
  )

  shapleys_rf <- predict_parts(
    explainer = xplainer_rf,
    type = "shap",
    new_observation = train_data_3,
    B = 25
  ) # the B=25 is the number of random orderings of the explanatory variables

  # shapleys_rf <- predict_parts(explainer = xplainer_rf, type='shap',new_observation =  train_data_3, B=25) # the B=25 is the number of random orderings of the explanatory variables

  pd_rf$agr_profiles$`_label_` <- "partial dependence"
  ld_rf$agr_profiles$`_label_` <- "local dependence"
  ale_rf$agr_profiles$`_label_` <- "accumulated local"

  partials <- plot(
    pd_rf,
    ale_rf
  )

  partials <- partials +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )
  partials

  # Actual changes

  actuals <- combined %>%
    dplyr::select(Patient_Factor, shared_actual$var) %>%
    pivot_longer(
      cols = 2:(number_features+1),
      names_to = "Feature",
      values_to = "Level"
    ) %>%
    ggplot(aes(x = Patient_Factor, y = Level, fill = Patient_Factor)) +
    # geom_violin()+
    geom_boxplot() +
    geom_jitter(size = 0.5, alpha = 0.7) +
    facet_wrap(~Feature,
      scale = "free", ncol = 5,
      labeller = label_wrap_gen(multi_line = TRUE)
    ) +
    stat_compare_means(
      method = "wilcox.test",
      label = "p.format",
      vjust = 0,
      hjust = -0.5,
      comparisons = list(
        c("Control", "Stage1"),
        c("Control", "Stage2"),
        c("Stage1", "Stage2")
      )
    ) +
    theme(
      legend.position = "none",
      strip.text.x = element_text(
        face = "bold",
        size = 12
      ),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12)
    )
  return(list(large_fs, ft_actual, trainn_2, met_table, met_roc, met_modell, final_roc_met, partials, actuals))
}
