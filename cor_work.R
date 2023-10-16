#' Takes the metadata from the file and performs correlation analysis based on the disease measure the user chooses
#' @param file the csv file that contains the omics data. Rownames should be in format, condition_1. First column should include the binary sample conditions.
#' @param metadata another csv file that contains the patient metadata, including patient factors and disease measures that are to be investigated
#' @param patient_factor the patient variable against which the omic features will be correlated.
#' @return Returns [1] table of correlation statistics [2] scatter plots showing the features significantly correlated with the patient factor or disease measure.
#'     Regression line and appropriate statistical tests are also included in the plot
#' @examples
#' ## Read data from a CSV file
#' file <- read.csv(system.file("extdata", "example_data.csv", package = "markerHuntR"))
#' metadata <- read.csv(system.file("extdata", "example_patient_metadata.csv", package = "markerHuntR"))

#' cor_plot <- cor_work(file, metadata, patient_factor)
#'
#' ## Alternatively, provide a data frame directly
#' file <- data.frame(...)
#' cor_plot <- cor_work(file, metadata, patient_factor)
#'
#' @export
#' @import tidyverse
#' @import ggplot2
#' @import ggpubr
#' @import purrr
#' @import factoextra
#' @import patchwork
#' @import discrim
#' @import grid
#' @import gridExtra
#' @import data.table
#' @import ggrepel

cor_work <- function(file, metadata, patient_factor) {
  if (is.character(file)) {
    ## Read data from a CSV file if a file path is provided
    peak_df <- read.csv(file, row.names = 1, check.names = FALSE)
  } else if (is.data.frame(file)) {
    ## Use the provided data frame as is
    peak_df <- file
  } else {
    stop("Invalid 'file' argument. Provide either a file path or a data frame.")
  }

  if (is.character(metadata)) {
    ## Read data from a CSV file if a file path is provided
    metadata_df <- read.csv(metadata, row.names = 1, check.names = FALSE)
  } else if (is.data.frame(metadata)) {
    ## Use the provided data frame as is
    metadata_df <- metadata
  } else {
    stop("Invalid 'metadata' argument. Provide either a file path or a data frame.")
  }
  ## Check if patient_factor is a valid column in metadata_df
  if (!patient_factor %in% colnames(metadata_df)) {
    stop("The 'patient_factor' variable is not found in the metadata data frame.")
  }

  ### Data ----
  # insert data
  file <- read.csv('/Users/cameronbest/MRC_DTP_in_Precision_Medicine/Project/RA/TaSER_DAS28/20230926_taser_peaks.csv',
                   check.names = F, row.names = 1)
  peak_df <- file

  metadata_df <- read.csv('/Users/cameronbest/MRC_DTP_in_Precision_Medicine/Project/RA/TaSER_DAS28/20230926_taser_metadata.csv',
                          check.names = F, row.names = 1)
  metadata <- metadata_df

  # Define the patient_factor based on user input
  # Replace 'Disease_Activity' with the actual user-defined factor variable
  patient_factor <- "Disease_Activity"

  peak_df[1:5,1:5]
  metadata[1:5,]

  ## Wrangle
  combined <- dplyr::left_join(metadata, peak_df, by=c('SampleID'))
  combined[1:5,1:5]

  combined <- combined %>%
    dplyr::select(1, !!rlang::sym(patient_factor), names(peak_df)[-c(1:2)])

  names(combined)[2] <- 'Patient_Factor'

  # format for correlation analysis
  combined_long <- combined %>%
    pivot_longer(cols=3:ncol(combined),
                 names_to='Feature',
                 values_to='Level')

  # perform correlation analysis
  ints_nested <- combined_long %>%
    group_by(Feature) %>%
    nest()

  ints_lm <- ints_nested %>%
    mutate(
      model = purrr::map(data, ~ {
        lm(formula = Level ~ Patient_Factor, data = .x)
      })
    )

  model_coef_nested <- ints_lm %>%
    mutate(coef = purrr::map(model, ~ tidy(.x)))
  model_coef <- model_coef_nested %>%
    unnest(coef)

  model_perf_nested <- ints_lm %>%
    mutate(fit = purrr::map(model, ~ glance(.x)))
  model_perf <- model_perf_nested %>%
    unnest(fit)

  best_fit <- model_perf %>%
    top_n(n = 4, wt = r.squared)
  bestest_fit <- with(model_perf, model_perf[order(-r.squared),])
  best_augmented <- bestest_fit %>%
    mutate(augmented = purrr::map(model, ~ augment(.x))) %>%
    unnest(augmented)
  best_augmented$adj_p <- p.adjust(best_augmented$p.value, method = 'BH')

  cor_table <- best_augmented[, -c(2, 3)]

  cor_plot <- best_augmented %>%
    filter(adj_p < 0.5) %>%
    ggplot(aes(x = Patient_Factor, y=Level)) +
    geom_point(size=2, alpha=0.7) +
    stat_cor(vjust=1, hjust=0,
             size=5)+
    geom_smooth(method='lm',
                colour='red')+
    facet_wrap(~Feature, scales = "free_y")+
    theme(
      strip.text.x= element_text(#face = "bold",
        size=12),
      axis.title = element_text(size=16),
      axis.text = element_text(size=14))
  cor_plot

  return(list(cor_table = cor_table,
              cor_plot = cor_plot))
}
