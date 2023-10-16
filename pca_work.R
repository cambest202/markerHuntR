#' Performs principal components analysis from the omics datafile imported
#' @param file the csv file that contains the omics data. Rownames should be in format, condition_1. First column should include the binary sample conditions.
#' @param metadata another csv file that contains the patient metadata, including patient factors and disease measures that are to be investigated
#' @param patient_factor the patient variable against which the omic features will be correlated.#'
#' @return Returns a PCA plot with the samples labelled based on the binary conditions from the dataset.
#' @examples
#' ## Read data from a CSV file
#' file <- read.csv(system.file("extdata", "example_data.csv", package = "markerHuntR"))
#' pca_plot <- pca_work(file, metadata, patient_factor)
#'
#' ## Alternatively, provide a data frame directly
#' file <- data.frame(...)
#' pca_plot <- pca_work(file, metadata, patient_factor)
#'
#' @export
#' @import tidyverse
#' @import ggplot2
#' @import ggpubr
#' @import parallel
#' @import doParallel
#' @import purrr
#' @import patchwork
#' @import data.table
#' @import ggrepel

pca_work <- function(file, metadata, patient_factor) {
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
  # file <- read.csv('/Users/cameronbest/MRC_DTP_in_Precision_Medicine/Project/RA/TaSER_DAS28/20230926_taser_peaks.csv',
  #                 check.names = F, row.names = 1)
  peak_df <- file

  # metadata_df <- read.csv('/Users/cameronbest/MRC_DTP_in_Precision_Medicine/Project/RA/TaSER_DAS28/20230926_taser_metadata.csv',
  #                      check.names = F, row.names = 1)
  metadata <- metadata_df

  peak_df[1:5,1:5]
  metadata[1:5,]

  ## Wrangle
  combined <- dplyr::left_join(metadata, peak_df, by=c('SampleID'))
  combined$patient_factor <- factor(combined[, patient_factor])
  combined$Condition <- as.factor(combined$Condition)
  combined[1:5,1:5]

  combined <- combined %>%
    dplyr::select(1, !!rlang::sym('Condition'), names(peak_df)[-c(1:2)])

  names(combined)[2] <- 'Patient_Factor'
  rownames(combined) <- paste0(combined$Patient_Factor, "_", 1:nrow(combined))
  combined$Patient_Factor <- as.factor(combined$Patient_Factor)

  condition_1 <- levels(combined$Patient_Factor)[1]
  condition_2 <- levels(combined$Patient_Factor)[2]
  condition_1;condition_2

  ### PCA
  pca_mets <- peak_df[, -c(1,2)]
  scaled_intensities <- scale(pca_mets)
  scaled_intensities[do.call(cbind, lapply(scaled_intensities, is.nan))] <- 0
  scaled_intensities <- as.data.frame(scaled_intensities)
  pca_data <- prcomp(scaled_intensities)
  pca_coord <- data.frame(pca_data$x)
  var_explained <- pca_data$sdev^2 / sum(pca_data$sdev^2)

  pca_coord$patient_factor <- combined$Patient_Factor
  pca_coord <- pca_coord %>%
    relocate(patient_factor)

  pca_plot <- pca_coord %>%
    ggplot() +
    geom_point(
      size = 3, alpha = 0.7,
      aes(x = PC1, y = PC2, colour = patient_factor, fill = patient_factor)
    ) +
    labs(
      x = paste0("PC1: ", round(var_explained[1] * 100, 1), "%"),
      y = paste0("PC2: ", round(var_explained[2] * 100, 1), "%"),
      colour = "Condition",
      fill = "Condition"
    ) +
    geom_hline(
      yintercept = 0,
      colour = "navy",
      linetype = "dashed"
    ) +
    geom_vline(
      xintercept = 0,
      colour = "navy",
      linetype = "dashed"
    ) +
    theme_minimal() +
    theme(
      legend.key.size = unit(01, "cm"), # change legend key size
      legend.key.height = unit(1, "cm"), # change legend key height
      legend.key.width = unit(1, "cm"), # change legend key width
      legend.title = element_blank(), # change legend title font size
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 16)
    )
  pca_plot
  return(pca_plot)
}
