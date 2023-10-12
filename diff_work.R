#' Performs differential analysis of the omics file imported
#' @param file the csv file that contains the omics data. Rownames should be in format, condition_1. First column should include the binary sample conditions.
#' @param metadata another csv file that contains the patient metadata, including patient factors and disease measures that are to be investigated
#' @param patient_factor the patient variable against which the omic features will be correlated.#'
#' @return Returns a list of objects. [1] returns a toptable from limma package with statistics. [2] returns a volcano plot showing log FC vs -log adjusted p-value
#' @examples
#' ## Read data from a CSV file
#' file <- read.csv(system.file("extdata", "example_data.csv", package = "markerHuntR"))
#' diff_plot <- diff_work(file, metadata, patient_factor)
#'
#' ## Alternatively, provide a data frame directly
#' file <- data.frame(...)
#' diff_plot <- diff_work(file, metadata, patient_factor)
#'
#' @export
#' @import tidyverse
#' @import ggplot2
#' @import ggpubr
#' @import limma
#' @import purrr
#' @import data.table
#' @import ggrepel

diff_work <- function(file, metadata, patient_factor) {
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
  metadata <- metadata

  # Define the patient_factor based on user input
  # Replace 'Disease_Activity' with the actual user-defined factor variable
  patient_factor <- "Disease_Activity"

  number_cols_meta <- length(metadata)

  peak_df[1:5,1:5]
  metadata[1:5,]

  ## Wrangle
  combined <- left_join(metadata, peak_df, by=c('SampleID'))
  combined$patient_factor <- factor(combined[, patient_factor])
  combined[1:5,1:5]

  combined <- combined %>%
    dplyr::select(1, !!rlang::sym(patient_factor), names(peak_df)[-c(1:2)])

  names(combined)[2] <- 'Patient_Factor'
  rownames(combined) <- paste0(combined$Patient_Factor, "_", 1:nrow(combined))
  combined$Patient_Factor <- as.factor(combined$Patient_Factor)

  condition_1 <- levels(combined$Patient_Factor)[1]
  condition_2 <- levels(combined$Patient_Factor)[2]
  condition_1;condition_2

  rownames(combined) <- paste0(combined$Patient_Factor, "_", 1:nrow(combined))

  # reformat for limma
  pre_limma <- combined[,-c(1,2)]

  pre_limma <- as.data.frame(t(pre_limma))
  colnames(pre_limma) <- gsub("(.*)_.*", "\\1", colnames(pre_limma))

  Group <- factor(colnames(pre_limma), levels = c(`condition_1`, `condition_2`))
  design <- model.matrix(~Group)
  colnames(design) <- c("condition_1", "condition_1vscondition_2")
  eset <- as.matrix(pre_limma)
  fit <- limma::lmFit(eset, design)
  fit <- eBayes(fit)
  toptable <- topTable(fit, coef = "condition_1vscondition_2", adjust = "BH", number = nrow(pre_limma))
  toptable <- as.data.frame(toptable)
  toptable$Feature <- rownames(toptable)
  toptable <- toptable[, c(ncol(toptable), 1:(ncol(toptable) - 1))]
  toptable$Sig <- 0
  toptable$Sig <- ifelse(toptable$adj.P.Val < 0.05, "< 0.05", "> 0.05")
  toptable$Sig_Names <- 0
  toptable$Sig_Names <- ifelse(toptable$Sig == "< 0.05", toptable$Feature, "")

  toptable$Sig <- 0
  toptable$Sig <- ifelse(toptable$adj.P.Val < 0.05, "< 0.05", "> 0.05")
  trunc_limma <- head(toptable, 20)

  volcano_plot <- toptable %>%
    filter(toptable$adj.P.Val != "NA") %>%
    mutate(Sig = ifelse(toptable$logFC > 0 & toptable$adj.P.Val < 0.05, yes = "Up",
      ifelse(toptable$logFC < 0 & toptable$adj.P.Val < 0.05,
        yes = "Down",
        no = "Below Cut-Off"
      )
    )) %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(
      size = 3, alpha = 0.7,
      aes(
        colour = Sig,
        group = Sig
      )
    ) +
    theme_minimal() +
    labs(
      x = "LogFC",
      y = "-Log p-value",
      colour = ""
    ) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text(aes(x = 2.5, y = -log10(0.05)),
      label = "Adjusted p-value = 0.05",
      vjust = 2.1, colour = "dark grey"
    ) +
    theme(
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 14)
    ) +
    scale_colour_manual(values = c("grey", "indianred1", "steelblue1"))
  volcano_plot
 result_list <- list(trunc_limma = trunc_limma, volcano_plot = volcano_plot)
 return(result_list)
}
