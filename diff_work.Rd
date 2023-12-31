% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_work.R
\name{diff_work}
\alias{diff_work}
\title{Performs differential analysis of the omics file imported}
\usage{
diff_work(file, metadata, patient_factor)
}
\arguments{
\item{file}{the csv file that contains the omics data. Rownames should be in format, condition_1. First column should include the binary sample conditions.}

\item{metadata}{another csv file that contains the patient metadata, including patient factors and disease measures that are to be investigated}

\item{patient_factor}{the patient variable against which the omic features will be correlated.#'}
}
\value{
Returns a list of objects. \link{1} returns a toptable from limma package with statistics. \link{2} returns a volcano plot showing log FC vs -log adjusted p-value
}
\description{
Performs differential analysis of the omics file imported
}
\examples{
## Read data from a CSV file
file <- read.csv(system.file("extdata", "example_data.csv", package = "markerHuntR"))
diff_plot <- diff_work(file, metadata, patient_factor)

## Alternatively, provide a data frame directly
file <- data.frame(...)
diff_plot <- diff_work(file, metadata, patient_factor)

}
