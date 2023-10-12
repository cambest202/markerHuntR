#' Several helpful functions used throughout analyses
#' @name utils
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
#' @import flextable
#' @import DALEX
#' @import rlang

`%notin%` <- Negate(`%in%`)
`%notlike%` <-  Negate(data.table::`%like%`)

theme_set(theme_minimal())
back_2_front <- function(df){
  df <- df[,c(ncol(df),1:(ncol(df)-1))]
  return(df)
}
flextable_only <- function(table){
  table_flex <- flextable::flextable(table)%>%
    flextable::fontsize(size = 8)%>%
    flextable::autofit()%>%
    flextable:: height( height = .5)%>%
    flextable::theme_vanilla()
}
