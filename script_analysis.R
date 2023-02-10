##### setup #####

library(targets)
library(tarchetypes)

source("./R/functions_analysis.R")

tar_option_set(packages = c("qs", 
                            "stats", 
                            "igraph", 
                            "Matrix", 
                            "S4Vectors", 
                            "tidyverse", 
                            "rmarkdown", 
                            "tidymodels", 
                            "SummarizedExperiment", 
                            "SingleCellExperiment"), 
               error = "continue", 
               memory = "transient", 
               garbage_collection = TRUE, 
               format = "qs")

# targets
list(
  ##### MODEL METRICS #####
  ### FILES ###
  tar_target(monocle3_metrics_single_file, "Data/monocle3_Metrics_Single_Subject.Rds", format = "file"),
  tar_target(monocle3_metrics_multi_file, "Data/monocle3_Metrics_Multi_Subject.Rds", format = "file"),
  tar_target(tradeSeq_metrics_single_file, "Data/tradeSeq_Metrics_Single_Subject.Rds", format = "file"),
  tar_target(tradeSeq_metrics_multi_file, "Data/tradeSeq_Metrics_Multi_Subject.Rds", format = "file"),
  tar_target(scLANE_GLM_metrics_single_file, "Data/scLANE_GLM_Metrics_Single_Subject.Rds", format = "file"),
  tar_target(scLANE_GLM_metrics_multi_file, "Data/scLANE_GLM_Metrics_Multi_Subject.Rds", format = "file"),
  tar_target(scLANE_GEE_metrics_multi_file, "Data/scLANE_GEE_Metrics_Multi_Subject.Rds", format = "file"),
  tar_target(Lamian_metrics_multi_file, "Data/Lamian_Metrics_Multi_Subject.Rds", format = "file"),
  
  ### OBJECTS ###
  tar_target(monocle3_metrics_single, readRDS(monocle3_metrics_single_file)),
  tar_target(monocle3_metrics_multi, readRDS(monocle3_metrics_multi_file)),
  tar_target(tradeSeq_metrics_single, readRDS(tradeSeq_metrics_single_file)),
  tar_target(tradeSeq_metrics_multi, readRDS(tradeSeq_metrics_multi_file)),
  tar_target(scLANE_GLM_metrics_single, readRDS(scLANE_GLM_metrics_single_file)),
  tar_target(scLANE_GLM_metrics_multi, readRDS(scLANE_GLM_metrics_multi_file)),
  tar_target(scLANE_GEE_metrics_multi, readRDS(scLANE_GEE_metrics_multi_file)),
  tar_target(Lamian_metrics_multi, readRDS(Lamian_metrics_multi_file)),
  
  ##### DATA CLEANING #####
  tar_target(metric_list, list(monocle3_metrics_single,
                               monocle3_metrics_multi, 
                               tradeSeq_metrics_single, 
                               tradeSeq_metrics_multi, 
                               scLANE_GLM_metrics_single, 
                               scLANE_GLM_metrics_multi, 
                               scLANE_GEE_metrics_multi, 
                               Lamian_metrics_multi)), 
  tar_target(metric_table, create_metric_table(metric.list = metric_list)), 
  
  ##### ANALYSIS #####
  tar_target(Figure_1, tar_render("Reports/Figure_1.Rmd"))
)