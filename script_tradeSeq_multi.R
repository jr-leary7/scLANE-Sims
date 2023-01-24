##### setup #####

library(future)
library(targets)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 100000 * 1024^2)  # 100GB max size -- necessary for loading in all sim datasets

source("./R/functions_tradeSeq.R")

tar_option_set(packages = c("mgcv", 
                            "pryr", 
                            "tidyr", 
                            "stats", 
                            "broom", 
                            "scran", 
                            "gamlss", 
                            "igraph",
                            "scater", 
                            "scRNAseq", 
                            "parallel", 
                            "tradeSeq", 
                            "S4Vectors", 
                            "tidyverse", 
                            "doParallel", 
                            "BiocParallel", 
                            "SummarizedExperiment", 
                            "SingleCellExperiment"),
               error = "continue", 
               retrieval = "main", 
               storage = "main", 
               memory = "transient",
               garbage_collection = TRUE, 
               format = "qs")

##### UPSTREAM TARGETS #####
### BRAIN -- MULTI-SUBJECT 
# 100 cells 
brain_sim_DEG_10_CELLS_100_balanced <- tar_read(brain_sim_DEG_10_CELLS_100_balanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_100_balanced <- tar_read(brain_sim_DEG_20_CELLS_100_balanced, store = "store_simulation")
brain_sim_DEG_10_CELLS_100_unbalanced <- tar_read(brain_sim_DEG_10_CELLS_100_unbalanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_100_unbalanced <- tar_read(brain_sim_DEG_20_CELLS_100_unbalanced, store = "store_simulation")
# 500 cells 
brain_sim_DEG_10_CELLS_500_balanced <- tar_read(brain_sim_DEG_10_CELLS_500_balanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_500_balanced <- tar_read(brain_sim_DEG_20_CELLS_500_balanced, store = "store_simulation")
brain_sim_DEG_10_CELLS_500_unbalanced <- tar_read(brain_sim_DEG_10_CELLS_500_unbalanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_500_unbalanced <- tar_read(brain_sim_DEG_20_CELLS_500_unbalanced, store = "store_simulation")
# 1000 cells 
brain_sim_DEG_10_CELLS_1000_balanced <- tar_read(brain_sim_DEG_10_CELLS_1000_balanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_1000_balanced <- tar_read(brain_sim_DEG_20_CELLS_1000_balanced, store = "store_simulation")
brain_sim_DEG_10_CELLS_1000_unbalanced <- tar_read(brain_sim_DEG_10_CELLS_1000_unbalanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_1000_unbalanced <- tar_read(brain_sim_DEG_20_CELLS_1000_unbalanced, store = "store_simulation")
# 2500 cells 
brain_sim_DEG_10_CELLS_2500_balanced <- tar_read(brain_sim_DEG_10_CELLS_2500_balanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_2500_balanced <- tar_read(brain_sim_DEG_20_CELLS_2500_balanced, store = "store_simulation")
brain_sim_DEG_10_CELLS_2500_unbalanced <- tar_read(brain_sim_DEG_10_CELLS_2500_unbalanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_2500_unbalanced <- tar_read(brain_sim_DEG_20_CELLS_2500_unbalanced, store = "store_simulation")
# 5000 cells 
brain_sim_DEG_10_CELLS_5000_balanced <- tar_read(brain_sim_DEG_10_CELLS_5000_balanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_5000_balanced <- tar_read(brain_sim_DEG_20_CELLS_5000_balanced, store = "store_simulation")
brain_sim_DEG_10_CELLS_5000_unbalanced <- tar_read(brain_sim_DEG_10_CELLS_5000_unbalanced, store = "store_simulation")
brain_sim_DEG_20_CELLS_5000_unbalanced <- tar_read(brain_sim_DEG_20_CELLS_5000_unbalanced, store = "store_simulation")

### PANCREAS -- MULTI-SUBJECT 
# 100 cells 
panc_sim_DEG_10_CELLS_100_balanced <- tar_read(panc_sim_DEG_10_CELLS_100_balanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_100_balanced <- tar_read(panc_sim_DEG_20_CELLS_100_balanced, store = "store_simulation")
panc_sim_DEG_10_CELLS_100_unbalanced <- tar_read(panc_sim_DEG_10_CELLS_100_unbalanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_100_unbalanced <- tar_read(panc_sim_DEG_20_CELLS_100_unbalanced, store = "store_simulation")
# 500 cells 
panc_sim_DEG_10_CELLS_500_balanced <- tar_read(panc_sim_DEG_10_CELLS_500_balanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_500_balanced <- tar_read(panc_sim_DEG_20_CELLS_500_balanced, store = "store_simulation")
panc_sim_DEG_10_CELLS_500_unbalanced <- tar_read(panc_sim_DEG_10_CELLS_500_unbalanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_500_unbalanced <- tar_read(panc_sim_DEG_20_CELLS_500_unbalanced, store = "store_simulation")
# 1000 cells 
panc_sim_DEG_10_CELLS_1000_balanced <- tar_read(panc_sim_DEG_10_CELLS_1000_balanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_1000_balanced <- tar_read(panc_sim_DEG_20_CELLS_1000_balanced, store = "store_simulation")
panc_sim_DEG_10_CELLS_1000_unbalanced <- tar_read(panc_sim_DEG_10_CELLS_1000_unbalanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_1000_unbalanced <- tar_read(panc_sim_DEG_20_CELLS_1000_unbalanced, store = "store_simulation")
# 2500 cells 
panc_sim_DEG_10_CELLS_2500_balanced <- tar_read(panc_sim_DEG_10_CELLS_2500_balanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_2500_balanced <- tar_read(panc_sim_DEG_20_CELLS_2500_balanced, store = "store_simulation")
panc_sim_DEG_10_CELLS_2500_unbalanced <- tar_read(panc_sim_DEG_10_CELLS_2500_unbalanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_2500_unbalanced <- tar_read(panc_sim_DEG_20_CELLS_2500_unbalanced, store = "store_simulation")
# 5000 cells 
panc_sim_DEG_10_CELLS_5000_balanced <- tar_read(panc_sim_DEG_10_CELLS_5000_balanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_5000_balanced <- tar_read(panc_sim_DEG_20_CELLS_5000_balanced, store = "store_simulation")
panc_sim_DEG_10_CELLS_5000_unbalanced <- tar_read(panc_sim_DEG_10_CELLS_5000_unbalanced, store = "store_simulation")
panc_sim_DEG_20_CELLS_5000_unbalanced <- tar_read(panc_sim_DEG_20_CELLS_5000_unbalanced, store = "store_simulation")

### ENDOCRINOGENESIS -- MULTI-SUBJECT 
# 100 cells 
endo_sim_DEG_10_CELLS_100_balanced <- tar_read(endo_sim_DEG_10_CELLS_100_balanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_100_balanced <- tar_read(endo_sim_DEG_20_CELLS_100_balanced, store = "store_simulation")
endo_sim_DEG_10_CELLS_100_unbalanced <- tar_read(endo_sim_DEG_10_CELLS_100_unbalanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_100_unbalanced <- tar_read(endo_sim_DEG_20_CELLS_100_unbalanced, store = "store_simulation")
# 500 cells 
endo_sim_DEG_10_CELLS_500_balanced <- tar_read(endo_sim_DEG_10_CELLS_500_balanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_500_balanced <- tar_read(endo_sim_DEG_20_CELLS_500_balanced, store = "store_simulation")
endo_sim_DEG_10_CELLS_500_unbalanced <- tar_read(endo_sim_DEG_10_CELLS_500_unbalanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_500_unbalanced <- tar_read(endo_sim_DEG_20_CELLS_500_unbalanced, store = "store_simulation")
# 1000 cells 
endo_sim_DEG_10_CELLS_1000_balanced <- tar_read(endo_sim_DEG_10_CELLS_1000_balanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_1000_balanced <- tar_read(endo_sim_DEG_20_CELLS_1000_balanced, store = "store_simulation")
endo_sim_DEG_10_CELLS_1000_unbalanced <- tar_read(endo_sim_DEG_10_CELLS_1000_unbalanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_1000_unbalanced <- tar_read(endo_sim_DEG_20_CELLS_1000_unbalanced, store = "store_simulation")
# 2500 cells 
endo_sim_DEG_10_CELLS_2500_balanced <- tar_read(endo_sim_DEG_10_CELLS_2500_balanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_2500_balanced <- tar_read(endo_sim_DEG_20_CELLS_2500_balanced, store = "store_simulation")
endo_sim_DEG_10_CELLS_2500_unbalanced <- tar_read(endo_sim_DEG_10_CELLS_2500_unbalanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_2500_unbalanced <- tar_read(endo_sim_DEG_20_CELLS_2500_unbalanced, store = "store_simulation")
# 5000 cells 
endo_sim_DEG_10_CELLS_5000_balanced <- tar_read(endo_sim_DEG_10_CELLS_5000_balanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_5000_balanced <- tar_read(endo_sim_DEG_20_CELLS_5000_balanced, store = "store_simulation")
endo_sim_DEG_10_CELLS_5000_unbalanced <- tar_read(endo_sim_DEG_10_CELLS_5000_unbalanced, store = "store_simulation")
endo_sim_DEG_20_CELLS_5000_unbalanced <- tar_read(endo_sim_DEG_20_CELLS_5000_unbalanced, store = "store_simulation")

# tradeSeq model targets 
list(
  ##### BRAIN #####
  ### MULTI-SUBJECT 
  # 100 cells
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_100_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_100_balanced,
                                                                              n.genes.sample = 2000, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 100, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_100_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_100_balanced,
                                                                              n.genes.sample = 1200, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 100, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_100_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_100_unbalanced,
                                                                                n.genes.sample = 1800, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 100, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_100_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_100_unbalanced,
                                                                                n.genes.sample = 3100, 
                                                                                n.iter = 1,
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 100, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  # 500 cells
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_500_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_500_balanced,
                                                                              n.genes.sample = 3000, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 500, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_500_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_500_balanced,
                                                                              n.genes.sample = 2600, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 500, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_500_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_500_unbalanced,
                                                                                n.genes.sample = 2200, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 500, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_500_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_500_unbalanced,
                                                                                n.genes.sample = 1500, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 500, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  # 1000 cells
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_1000_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_1000_balanced,
                                                                               n.genes.sample = 1400, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 1000, Allocation = "balanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_1000_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_1000_balanced,
                                                                               n.genes.sample = 2000, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 1000, Allocation = "balanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_1000_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_1000_unbalanced,
                                                                                 n.genes.sample = 3100, 
                                                                                 param.list = list(Prop_Dyn_Genes = 0.1, Cells = 1000, Allocation = "unbalanced", Method = "GAM"),
                                                                                 n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_1000_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_1000_unbalanced,
                                                                                 n.genes.sample = 1900, 
                                                                                 param.list = list(Prop_Dyn_Genes = 0.2, Cells = 1000, Allocation = "unbalanced", Method = "GAM"),
                                                                                 n.cores = 4)), 
  # 2500 cells
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_2500_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_2500_balanced,
                                                                               n.genes.sample = 1900, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 2500, Allocation = "balanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_2500_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_2500_balanced,
                                                                               n.genes.sample = 1700, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 2500, Allocation = "balanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_2500_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_2500_unbalanced,
                                                                                 n.genes.sample = 2800, 
                                                                                 param.list = list(Prop_Dyn_Genes = 0.1, Cells = 2500, Allocation = "unbalanced", Method = "GAM"),
                                                                                 n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_2500_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_2500_unbalanced,
                                                                                 n.genes.sample = 1600,
                                                                                 param.list = list(Prop_Dyn_Genes = 0.2, Cells = 2500, Allocation = "unbalanced", Method = "GAM"),
                                                                                 n.cores = 4)), 
  # 5000 cells
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_5000_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_5000_balanced,
                                                                               n.genes.sample = 2200, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 5000, Allocation = "balanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_5000_balanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_5000_balanced,
                                                                               n.genes.sample = 1100, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 5000, Allocation = "balanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_10_CELLS_5000_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_10_CELLS_5000_unbalanced,
                                                                                 n.genes.sample = 1800, 
                                                                                 param.list = list(Prop_Dyn_Genes = 0.1, Cells = 5000, Allocation = "unbalanced", Method = "GAM"),
                                                                                 n.cores = 4)),
  tar_target(tradeSeq_res_brain_DEG_20_CELLS_5000_unbalanced, run_tradeSeq_multi(sim.data = brain_sim_DEG_20_CELLS_5000_unbalanced,
                                                                                 n.genes.sample = 2500, 
                                                                                 param.list = list(Prop_Dyn_Genes = 0.2, Cells = 5000, Allocation = "unbalanced", Method = "GAM"),
                                                                                 n.cores = 4)), 
  
  ##### PANCREAS #####
  ### MULTI-SUBJECT 
  # 100 cells
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_100_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_100_balanced,
                                                                             n.genes.sample = 2300, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.1, Cells = 100, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_100_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_100_balanced,
                                                                             n.genes.sample = 1200, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.2, Cells = 100, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_100_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_100_unbalanced,
                                                                               n.genes.sample = 1800, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 100, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_100_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_100_unbalanced,
                                                                               n.genes.sample = 1900, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 100, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)),
  # 500 cells
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_500_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_500_balanced,
                                                                             n.genes.sample = 3100, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.1, Cells = 500, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_500_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_500_balanced,
                                                                             n.genes.sample = 2600, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.2, Cells = 500, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_500_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_500_unbalanced,
                                                                               n.genes.sample = 2200, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 500, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_500_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_500_unbalanced,
                                                                               n.genes.sample = 1500, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 500, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)), 
  # 1000 cells
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_1000_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_1000_balanced,
                                                                              n.genes.sample = 1400, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 1000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_1000_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_1000_balanced,
                                                                              n.genes.sample = 2000, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 1000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_1000_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_1000_unbalanced,
                                                                                n.genes.sample = 3100, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 1000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_1000_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_1000_unbalanced,
                                                                                n.genes.sample = 1900, 
                                                                                n.iter = 1,
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 1000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  # 2500 cells
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_2500_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_2500_balanced,
                                                                              n.genes.sample = 1900, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 2500, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_2500_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_2500_balanced,
                                                                              n.genes.sample = 1700, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 2500, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_2500_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_2500_unbalanced,
                                                                                n.genes.sample = 2800, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 2500, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_2500_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_2500_unbalanced,
                                                                                n.genes.sample = 1600, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 2500, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  # 5000 cells
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_5000_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_5000_balanced,
                                                                              n.genes.sample = 2200, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 5000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_5000_balanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_5000_balanced,
                                                                              n.genes.sample = 1300, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 5000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_10_CELLS_5000_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_10_CELLS_5000_unbalanced,
                                                                                n.genes.sample = 1800, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 5000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_panc_DEG_20_CELLS_5000_unbalanced, run_tradeSeq_multi(sim.data = panc_sim_DEG_20_CELLS_5000_unbalanced,
                                                                                n.genes.sample = 2500, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 5000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  
  ##### ENDOCRINOGENESIS #####
  ### MULTI-SUBJECT 
  # 100 cells
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_100_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_100_balanced,
                                                                             n.genes.sample = 2100, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.1, Cells = 100, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_100_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_100_balanced,
                                                                             n.genes.sample = 1500, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.2, Cells = 100, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_100_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_100_unbalanced,
                                                                               n.genes.sample = 1700, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 100, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_100_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_100_unbalanced,
                                                                               n.genes.sample = 2400, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 100, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)),
  # 500 cells
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_500_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_500_balanced,
                                                                             n.genes.sample = 1900, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.1, Cells = 500, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_500_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_500_balanced,
                                                                             n.genes.sample = 2000, 
                                                                             param.list = list(Prop_Dyn_Genes = 0.2, Cells = 500, Allocation = "balanced", Method = "GAM"),
                                                                             n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_500_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_500_unbalanced,
                                                                               n.genes.sample = 1400, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.1, Cells = 500, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_500_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_500_unbalanced,
                                                                               n.genes.sample = 1800, 
                                                                               param.list = list(Prop_Dyn_Genes = 0.2, Cells = 500, Allocation = "unbalanced", Method = "GAM"),
                                                                               n.cores = 4)), 
  # 1000 cells
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_1000_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_1000_balanced,
                                                                              n.genes.sample = 1400, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 1000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_1000_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_1000_balanced,
                                                                              n.genes.sample = 2200, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 1000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_1000_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_1000_unbalanced,
                                                                                n.genes.sample = 2800, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 1000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_1000_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_1000_unbalanced,
                                                                                n.genes.sample = 1600, 
                                                                                n.iter = 1,
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 1000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  # 2500 cells
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_2500_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_2500_balanced,
                                                                              n.genes.sample = 2200, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 2500, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_2500_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_2500_balanced,
                                                                              n.genes.sample = 1700, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 2500, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_2500_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_2500_unbalanced,
                                                                                n.genes.sample = 1600, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 2500, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_2500_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_2500_unbalanced,
                                                                                n.genes.sample = 2600, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 2500, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  # 5000 cells
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_5000_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_5000_balanced,
                                                                              n.genes.sample = 1800, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.1, Cells = 5000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_5000_balanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_5000_balanced,
                                                                              n.genes.sample = 1500, 
                                                                              param.list = list(Prop_Dyn_Genes = 0.2, Cells = 5000, Allocation = "balanced", Method = "GAM"),
                                                                              n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_10_CELLS_5000_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_10_CELLS_5000_unbalanced,
                                                                                n.genes.sample = 2000, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.1, Cells = 5000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)),
  tar_target(tradeSeq_res_endo_DEG_20_CELLS_5000_unbalanced, run_tradeSeq_multi(sim.data = endo_sim_DEG_20_CELLS_5000_unbalanced,
                                                                                n.genes.sample = 2400, 
                                                                                param.list = list(Prop_Dyn_Genes = 0.2, Cells = 5000, Allocation = "unbalanced", Method = "GAM"),
                                                                                n.cores = 4)), 
  
  ##### QC #####
  tar_render(QC_report, "Reports/Model_Results_tradeSeq_Multi.Rmd")
)
