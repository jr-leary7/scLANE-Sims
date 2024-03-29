##### setup #####

library(future)
library(targets)
library(tarchetypes)
library(future.callr)

future::plan(future.callr::callr)
options(future.globals.maxSize = 100000 * 1024^2)  # 100GB max size -- necessary for loading in all sim datasets

source("./R/functions_scLANE_GLM.R")

tar_option_set(packages = c("qs", 
                            "glm2", 
                            "geeM", 
                            "mgcv", 
                            "pryr", 
                            "MASS", 
                            "tidyr", 
                            "stats", 
                            "broom", 
                            "scran", 
                            "gamlss", 
                            "igraph",
                            "scater", 
                            "scLANE", 
                            "foreach",  
                            "scRNAseq", 
                            "parallel", 
                            "bigstatsr", 
                            "S4Vectors", 
                            "tidyverse", 
                            "doParallel", 
                            "SummarizedExperiment", 
                            "SingleCellExperiment"), 
               imports = c("scLANE"), 
               error = "continue", 
               memory = "transient",
               retrieval = "worker", 
               storage = "worker", 
               garbage_collection = TRUE, 
               format = "qs")

# targets 
list(
  ##### UPSTREAM TARGETS #####
  ##### BRAIN #####
  ### FILES ###
  # 100 cells
  tar_target(brain_sim_DEG_01_CELLS_100_file, "store_simulation/objects/brain_sim_DEG_01_CELLS_100", format = "file"),
  tar_target(brain_sim_DEG_05_CELLS_100_file, "store_simulation/objects/brain_sim_DEG_05_CELLS_100", format = "file"),
  tar_target(brain_sim_DEG_10_CELLS_100_file, "store_simulation/objects/brain_sim_DEG_10_CELLS_100", format = "file"),
  tar_target(brain_sim_DEG_20_CELLS_100_file, "store_simulation/objects/brain_sim_DEG_20_CELLS_100", format = "file"),
  # 500 cells
  tar_target(brain_sim_DEG_01_CELLS_500_file, "store_simulation/objects/brain_sim_DEG_01_CELLS_500", format = "file"),
  tar_target(brain_sim_DEG_05_CELLS_500_file, "store_simulation/objects/brain_sim_DEG_05_CELLS_500", format = "file"),
  tar_target(brain_sim_DEG_10_CELLS_500_file, "store_simulation/objects/brain_sim_DEG_10_CELLS_500", format = "file"),
  tar_target(brain_sim_DEG_20_CELLS_500_file, "store_simulation/objects/brain_sim_DEG_20_CELLS_500", format = "file"),
  # 1000 cells
  tar_target(brain_sim_DEG_01_CELLS_1000_file, "store_simulation/objects/brain_sim_DEG_01_CELLS_1000", format = "file"),
  tar_target(brain_sim_DEG_05_CELLS_1000_file, "store_simulation/objects/brain_sim_DEG_05_CELLS_1000", format = "file"),
  tar_target(brain_sim_DEG_10_CELLS_1000_file, "store_simulation/objects/brain_sim_DEG_10_CELLS_1000", format = "file"),
  tar_target(brain_sim_DEG_20_CELLS_1000_file, "store_simulation/objects/brain_sim_DEG_20_CELLS_1000", format = "file"),
  # 2500 cells
  tar_target(brain_sim_DEG_01_CELLS_2500_file, "store_simulation/objects/brain_sim_DEG_01_CELLS_2500", format = "file"),
  tar_target(brain_sim_DEG_05_CELLS_2500_file, "store_simulation/objects/brain_sim_DEG_05_CELLS_2500", format = "file"),
  tar_target(brain_sim_DEG_10_CELLS_2500_file, "store_simulation/objects/brain_sim_DEG_10_CELLS_2500", format = "file"),
  tar_target(brain_sim_DEG_20_CELLS_2500_file, "store_simulation/objects/brain_sim_DEG_20_CELLS_2500", format = "file"),
  # 5000 cells
  tar_target(brain_sim_DEG_01_CELLS_5000_file, "store_simulation/objects/brain_sim_DEG_01_CELLS_5000", format = "file"),
  tar_target(brain_sim_DEG_05_CELLS_5000_file, "store_simulation/objects/brain_sim_DEG_05_CELLS_5000", format = "file"),
  tar_target(brain_sim_DEG_10_CELLS_5000_file, "store_simulation/objects/brain_sim_DEG_10_CELLS_5000", format = "file"),
  tar_target(brain_sim_DEG_20_CELLS_5000_file, "store_simulation/objects/brain_sim_DEG_20_CELLS_5000", format = "file"),
  
  ### OBJECTS ###
  # 100 cells 
  tar_target(brain_sim_DEG_01_CELLS_100, qs::qread(brain_sim_DEG_01_CELLS_100_file)), 
  tar_target(brain_sim_DEG_05_CELLS_100, qs::qread(brain_sim_DEG_05_CELLS_100_file)), 
  tar_target(brain_sim_DEG_10_CELLS_100, qs::qread(brain_sim_DEG_10_CELLS_100_file)), 
  tar_target(brain_sim_DEG_20_CELLS_100, qs::qread(brain_sim_DEG_20_CELLS_100_file)), 
  # 500 cells 
  tar_target(brain_sim_DEG_01_CELLS_500, qs::qread(brain_sim_DEG_01_CELLS_500_file)), 
  tar_target(brain_sim_DEG_05_CELLS_500, qs::qread(brain_sim_DEG_05_CELLS_500_file)), 
  tar_target(brain_sim_DEG_10_CELLS_500, qs::qread(brain_sim_DEG_10_CELLS_500_file)), 
  tar_target(brain_sim_DEG_20_CELLS_500, qs::qread(brain_sim_DEG_20_CELLS_500_file)), 
  # 1000 cells 
  tar_target(brain_sim_DEG_01_CELLS_1000, qs::qread(brain_sim_DEG_01_CELLS_1000_file)), 
  tar_target(brain_sim_DEG_05_CELLS_1000, qs::qread(brain_sim_DEG_05_CELLS_1000_file)), 
  tar_target(brain_sim_DEG_10_CELLS_1000, qs::qread(brain_sim_DEG_10_CELLS_1000_file)), 
  tar_target(brain_sim_DEG_20_CELLS_1000, qs::qread(brain_sim_DEG_20_CELLS_1000_file)), 
  # 2500 cells
  tar_target(brain_sim_DEG_01_CELLS_2500, qs::qread(brain_sim_DEG_01_CELLS_2500_file)), 
  tar_target(brain_sim_DEG_05_CELLS_2500, qs::qread(brain_sim_DEG_05_CELLS_2500_file)), 
  tar_target(brain_sim_DEG_10_CELLS_2500, qs::qread(brain_sim_DEG_10_CELLS_2500_file)), 
  tar_target(brain_sim_DEG_20_CELLS_2500, qs::qread(brain_sim_DEG_20_CELLS_2500_file)), 
  # 5000 cells 
  tar_target(brain_sim_DEG_01_CELLS_5000, qs::qread(brain_sim_DEG_01_CELLS_5000_file)), 
  tar_target(brain_sim_DEG_05_CELLS_5000, qs::qread(brain_sim_DEG_05_CELLS_5000_file)), 
  tar_target(brain_sim_DEG_10_CELLS_5000, qs::qread(brain_sim_DEG_10_CELLS_5000_file)), 
  tar_target(brain_sim_DEG_20_CELLS_5000, qs::qread(brain_sim_DEG_20_CELLS_5000_file)), 
  
  ##### PANCREAS #####
  ### FILES ###
  # 100 cells
  tar_target(panc_sim_DEG_01_CELLS_100_file, "store_simulation/objects/panc_sim_DEG_01_CELLS_100", format = "file"),
  tar_target(panc_sim_DEG_05_CELLS_100_file, "store_simulation/objects/panc_sim_DEG_05_CELLS_100", format = "file"),
  tar_target(panc_sim_DEG_10_CELLS_100_file, "store_simulation/objects/panc_sim_DEG_10_CELLS_100", format = "file"),
  tar_target(panc_sim_DEG_20_CELLS_100_file, "store_simulation/objects/panc_sim_DEG_20_CELLS_100", format = "file"),
  # 500 cells
  tar_target(panc_sim_DEG_01_CELLS_500_file, "store_simulation/objects/panc_sim_DEG_01_CELLS_500", format = "file"),
  tar_target(panc_sim_DEG_05_CELLS_500_file, "store_simulation/objects/panc_sim_DEG_05_CELLS_500", format = "file"),
  tar_target(panc_sim_DEG_10_CELLS_500_file, "store_simulation/objects/panc_sim_DEG_10_CELLS_500", format = "file"),
  tar_target(panc_sim_DEG_20_CELLS_500_file, "store_simulation/objects/panc_sim_DEG_20_CELLS_500", format = "file"),
  # 1000 cells
  tar_target(panc_sim_DEG_01_CELLS_1000_file, "store_simulation/objects/panc_sim_DEG_01_CELLS_1000", format = "file"),
  tar_target(panc_sim_DEG_05_CELLS_1000_file, "store_simulation/objects/panc_sim_DEG_05_CELLS_1000", format = "file"),
  tar_target(panc_sim_DEG_10_CELLS_1000_file, "store_simulation/objects/panc_sim_DEG_10_CELLS_1000", format = "file"),
  tar_target(panc_sim_DEG_20_CELLS_1000_file, "store_simulation/objects/panc_sim_DEG_20_CELLS_1000", format = "file"),
  # 2500 cells
  tar_target(panc_sim_DEG_01_CELLS_2500_file, "store_simulation/objects/panc_sim_DEG_01_CELLS_2500", format = "file"),
  tar_target(panc_sim_DEG_05_CELLS_2500_file, "store_simulation/objects/panc_sim_DEG_05_CELLS_2500", format = "file"),
  tar_target(panc_sim_DEG_10_CELLS_2500_file, "store_simulation/objects/panc_sim_DEG_10_CELLS_2500", format = "file"),
  tar_target(panc_sim_DEG_20_CELLS_2500_file, "store_simulation/objects/panc_sim_DEG_20_CELLS_2500", format = "file"),
  # 5000 cells
  tar_target(panc_sim_DEG_01_CELLS_5000_file, "store_simulation/objects/panc_sim_DEG_01_CELLS_5000", format = "file"),
  tar_target(panc_sim_DEG_05_CELLS_5000_file, "store_simulation/objects/panc_sim_DEG_05_CELLS_5000", format = "file"),
  tar_target(panc_sim_DEG_10_CELLS_5000_file, "store_simulation/objects/panc_sim_DEG_10_CELLS_5000", format = "file"),
  tar_target(panc_sim_DEG_20_CELLS_5000_file, "store_simulation/objects/panc_sim_DEG_20_CELLS_5000", format = "file"),
  
  ### OBJECTS ###
  # 100 cells 
  tar_target(panc_sim_DEG_01_CELLS_100, qs::qread(panc_sim_DEG_01_CELLS_100_file)), 
  tar_target(panc_sim_DEG_05_CELLS_100, qs::qread(panc_sim_DEG_05_CELLS_100_file)), 
  tar_target(panc_sim_DEG_10_CELLS_100, qs::qread(panc_sim_DEG_10_CELLS_100_file)), 
  tar_target(panc_sim_DEG_20_CELLS_100, qs::qread(panc_sim_DEG_20_CELLS_100_file)), 
  # 500 cells 
  tar_target(panc_sim_DEG_01_CELLS_500, qs::qread(panc_sim_DEG_01_CELLS_500_file)), 
  tar_target(panc_sim_DEG_05_CELLS_500, qs::qread(panc_sim_DEG_05_CELLS_500_file)), 
  tar_target(panc_sim_DEG_10_CELLS_500, qs::qread(panc_sim_DEG_10_CELLS_500_file)), 
  tar_target(panc_sim_DEG_20_CELLS_500, qs::qread(panc_sim_DEG_20_CELLS_500_file)), 
  # 1000 cells 
  tar_target(panc_sim_DEG_01_CELLS_1000, qs::qread(panc_sim_DEG_01_CELLS_1000_file)), 
  tar_target(panc_sim_DEG_05_CELLS_1000, qs::qread(panc_sim_DEG_05_CELLS_1000_file)), 
  tar_target(panc_sim_DEG_10_CELLS_1000, qs::qread(panc_sim_DEG_10_CELLS_1000_file)), 
  tar_target(panc_sim_DEG_20_CELLS_1000, qs::qread(panc_sim_DEG_20_CELLS_1000_file)), 
  # 2500 cells
  tar_target(panc_sim_DEG_01_CELLS_2500, qs::qread(panc_sim_DEG_01_CELLS_2500_file)), 
  tar_target(panc_sim_DEG_05_CELLS_2500, qs::qread(panc_sim_DEG_05_CELLS_2500_file)), 
  tar_target(panc_sim_DEG_10_CELLS_2500, qs::qread(panc_sim_DEG_10_CELLS_2500_file)), 
  tar_target(panc_sim_DEG_20_CELLS_2500, qs::qread(panc_sim_DEG_20_CELLS_2500_file)), 
  # 5000 cells 
  tar_target(panc_sim_DEG_01_CELLS_5000, qs::qread(panc_sim_DEG_01_CELLS_5000_file)), 
  tar_target(panc_sim_DEG_05_CELLS_5000, qs::qread(panc_sim_DEG_05_CELLS_5000_file)), 
  tar_target(panc_sim_DEG_10_CELLS_5000, qs::qread(panc_sim_DEG_10_CELLS_5000_file)), 
  tar_target(panc_sim_DEG_20_CELLS_5000, qs::qread(panc_sim_DEG_20_CELLS_5000_file)), 
  
  ##### ENDOCRINOGENESIS #####
  ### FILES ###
  # 100 cells
  tar_target(endo_sim_DEG_01_CELLS_100_file, "store_simulation/objects/endo_sim_DEG_01_CELLS_100", format = "file"),
  tar_target(endo_sim_DEG_05_CELLS_100_file, "store_simulation/objects/endo_sim_DEG_05_CELLS_100", format = "file"),
  tar_target(endo_sim_DEG_10_CELLS_100_file, "store_simulation/objects/endo_sim_DEG_10_CELLS_100", format = "file"),
  tar_target(endo_sim_DEG_20_CELLS_100_file, "store_simulation/objects/endo_sim_DEG_20_CELLS_100", format = "file"),
  # 500 cells
  tar_target(endo_sim_DEG_01_CELLS_500_file, "store_simulation/objects/endo_sim_DEG_01_CELLS_500", format = "file"),
  tar_target(endo_sim_DEG_05_CELLS_500_file, "store_simulation/objects/endo_sim_DEG_05_CELLS_500", format = "file"),
  tar_target(endo_sim_DEG_10_CELLS_500_file, "store_simulation/objects/endo_sim_DEG_10_CELLS_500", format = "file"),
  tar_target(endo_sim_DEG_20_CELLS_500_file, "store_simulation/objects/endo_sim_DEG_20_CELLS_500", format = "file"),
  # 1000 cells
  tar_target(endo_sim_DEG_01_CELLS_1000_file, "store_simulation/objects/endo_sim_DEG_01_CELLS_1000", format = "file"),
  tar_target(endo_sim_DEG_05_CELLS_1000_file, "store_simulation/objects/endo_sim_DEG_05_CELLS_1000", format = "file"),
  tar_target(endo_sim_DEG_10_CELLS_1000_file, "store_simulation/objects/endo_sim_DEG_10_CELLS_1000", format = "file"),
  tar_target(endo_sim_DEG_20_CELLS_1000_file, "store_simulation/objects/endo_sim_DEG_20_CELLS_1000", format = "file"),
  # 2500 cells
  tar_target(endo_sim_DEG_01_CELLS_2500_file, "store_simulation/objects/endo_sim_DEG_01_CELLS_2500", format = "file"),
  tar_target(endo_sim_DEG_05_CELLS_2500_file, "store_simulation/objects/endo_sim_DEG_05_CELLS_2500", format = "file"),
  tar_target(endo_sim_DEG_10_CELLS_2500_file, "store_simulation/objects/endo_sim_DEG_10_CELLS_2500", format = "file"),
  tar_target(endo_sim_DEG_20_CELLS_2500_file, "store_simulation/objects/endo_sim_DEG_20_CELLS_2500", format = "file"),
  # 5000 cells
  tar_target(endo_sim_DEG_01_CELLS_5000_file, "store_simulation/objects/endo_sim_DEG_01_CELLS_5000", format = "file"),
  tar_target(endo_sim_DEG_05_CELLS_5000_file, "store_simulation/objects/endo_sim_DEG_05_CELLS_5000", format = "file"),
  tar_target(endo_sim_DEG_10_CELLS_5000_file, "store_simulation/objects/endo_sim_DEG_10_CELLS_5000", format = "file"),
  tar_target(endo_sim_DEG_20_CELLS_5000_file, "store_simulation/objects/endo_sim_DEG_20_CELLS_5000", format = "file"),
  
  ### OBJECTS ###
  # 100 cells 
  tar_target(endo_sim_DEG_01_CELLS_100, qs::qread(endo_sim_DEG_01_CELLS_100_file)), 
  tar_target(endo_sim_DEG_05_CELLS_100, qs::qread(endo_sim_DEG_05_CELLS_100_file)), 
  tar_target(endo_sim_DEG_10_CELLS_100, qs::qread(endo_sim_DEG_10_CELLS_100_file)), 
  tar_target(endo_sim_DEG_20_CELLS_100, qs::qread(endo_sim_DEG_20_CELLS_100_file)), 
  # 500 cells 
  tar_target(endo_sim_DEG_01_CELLS_500, qs::qread(endo_sim_DEG_01_CELLS_500_file)), 
  tar_target(endo_sim_DEG_05_CELLS_500, qs::qread(endo_sim_DEG_05_CELLS_500_file)), 
  tar_target(endo_sim_DEG_10_CELLS_500, qs::qread(endo_sim_DEG_10_CELLS_500_file)), 
  tar_target(endo_sim_DEG_20_CELLS_500, qs::qread(endo_sim_DEG_20_CELLS_500_file)), 
  # 1000 cells 
  tar_target(endo_sim_DEG_01_CELLS_1000, qs::qread(endo_sim_DEG_01_CELLS_1000_file)), 
  tar_target(endo_sim_DEG_05_CELLS_1000, qs::qread(endo_sim_DEG_05_CELLS_1000_file)), 
  tar_target(endo_sim_DEG_10_CELLS_1000, qs::qread(endo_sim_DEG_10_CELLS_1000_file)), 
  tar_target(endo_sim_DEG_20_CELLS_1000, qs::qread(endo_sim_DEG_20_CELLS_1000_file)), 
  # 2500 cells
  tar_target(endo_sim_DEG_01_CELLS_2500, qs::qread(endo_sim_DEG_01_CELLS_2500_file)), 
  tar_target(endo_sim_DEG_05_CELLS_2500, qs::qread(endo_sim_DEG_05_CELLS_2500_file)), 
  tar_target(endo_sim_DEG_10_CELLS_2500, qs::qread(endo_sim_DEG_10_CELLS_2500_file)), 
  tar_target(endo_sim_DEG_20_CELLS_2500, qs::qread(endo_sim_DEG_20_CELLS_2500_file)), 
  # 5000 cells 
  tar_target(endo_sim_DEG_01_CELLS_5000, qs::qread(endo_sim_DEG_01_CELLS_5000_file)), 
  tar_target(endo_sim_DEG_05_CELLS_5000, qs::qread(endo_sim_DEG_05_CELLS_5000_file)), 
  tar_target(endo_sim_DEG_10_CELLS_5000, qs::qread(endo_sim_DEG_10_CELLS_5000_file)), 
  tar_target(endo_sim_DEG_20_CELLS_5000, qs::qread(endo_sim_DEG_20_CELLS_5000_file)), 
  
  ##### MODELS #####
  ##### BRAIN #####
  ### SINGLE-SUBJECT
  # 100 cells
  tar_target(scLANE_GLM_brain_DEG_01_CELLS_100, run_scLANE_GLM(sim.data = brain_sim_DEG_01_CELLS_100, 
                                                               n.genes.sample = 3000, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 100, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_05_CELLS_100, run_scLANE_GLM(sim.data = brain_sim_DEG_05_CELLS_100, 
                                                               n.genes.sample = 4000, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 100, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_10_CELLS_100, run_scLANE_GLM(sim.data = brain_sim_DEG_10_CELLS_100,
                                                               n.genes.sample = 2000, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 100, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_20_CELLS_100, run_scLANE_GLM(sim.data = brain_sim_DEG_20_CELLS_100,
                                                               n.genes.sample = 3500, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 100, Method = "GLM"),
                                                               n.cores = 4)),
  # 500 cells
  tar_target(scLANE_GLM_brain_DEG_01_CELLS_500, run_scLANE_GLM(sim.data = brain_sim_DEG_01_CELLS_500,
                                                               n.genes.sample = 2700, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_05_CELLS_500, run_scLANE_GLM(sim.data = brain_sim_DEG_05_CELLS_500,
                                                               n.genes.sample = 3200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_10_CELLS_500, run_scLANE_GLM(sim.data = brain_sim_DEG_10_CELLS_500,
                                                               n.genes.sample = 4200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_20_CELLS_500, run_scLANE_GLM(sim.data = brain_sim_DEG_20_CELLS_500,
                                                               n.genes.sample = 1000, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 500, Method = "GLM"),
                                                               n.cores = 4)),
  # 1000 cells
  tar_target(scLANE_GLM_brain_DEG_01_CELLS_1000, run_scLANE_GLM(sim.data = brain_sim_DEG_01_CELLS_1000,
                                                                n.genes.sample = 1100, 
                                                                param.list = list(Prop_Dyn_Genes = 0.01, Cells = 1000, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_05_CELLS_1000, run_scLANE_GLM(sim.data = brain_sim_DEG_05_CELLS_1000,
                                                                n.genes.sample = 3600, 
                                                                param.list = list(Prop_Dyn_Genes = 0.05, Cells = 1000, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_10_CELLS_1000, run_scLANE_GLM(sim.data = brain_sim_DEG_10_CELLS_1000,
                                                                n.genes.sample = 2700, 
                                                                param.list = list(Prop_Dyn_Genes = 0.10, Cells = 1000, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_20_CELLS_1000, run_scLANE_GLM(sim.data = brain_sim_DEG_20_CELLS_1000,
                                                                n.genes.sample = 1400, 
                                                                param.list = list(Prop_Dyn_Genes = 0.20, Cells = 1000, Method = "GLM"),
                                                                n.cores = 4)),
  # 2500 cells
  tar_target(scLANE_GLM_brain_DEG_01_CELLS_2500, run_scLANE_GLM(sim.data = brain_sim_DEG_01_CELLS_2500,
                                                                n.genes.sample = 2300, 
                                                                param.list = list(Prop_Dyn_Genes = 0.01, Cells = 2500, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_05_CELLS_2500, run_scLANE_GLM(sim.data = brain_sim_DEG_05_CELLS_2500,
                                                                n.genes.sample = 1000, 
                                                                param.list = list(Prop_Dyn_Genes = 0.05, Cells = 2500, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_10_CELLS_2500, run_scLANE_GLM(sim.data = brain_sim_DEG_10_CELLS_2500,
                                                                n.genes.sample = 2400, 
                                                                param.list = list(Prop_Dyn_Genes = 0.10, Cells = 2500, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_20_CELLS_2500, run_scLANE_GLM(sim.data = brain_sim_DEG_20_CELLS_2500,
                                                                n.genes.sample = 3200, 
                                                                param.list = list(Prop_Dyn_Genes = 0.20, Cells = 2500, Method = "GLM"),
                                                                n.cores = 4)),
  # 5000 cells
  tar_target(scLANE_GLM_brain_DEG_01_CELLS_5000, run_scLANE_GLM(sim.data = brain_sim_DEG_01_CELLS_5000,
                                                                n.genes.sample = 3500, 
                                                                param.list = list(Prop_Dyn_Genes = 0.01, Cells = 5000, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_05_CELLS_5000, run_scLANE_GLM(sim.data = brain_sim_DEG_05_CELLS_5000,
                                                                n.genes.sample = 1100, 
                                                                param.list = list(Prop_Dyn_Genes = 0.05, Cells = 5000, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_10_CELLS_5000, run_scLANE_GLM(sim.data = brain_sim_DEG_10_CELLS_5000,
                                                                n.genes.sample = 3100, 
                                                                param.list = list(Prop_Dyn_Genes = 0.10, Cells = 5000, Method = "GLM"),
                                                                n.cores = 4)),
  tar_target(scLANE_GLM_brain_DEG_20_CELLS_5000, run_scLANE_GLM(sim.data = brain_sim_DEG_20_CELLS_5000,
                                                                n.genes.sample = 2200, 
                                                                param.list = list(Prop_Dyn_Genes = 0.20, Cells = 5000, Method = "GLM"),
                                                                n.cores = 4)), 
  
  ##### PANCREAS #####
  ### SINGLE-SUBJECT
  # 100 cells
  tar_target(scLANE_GLM_panc_DEG_01_CELLS_100, run_scLANE_GLM(sim.data = panc_sim_DEG_01_CELLS_100, 
                                                              n.genes.sample = 2100, 
                                                              param.list = list(Prop_Dyn_Genes = 0.01, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_05_CELLS_100, run_scLANE_GLM(sim.data = panc_sim_DEG_05_CELLS_100, 
                                                              n.genes.sample = 3500, 
                                                              param.list = list(Prop_Dyn_Genes = 0.05, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_10_CELLS_100, run_scLANE_GLM(sim.data = panc_sim_DEG_10_CELLS_100,
                                                              n.genes.sample = 1900, 
                                                              param.list = list(Prop_Dyn_Genes = 0.10, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_20_CELLS_100, run_scLANE_GLM(sim.data = panc_sim_DEG_20_CELLS_100,
                                                              n.genes.sample = 3300, 
                                                              param.list = list(Prop_Dyn_Genes = 0.20, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  # 500 cells
  tar_target(scLANE_GLM_panc_DEG_01_CELLS_500, run_scLANE_GLM(sim.data = panc_sim_DEG_01_CELLS_500,
                                                              n.genes.sample = 2700, 
                                                              param.list = list(Prop_Dyn_Genes = 0.01, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_05_CELLS_500, run_scLANE_GLM(sim.data = panc_sim_DEG_05_CELLS_500,
                                                              n.genes.sample = 3200, 
                                                              param.list = list(Prop_Dyn_Genes = 0.05, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_10_CELLS_500, run_scLANE_GLM(sim.data = panc_sim_DEG_10_CELLS_500,
                                                              n.genes.sample = 3300, 
                                                              param.list = list(Prop_Dyn_Genes = 0.10, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_20_CELLS_500, run_scLANE_GLM(sim.data = panc_sim_DEG_20_CELLS_500,
                                                              n.genes.sample = 1200, 
                                                              param.list = list(Prop_Dyn_Genes = 0.20, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  # 1000 cells
  tar_target(scLANE_GLM_panc_DEG_01_CELLS_1000, run_scLANE_GLM(sim.data = panc_sim_DEG_01_CELLS_1000,
                                                               n.genes.sample = 1100, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_05_CELLS_1000, run_scLANE_GLM(sim.data = panc_sim_DEG_05_CELLS_1000,
                                                               n.genes.sample = 3600, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_10_CELLS_1000, run_scLANE_GLM(sim.data = panc_sim_DEG_10_CELLS_1000,
                                                               n.genes.sample = 2200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_20_CELLS_1000, run_scLANE_GLM(sim.data = panc_sim_DEG_20_CELLS_1000,
                                                               n.genes.sample = 1400, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  # 2500 cells
  tar_target(scLANE_GLM_panc_DEG_01_CELLS_2500, run_scLANE_GLM(sim.data = panc_sim_DEG_01_CELLS_2500,
                                                               n.genes.sample = 2400, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_05_CELLS_2500, run_scLANE_GLM(sim.data = panc_sim_DEG_05_CELLS_2500,
                                                               n.genes.sample = 1100, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_10_CELLS_2500, run_scLANE_GLM(sim.data = panc_sim_DEG_10_CELLS_2500,
                                                               n.genes.sample = 2400, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_20_CELLS_2500, run_scLANE_GLM(sim.data = panc_sim_DEG_20_CELLS_2500,
                                                               n.genes.sample = 3200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  # 5000 cells 
  tar_target(scLANE_GLM_panc_DEG_01_CELLS_5000, run_scLANE_GLM(sim.data = panc_sim_DEG_01_CELLS_5000,
                                                               n.genes.sample = 3100, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_05_CELLS_5000, run_scLANE_GLM(sim.data = panc_sim_DEG_05_CELLS_5000,
                                                               n.genes.sample = 1100, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_10_CELLS_5000, run_scLANE_GLM(sim.data = panc_sim_DEG_10_CELLS_5000,
                                                               n.genes.sample = 2600, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_panc_DEG_20_CELLS_5000, run_scLANE_GLM(sim.data = panc_sim_DEG_20_CELLS_5000,
                                                               n.genes.sample = 2200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)), 
  
  ##### ENDOCRINOGENESIS #####
  ### SINGLE-SUBJECT
  # 100 cells
  tar_target(scLANE_GLM_endo_DEG_01_CELLS_100, run_scLANE_GLM(sim.data = endo_sim_DEG_01_CELLS_100, 
                                                              n.genes.sample = 1100, 
                                                              param.list = list(Prop_Dyn_Genes = 0.01, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_05_CELLS_100, run_scLANE_GLM(sim.data = endo_sim_DEG_05_CELLS_100, 
                                                              n.genes.sample = 2200, 
                                                              param.list = list(Prop_Dyn_Genes = 0.05, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_10_CELLS_100, run_scLANE_GLM(sim.data = endo_sim_DEG_10_CELLS_100,
                                                              n.genes.sample = 1900, 
                                                              param.list = list(Prop_Dyn_Genes = 0.10, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_20_CELLS_100, run_scLANE_GLM(sim.data = endo_sim_DEG_20_CELLS_100,
                                                              n.genes.sample = 2400, 
                                                              param.list = list(Prop_Dyn_Genes = 0.20, Cells = 100, Method = "GLM"),
                                                              n.cores = 4)),
  # 500 cells
  tar_target(scLANE_GLM_endo_DEG_01_CELLS_500, run_scLANE_GLM(sim.data = endo_sim_DEG_01_CELLS_500,
                                                              n.genes.sample = 2700, 
                                                              param.list = list(Prop_Dyn_Genes = 0.01, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_05_CELLS_500, run_scLANE_GLM(sim.data = endo_sim_DEG_05_CELLS_500,
                                                              n.genes.sample = 1000, 
                                                              param.list = list(Prop_Dyn_Genes = 0.05, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_10_CELLS_500, run_scLANE_GLM(sim.data = endo_sim_DEG_10_CELLS_500,
                                                              n.genes.sample = 2500, 
                                                              param.list = list(Prop_Dyn_Genes = 0.10, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_20_CELLS_500, run_scLANE_GLM(sim.data = endo_sim_DEG_20_CELLS_500,
                                                              n.genes.sample = 2000, 
                                                              param.list = list(Prop_Dyn_Genes = 0.20, Cells = 500, Method = "GLM"),
                                                              n.cores = 4)),
  # 1000 cells
  tar_target(scLANE_GLM_endo_DEG_01_CELLS_1000, run_scLANE_GLM(sim.data = endo_sim_DEG_01_CELLS_1000,
                                                               n.genes.sample = 1100, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_05_CELLS_1000, run_scLANE_GLM(sim.data = endo_sim_DEG_05_CELLS_1000,
                                                               n.genes.sample = 3600, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_10_CELLS_1000, run_scLANE_GLM(sim.data = endo_sim_DEG_10_CELLS_1000,
                                                               n.genes.sample = 2200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_20_CELLS_1000, run_scLANE_GLM(sim.data = endo_sim_DEG_20_CELLS_1000,
                                                               n.genes.sample = 1400, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 1000, Method = "GLM"),
                                                               n.cores = 4)),
  # 2500 cells
  tar_target(scLANE_GLM_endo_DEG_01_CELLS_2500, run_scLANE_GLM(sim.data = endo_sim_DEG_01_CELLS_2500,
                                                               n.genes.sample = 1700, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_05_CELLS_2500, run_scLANE_GLM(sim.data = endo_sim_DEG_05_CELLS_2500,
                                                               n.genes.sample = 2200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_10_CELLS_2500, run_scLANE_GLM(sim.data = endo_sim_DEG_10_CELLS_2500,
                                                               n.genes.sample = 3000, 
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_20_CELLS_2500, run_scLANE_GLM(sim.data = endo_sim_DEG_20_CELLS_2500,
                                                               n.genes.sample = 3200, 
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 2500, Method = "GLM"),
                                                               n.cores = 4)),
  # 5000 cells 
  tar_target(scLANE_GLM_endo_DEG_01_CELLS_5000, run_scLANE_GLM(sim.data = endo_sim_DEG_01_CELLS_5000,
                                                               n.genes.sample = 1400, 
                                                               param.list = list(Prop_Dyn_Genes = 0.01, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_05_CELLS_5000, run_scLANE_GLM(sim.data = endo_sim_DEG_05_CELLS_5000,
                                                               n.genes.sample = 1500, 
                                                               param.list = list(Prop_Dyn_Genes = 0.05, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_10_CELLS_5000, run_scLANE_GLM(sim.data = endo_sim_DEG_10_CELLS_5000,
                                                               n.genes.sample = 2100,
                                                               param.list = list(Prop_Dyn_Genes = 0.10, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)),
  tar_target(scLANE_GLM_endo_DEG_20_CELLS_5000, run_scLANE_GLM(sim.data = endo_sim_DEG_20_CELLS_5000,
                                                               n.genes.sample = 3000,
                                                               param.list = list(Prop_Dyn_Genes = 0.20, Cells = 5000, Method = "GLM"),
                                                               n.cores = 4)), 
  
  ##### QC #####
  tar_render(QC_report, "Reports/Model_Results_scLANE_GLM_Single.Rmd")
)
