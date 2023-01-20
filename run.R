# set up config -- only run once 
# tar_config_set(script = "script_simulation.R", 
#                store = "store_simulation", 
#                project = "simulation", 
#                reporter_make = "verbose",
#                shortcut = FALSE)
# tar_config_set(script = "script_scLANE_GLM.R", 
#                store = "store_scLANE_GLM", 
#                project = "scLANE_GLM_models",
#                inherits = "simulation")
# tar_config_set(script = "script_scLANE_GEE.R", 
#                store = "store_scLANE_GEE", 
#                project = "scLANE_GEE_models",
#                inherits = "simulation")
# tar_config_set(script = "script_scLANE_GLMM.R", 
#                store = "store_scLANE_GLMM", 
#                project = "scLANE_GLMM_models",
#                inherits = "simulation")
# tar_config_set(script = "script_tradeSeq.R",
#                store = "store_tradeSeq",
#                project = "tradeSeq_models",
#                inherits = "simulation")
# tar_config_set(script = "script_Lamian.R",
#                store = "store_Lamian",
#                project = "Lamian_models",
#                inherits = "simulation")

# command: sbatch -t 80:00:00 -c 2 --mem=250G -J scLANE_sim --account=biostat-dept --qos=biostat-dept-b --wrap="module load R; Rscript run.R"

# setup
setwd("/blue/rbacher/j.leary/repos/scLANE_Analysis/")
library(targets)
library(tarchetypes)

# run simulation 
Sys.setenv(TAR_PROJECT = "simulation")
tar_make_future(workers = 6)

# run tradeSeq
Sys.setenv(TAR_PROJECT = "tradeSeq_models")
tar_make_future(workers = 6)

# run Lamian
Sys.setenv(TAR_PROJECT = "Lamian_models")
tar_make_future(workers = 6)

# run scLANE (w/ GLM backend)
Sys.setenv(TAR_PROJECT = "scLANE_GLM_models")
tar_make_future(workers = 6)

# run scLANE (w/ GEE backend)
# Sys.setenv(TAR_PROJECT = "scLANE_GEE_models")
# tar_make_future(workers = 6)

# run scLANE (w/ GLMM backend)
# Sys.setenv(TAR_PROJECT = "scLANE_GLMM_models")
# tar_make_future(workers = 6)
