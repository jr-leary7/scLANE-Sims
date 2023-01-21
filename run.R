###### CONFIG #####

# tar_config_set(script = "script_simulation.R", 
#                store = "store_simulation", 
#                project = "simulation", 
#                reporter_make = "verbose",
#                shortcut = FALSE)
# tar_config_set(script = "script_scLANE_GLM_single.R",
#                store = "store_scLANE_GLM_single",
#                project = "scLANE_GLM_models_single",
#                inherits = "simulation")
# tar_config_set(script = "script_scLANE_GLM_multi.R",
#                store = "store_scLANE_GLM_multi",
#                project = "scLANE_GLM_models_multi",
#                inherits = "simulation")
# tar_config_set(script = "script_scLANE_GEE_multi.R", 
#                store = "store_scLANE_GEE_multi", 
#                project = "scLANE_GEE_models_multi",
#                inherits = "simulation")
# tar_config_set(script = "script_scLANE_GLMM_multi.R", 
#                store = "store_scLANE_GLMM_multi", 
#                project = "scLANE_GLMM_models_multi",
#                inherits = "simulation")
# tar_config_set(script = "script_tradeSeq_single.R",
#                store = "store_tradeSeq_single",
#                project = "tradeSeq_models_single",
#                inherits = "simulation")
# tar_config_set(script = "script_tradeSeq_multi.R",
#                store = "store_tradeSeq_multi",
#                project = "tradeSeq_models_multi",
#                inherits = "simulation")
# tar_config_set(script = "script_Lamian_single.R",
#                store = "store_Lamian_single",
#                project = "Lamian_models_single",
#                inherits = "simulation")
# tar_config_set(script = "script_Lamian_multi.R",
#                store = "store_Lamian_multi",
#                project = "Lamian_models_multi",
#                inherits = "simulation")
# tar_config_set(script = "script_monocle3_single.R",
#                store = "store_monocle3_single",
#                project = "monocle3_models_single",
#                inherits = "simulation")
# tar_config_set(script = "script_monocle3_multi.R",
#                store = "store_monocle3_multi",
#                project = "monocle3_models_multi",
#                inherits = "simulation")

##### BATCH JOB #####

# sbatch -t 80:00:00 -c 2 --mem=250G -J scLANE_sim --account=biostat-dept --qos=biostat-dept-b --wrap="module load R; Rscript run.R"

##### PIPELINES #####

# setup
setwd("/blue/rbacher/j.leary/repos/scLANE-Sims/")
library(targets)
library(tarchetypes)

# generate simulated datasets
Sys.setenv(TAR_PROJECT = "simulation")
tar_make_future(workers = 6)

# run tradeSeq -- single-subject
Sys.setenv(TAR_PROJECT = "tradeSeq_models_single")
tar_make_future(workers = 4)

# run tradeSeq -- multi-subject
Sys.setenv(TAR_PROJECT = "tradeSeq_models_multi")
tar_make_future(workers = 4)

# run Lamian -- multi-subject
Sys.setenv(TAR_PROJECT = "Lamian_models_multi")
tar_make_future(workers = 4)

# run scLANE (GLM backend) -- single-subject
Sys.setenv(TAR_PROJECT = "scLANE_GLM_models_single")
tar_make_future(workers = 4)

# run scLANE (GLM backend) -- multi-subject
Sys.setenv(TAR_PROJECT = "scLANE_GLM_models_multi")
tar_make_future(workers = 4)
