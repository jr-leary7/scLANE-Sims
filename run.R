# command: sbatch -t 80:00:00 -c 6 --mem=180G -J sim --account=biostat-dept --qos=biostat-dept-b --wrap="module load R; Rscript run.R" -- 3 days, 8h total runtime possible

setwd("/blue/rbacher/j.leary/repos/scLANE-Sims")
library(targets)
library(tarchetypes)
tar_make_future(workers = 10)
