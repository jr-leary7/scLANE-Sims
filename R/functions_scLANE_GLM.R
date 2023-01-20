# run GLM scLANE on subsampled single-subject data 
run_scLANE_GLM <- function(sim.data = NULL,
                           n.genes.sample = 1000,
                           n.iter = 1,
                           param.list = NULL,
                           n.cores = 4,
                           scLANE.log = FALSE,
                           scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) || is.null(param.list)) { stop("You failed to provide necessary parameters to run_scLANE_GLM().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe -- {targets} takes care of seed
  p_dynamic <- mean(SummarizedExperiment::rowData(sim.data)[, 1] == "Dynamic")
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- SummarizedExperiment::rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- SummarizedExperiment::rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- SummarizedExperiment::colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  n.potential.basis.fns = 5,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, 
                                            p.adj.method = "holm", 
                                            fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "holm",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                 as.data.frame() %>%
                                                 dplyr::mutate(gene = rownames(.))), 
                                                by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}

# run GLM scLANE on subsampled multi-subject data 
run_scLANE_GLM_multi <- function(sim.data = NULL,
                                 n.genes.sample = 1000,
                                 n.iter = 1,
                                 param.list = NULL,
                                 n.cores = 4,
                                 scLANE.log = FALSE,
                                 scLANE.log.iter = 1000) {
  # check inputs
  if (is.null(sim.data) || is.null(param.list)) { stop("You failed to provide necessary parameters to run_scLANE_GLM_multi().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 9)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "testDynamic_results_raw", "testDynamic_results_tidy", "testSlope_results", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  td_res_raw_list <- vector("list", length = n.iter)
  td_res_tidy_list <- vector("list", length = n.iter)
  ts_res_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe -- {targets} takes care of seed
  p_dynamic <- SummarizedExperiment::rowData(sim.data) %>%
               as.data.frame() %>% 
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
               # dplyr::summarise(across(everything(), \(x) mean(x == "Dynamic")))
               tidyr::pivot_longer(cols = tidyselect::everything(), values_to = "geneStatus") %>%
               dplyr::summarise(P = mean(geneStatus == "Dynamic")) %>%
               dplyr::pull(P)
  n_dyn_genes <- ceiling(p_dynamic * n.genes.sample)
  n_norm_genes <- n.genes.sample - n_dyn_genes
  samp_dyn_genes <- SummarizedExperiment::rowData(sim.data) %>%
                    as.data.frame() %>%
                    dplyr::filter(geneStatus_overall == "Dynamic") %>%
                    dplyr::slice_sample(n = n_dyn_genes) %>%
                    rownames(.)
  samp_norm_genes <- SummarizedExperiment::rowData(sim.data) %>%
                     as.data.frame() %>%
                     dplyr::filter(geneStatus_overall == "NotDynamic") %>%
                     dplyr::slice_sample(n = n_norm_genes) %>%
                     rownames(.)
  samp_genes <- c(samp_dyn_genes, samp_norm_genes)
  sim.data <- sim.data[rownames(sim.data) %in% samp_genes, ]
  
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- as.matrix(t(SingleCellExperiment::counts(sim.data)))
  pt_df <- SummarizedExperiment::colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run scLANE
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        gene_stats <- testDynamic(expr.mat = sim_counts,
                                  pt = pt_df,
                                  n.potential.basis.fns = 5,
                                  parallel.exec = TRUE,
                                  n.cores = n.cores,
                                  approx.knot = TRUE, 
                                  track.time = TRUE,
                                  log.file = scLANE.log,
                                  log.iter = scLANE.log.iter)
        global_test_results <- getResultsDE(gene_stats, 
                                            p.adj.method = "holm", 
                                            fdr.cutoff = 0.01) %>%
                               dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = c("Gene" = "gene"))
        slope_test_results <- testSlope(test.dyn.results = gene_stats,
                                        p.adj.method = "holm",
                                        fdr.cutoff = 0.01) %>%
                              dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                 as.data.frame() %>%
                                                 dplyr::select(geneStatus_overall) %>% 
                                                 dplyr::mutate(gene = rownames(.))), 
                                                by = c("Gene" = "gene"))
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    td_res_raw_list[[i]] <- gene_stats
    td_res_tidy_list[[i]] <- global_test_results
    ts_res_list[[i]] <- slope_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(colnames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[, x], estimate = exp(y$Lineage_A$MARGE_Preds$marge_link_fit))
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- colnames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$testDynamic_results_raw <- td_res_raw_list
  res_list$testDynamic_results_tidy <- td_res_tidy_list
  res_list$testSlope_results <- ts_res_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}
