# run tradeSeq on subsampled single-subject data
run_tradeSeq <- function(sim.data = NULL,
                         n.genes.sample = 1000,
                         n.iter = 1,
                         param.list = NULL,
                         n.cores = 4) {
  # check inputs
  if (is.null(sim.data) || is.null(param.list)) { stop("You failed to provide necessary parameters to run_tradeSeq().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "tradeSeq_results_raw", "tradeSeq_results_tidy", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  ts_res_raw_list <- vector("list", length = n.iter)
  ts_res_tidy_list <- vector("list", length = n.iter)
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
  sim.data <- sim.data[, colSums(SingleCellExperiment::counts(sim.data)) > 0]
  sim_counts <- SingleCellExperiment::counts(sim.data)
  pt_df <- SummarizedExperiment::colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run tradeSeq
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        doParallel::registerDoParallel(n.cores)
        BiocParallel::register(BiocParallel::DoparParam())
        gene_stats <- tradeSeq::fitGAM(counts = sim_counts, 
                                       pseudotime = pt_df, 
                                       cellWeights = matrix(rep(1, nrow(pt_df)), ncol = 1), 
                                       nknots = 6, 
                                       sce = FALSE, 
                                       parallel = TRUE, 
                                       verbose = TRUE, 
                                       BPPARAM = BiocParallel::DoparParam())
        global_test_results <- tradeSeq::associationTest(models = gene_stats, global = TRUE) %>% 
                               dplyr::arrange(pvalue) %>% 
                               dplyr::mutate(gene = rownames(.), 
                                             pvalue_adj = stats::p.adjust(pvalue, method = "holm"), 
                                             gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>% 
                               dplyr::relocate(gene) %>% 
                               dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = "gene")
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    ts_res_raw_list[[i]] <- gene_stats
    ts_res_tidy_list[[i]] <- global_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(rownames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[x, ], estimate = y$fitted.values)
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- rownames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$tradeSeq_results_raw <- ts_res_raw_list
  res_list$tradeSeq_results_tidy <- ts_res_tidy_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}

# run tradeSeq on subsampled multi-subject data 
run_tradeSeq_multi <- function(sim.data = NULL,
                               n.genes.sample = 1000,
                               n.iter = 1, 
                               param.list = NULL,
                               n.cores = 4) {
  # check inputs
  if (is.null(sim.data) || is.null(param.list)) { stop("You failed to provide necessary parameters to run_tradeSeq_multi().") }
  if (n.iter <= 0) { stop("n.iter HAS to be positive, come on.") }
  if (n.cores <= 0) { stop("n.cores HAS to be positive, come on.") }
  if (n.genes.sample <= 0) { stop("n.genes.sample HAS to be positive, come on.") }
  
  # prepare results objects & sub-lists
  res_list <- vector("list", length = 8)
  names(res_list) <- c("sim_parameters", "start_time", "end_time", "time_diff", "mem_usage", "tradeSeq_results_raw", "tradeSeq_results_tidy", "RMSE_estimates")
  param_list <- vector("list", length = n.iter)
  start_time_list <- vector("list", length = n.iter)
  end_time_list <- vector("list", length = n.iter)
  time_diff_list <- vector("list", length = n.iter)
  mem_usage_list <- vector("list", length = n.iter)
  ts_res_raw_list <- vector("list", length = n.iter)
  ts_res_tidy_list <- vector("list", length = n.iter)
  rmse_list <- vector("list", length = n.iter)
  
  # prepare subsampled (preserves % dynamic genes) counts matrix & cell-ordering dataframe -- {targets} takes care of seed
  p_dynamic <- SummarizedExperiment::rowData(sim.data) %>%
               as.data.frame() %>% 
               dplyr::select(dplyr::contains("geneStatus_P")) %>%
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
  sim.data <- sim.data[, colSums(SingleCellExperiment::counts(sim.data)) > 0]
  
  # prepare counts matrix & cell-ordering dataframe
  sim_counts <- SingleCellExperiment::counts(sim.data)
  pt_df <- SummarizedExperiment::colData(sim.data) %>%
           as.data.frame() %>%
           dplyr::select(cell_time_normed) %>%
           dplyr::rename(PT = cell_time_normed)
  
  # run tradeSeq
  for (i in seq(n.iter)) {
    mem_usage_list[[i]] <- pryr::mem_change(
      {
        start_time <- Sys.time()
        doParallel::registerDoParallel(n.cores)
        BiocParallel::register(BiocParallel::DoparParam())
        gene_stats <- tradeSeq::fitGAM(counts = sim_counts, 
                                       pseudotime = pt_df, 
                                       cellWeights = matrix(rep(1, nrow(pt_df)), ncol = 1), 
                                       nknots = 6, 
                                       sce = FALSE, 
                                       parallel = TRUE, 
                                       verbose = TRUE, 
                                       BPPARAM = BiocParallel::DoparParam())
        global_test_results <- tradeSeq::associationTest(models = gene_stats, global = TRUE) %>% 
                               dplyr::arrange(pvalue) %>% 
                               dplyr::mutate(gene = rownames(.), 
                                             pvalue_adj = stats::p.adjust(pvalue, method = "bonferroni"), 
                                             gene_dynamic_overall = dplyr::case_when(pvalue_adj < 0.01 ~ 1, TRUE ~ 0)) %>% 
                               dplyr::relocate(gene) %>% 
                               dplyr::inner_join((SummarizedExperiment::rowData(sim.data) %>%
                                                  as.data.frame() %>%
                                                  dplyr::select(geneStatus_overall) %>% 
                                                  dplyr::mutate(gene = rownames(.))), 
                                                 by = "gene")
        end_time <- Sys.time()
        time_diff <- end_time - start_time
      }
    )
    start_time_list[[i]] <- start_time
    end_time_list[[i]] <- end_time
    time_diff_list[[i]] <- time_diff
    ts_res_raw_list[[i]] <- gene_stats
    ts_res_tidy_list[[i]] <- global_test_results
    param_list[[i]] <- param.list
    rmse_list[[i]] <- purrr::map2(rownames(sim_counts), 
                                  gene_stats, 
                                  function(x, y) {
                                    rmse <- try({ 
                                      yardstick::rmse_vec(truth = sim_counts[x, ], estimate = y$fitted.values)
                                    }, silent = TRUE) 
                                  })
    names(rmse_list[[i]]) <- rownames(sim_counts)
  }
  
  # set up results list & return
  res_list$sim_parameters <- param_list
  res_list$start_time <- start_time_list
  res_list$end_time <- end_time_list
  res_list$time_diff <- time_diff_list
  res_list$mem_usage <- mem_usage_list
  res_list$tradeSeq_results_raw <- ts_res_raw_list
  res_list$tradeSeq_results_tidy <- ts_res_tidy_list
  res_list$RMSE_estimates <- rmse_list
  return(res_list)
}
